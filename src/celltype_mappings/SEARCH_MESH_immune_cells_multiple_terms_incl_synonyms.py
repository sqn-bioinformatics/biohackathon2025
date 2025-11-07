import requests
import json
import time
import re


def normalize_label(s):
    # Lowercase, remove punctuation except hyphen and spaces, collapse spaces
    s = s.lower()
    s = re.sub(r"[^\w\s\-]", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s

def esearch(term, field_tag=None):
    q = f"{term}{field_tag}" if field_tag else term
    resp = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        params={"db": "mesh", "term": q, "retmode": "json"},
        timeout=15
    )
    resp.raise_for_status()
    return resp.json().get("esearchresult", {}).get("idlist", [])

def esummary(mesh_id):
    #Get esummary JSON for a MeSH numeric UID
    resp = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        params={"db": "mesh", "id": mesh_id, "retmode": "json"},
        timeout=15
    )
    resp.raise_for_status()
    return resp.json().get("result", {}).get(str(mesh_id), {})

cell_types = [
    "MAIT T-cell",
    "NK-cell",
    "regulatory T-cell", #changed from Treg in HPA
    "basophil",
    "classical monocyte",
    "eosinophil",
    "gamma-delta T-Cell Receptor",# changed from gdT-cell in HPA
    "intermediate monocyte",
    "memory B-cell",
    "memory CD4 T-cell",
    "memory CD8 T-cell",
    #"myeloid DC", #no good MESH terms for myeloid
    "naive B-cell",
    "naive CD4 T-cell",
    "naive CD8 T-cell",
    "neutrophil",
    "non-classical monocyte",
    "dendritic-cell" #changed from plasmacytoid DC to dendritic-cell, MESH term for Dendritic cell and plasmacytoid is an entry term (symonoym)
]

manual_map = {
    "memory CD4 T-cell": "68015496", # CD4-Positive T-Lymphocytes
    "naive CD4 T-cell": "68015496", # CD4-Positive T-Lymphocytes
    "MAIT T-cell": "2016683", # Mucosal-Associated Invariant T Cells
    "memory CD8 T-cell": "68018414", # CD8-Positive T-Lymphocytes
    "naive CD8 T-cell": "68018414" # CD8-Positive T-Lymphocytes
    #"plasmacytoid DC": "" # Dendritic Cells - check uid and add
}

results = {}
DELAY = 0.05  # seconds between requests

for term in cell_types:
    term_norm = normalize_label(term)
    mesh_hits = []

    # Search for exact whole-term descriptor search first
    try:
        exact_ids = esearch(term, field_tag="[MeSH Terms]")
    except Exception as e:
        print(f"esearch failed for exact term '{term}': {e}")
        exact_ids = []
    time.sleep(DELAY)

    if exact_ids:
        for mid in exact_ids:
            try:
                summary = esummary(mid)
            except Exception as e:
                print(f"esummary failed for id {mid}: {e}")
                continue

            name = summary.get("name") or (summary.get("ds_meshterms") or [None])[0]
            if not name:
                continue

            ds_terms = summary.get("ds_meshterms") or []
            # exact_match if canonical name or any synonym matches full term
            exact_match = term_norm == normalize_label(name) or any(
                normalize_label(t) == term_norm for t in ds_terms
            )

            mesh_hits.append({
                "name": name,
                "numeric_uid": str(mid),
                "exact_match": exact_match,
                "synonyms": ds_terms
            })
            time.sleep(DELAY)

        # prefer exact matches if any
        exact_matches = [h for h in mesh_hits if h.get("exact_match")]
        if exact_matches:
            mesh_hits = exact_matches
    else:
        # fallback to splitting into words
        words = term.split()
        collected = []
        for word in words:
            try:
                ids = esearch(word, field_tag="[MeSH Terms]")
            except Exception as e:
                print(f"esearch failed for word '{word}': {e}")
                ids = []
            time.sleep(DELAY)

            for mid in ids:
                try:
                    summary = esummary(mid)
                except Exception as e:
                    print(f"esummary failed for id {mid}: {e}")
                    continue

                name = summary.get("name") or (summary.get("ds_meshterms") or [None])[0]
                if not name:
                    continue

                ds_terms = summary.get("ds_meshterms") or []

                collected.append({
                    "name": name,
                    "numeric_uid": str(mid),
                    "exact_match": False,
                    "synonyms": ds_terms
                })
                time.sleep(DELAY)

        mesh_hits = collected

    # Insert manual map at the front if exists
    if term in manual_map:
        mid = manual_map[term]
        try:
            summary = esummary(mid)
            name = summary.get("name") or (summary.get("ds_meshterms") or [None])[0]
            ds_terms = summary.get("ds_meshterms") or []
            mesh_hits.insert(0, {  # insert at front to ensure priority
                "name": name,
                "numeric_uid": str(mid),
                "exact_match": True,
                "synonyms": ds_terms
            })
        except Exception as e:
            print(f"esummary failed for manual_map id {mid} ({term}): {e}")
        time.sleep(DELAY)

    # not found any hits
    if not mesh_hits:
        results[term] = [{"name": "Not found", "numeric_uid": "N/A", "exact_match": False, "synonyms": []}]
    else:
        # deduplicate by (name, numeric_uid)
        seen = set()
        unique = []
        for h in mesh_hits:
            key = (h["name"], h["numeric_uid"])
            if key not in seen:
                seen.add(key)
                unique.append(h)
        results[term] = unique

# Save JSON
with open("mesh_immune_cells_with_synonyms.json", "w", encoding="utf-8") as fh:
    json.dump(results, fh, indent=2, ensure_ascii=False)

print(json.dumps(results, indent=2, ensure_ascii=False))

