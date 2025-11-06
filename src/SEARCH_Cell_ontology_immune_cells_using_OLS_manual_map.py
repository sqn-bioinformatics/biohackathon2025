"""
Script to map cell-type names from schmiedel, monac, mann2017 and blood atlas datasets to CL ontology IDs
 Start with the best match from the OLS in a manual_map (used ChatGPT with original terms as input and checked them to fix incorrect mappings using:
 Ontology Lookup Service https://www.ebi.ac.uk/ols4/
 Then the script  uses the API to get synonyms (could also get sub-classes and parents for greater or less specificity)
 Note: currently as well as exact synonyms for some cell types also returns broader synonyms for example neuttophil and basophil also map to polymorphonuclear leucocyte

"""

import requests
import json
import time


manual_map = {
    "Basophil": "CL:0000767",                     # basophil
    "CD4_CentralMemory": "CL:0000904",            # central memory CD4-positive, alpha-beta T cell
    "CD4_EffectorMemory": "CL:0000905",           # effector memory CD4-positive, alpha-beta T cell
    "CD4_TEMRA": "CL:0001087",                    # terminally differentiated effector CD4-positive, alpha-beta T cell
    "CD4_Naive": "CL:0000895",                    # naive CD4-positive, alpha-beta T cell
    "CD8_CentralMemory": "CL:0000907",            # central memory CD8-positive, alpha-beta T cell
    "CD8_EffectorMemory": "CL:0000913",           # effector and central memory CD8-positive, alpha-beta T cell are not distinguished in Cell Ontology?
    "CD8_TEMRA": "CL:0001062",                    # terminally differentiated effector CD8-positive, alpha-beta T cell
    "CD8_Naive": "CL:0000900",                    # naive CD8-positive, alpha-beta T cell
    "Eosinophil": "CL:0000771",                   # eosinophil
    "T_GammaDelta": "CL:0000798",                 # gamma-delta T cell
    "T_MAIT": "CL:0000940",                       # mucosal associated invariant T cell
    "Monocyte_Nonclassical": "CL:0000875",        # non-classical monocyte
    "Monocyte_Classical": "CL:0000860",           # classical monocyte
    "Monocyte_Intermediate": "CL:0002393",        # intermediate monocyte
    "B_Naive": "CL:0000788",                      # naive B cell
    "Neutrophil": "CL:0000775",                   # neutrophil
    "NK_Dim": "CL:0000939",                       # CD56-dim natural killer cell
    "NK_Bright": "CL:0000938",                    # CD56-bright natural killer cell
    "Plasmacytoid_DC": "CL:0000784",              # plasmacytoid dendritic cell
    "B_Memory": "CL:0000787",                     # memory B cell
    "CD4_MemoryTfh": "CL:0002038",                # memory T follicular helper cell
    "CD4_MemoryTh1": "CL:0000545",                # memory Th1 cell
    "CD4_MemoryTh2": "CL:0000546",                # memory Th2 cell
    "CD4_MemoryTh17": "CL:0000899",               # memory Th17 cell
    "CD4_Treg": "CL:0000792",                     # CD4+ regulatory T cell - also CD25-positive -IS THIS THE CORRECT BEST MATCH MAPPING?
    "Myeloid_DC": "CL:0001057"                    # myeloid dendritic cell
}

results = {}
DELAY = 0.2  # 0.2 sec  delay between API calls

for name, cl_id in manual_map.items():
    url = f"https://www.ebi.ac.uk/ols/api/ontologies/cl/terms?obo_id={cl_id.replace(':', '_')}"
    print(f"Fetching {cl_id} for {name} ...")

    try:
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        data = resp.json()

        terms = data.get("_embedded", {}).get("terms", [])
        if not terms:
            print(f"No term data found for {cl_id} ({name})")
            results[name] = [{"name": "Not found", "cl_id": cl_id, "exact_match": False, "synonyms": []}]
            continue

        term = terms[0]
        label = term.get("label", "N/A")
        synonyms = term.get("synonyms", [])
        iri = term.get("iri", "")

        results[name] = [{
            "name": label,
            "cl_id": cl_id,
            "exact_match": True,
            "synonyms": synonyms,
            "iri": iri
        }]

        print(f"  {name} → {label} ({len(synonyms)} synonyms)")

    except Exception as e:
        print(f"Error fetching {cl_id} ({name}): {e}")
        results[name] = [{"name": "Not found", "cl_id": cl_id, "exact_match": False, "synonyms": []}]
    time.sleep(DELAY)

# Save JSON
output_file = "cell_ontology_immune_cells_with_synonyms.json"
with open(output_file, "w", encoding="utf-8") as fh:
    json.dump(results, fh, indent=2, ensure_ascii=False)

print(f"\n  Results saved to {output_file}")

