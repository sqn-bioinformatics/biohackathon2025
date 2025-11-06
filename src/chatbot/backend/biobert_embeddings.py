# ingest_epmc_to_chroma.py
import json, re, textwrap, uuid
from pathlib import Path
import chromadb
from chromadb.config import Settings

from transformers import AutoTokenizer, AutoModel
import torch

EMBED_MODEL = "dmis-lab/biobert-base-cased-v1.1"
tok = AutoTokenizer.from_pretrained(EMBED_MODEL)
enc = AutoModel.from_pretrained(EMBED_MODEL).eval()

def embed(texts: list[str]) -> list[list[float]]:
    with torch.no_grad():
        inp = tok(texts, padding=True, truncation=True, return_tensors="pt")
        out = enc(**inp).last_hidden_state
        mask = inp["attention_mask"].unsqueeze(-1)
        pooled = (out * mask).sum(1) / mask.sum(1).clamp(min=1)
        return pooled.cpu().tolist()

def clean(txt: str) -> str:
    return re.sub(r"\s+", " ", (txt or "").strip())

def chunk(text: str, max_tokens=350, overlap=60):
    # naive token-ish split; replace with tiktoken if desired
    words = text.split()
    step = max_tokens - overlap
    for i in range(0, len(words), step):
        yield " ".join(words[i:i+max_tokens])

client = chromadb.Client(Settings(persist_directory="./chroma", anonymized_telemetry=False))
col = client.get_or_create_collection("epmc_biomed", embedding_function=None)

data = json.loads(Path("PMC6952411_metadata.json").read_text())
title = clean(data.get("title"))
abstract = clean(data.get("abstractText"))
body = clean(data.get("body", {}).get("text", ""))  # <-- full text present in your file
year = int(data.get("journalInfo", {}).get("yearOfPublication") or 0)
journal = (data.get("journalInfo", {}).get("journal", {}) or {}).get("title")
pmid = data.get("pmid"); pmcid = data.get("pmcid"); doi = data.get("doi")
mesh = [m["descriptorName"] for m in data.get("meshHeadingList", {}).get("meshHeading", []) if "descriptorName" in m]
keywords = data.get("keywordList", {}).get("keyword", [])
chems = [c["name"] for c in data.get("chemicalList", {}).get("chemical", [])]

base_meta = {
    "pmid": pmid, "pmcid": pmcid, "doi": doi, "year": year, "journal": journal,
    "mesh_headings": mesh, "keywords": keywords, "chemicals": chems,
    "title": title,
}

# Section-aware: keep title+abstract as their own small chunks
records = []
if title:
    records.append(("title", title))
if abstract:
    records.append(("abstract", abstract))

# Body chunks
for ch in chunk(body, max_tokens=350, overlap=60):
    # try to detect section tags if present in text; fallback "body"
    sect = "introduction" if "introduction" in ch.lower() else "body"
    records.append((sect, ch))

ids, docs, metas = [], [], []
for sect, txt in records:
    cid = str(uuid.uuid4())
    ids.append(cid)
    docs.append(txt)
    m = base_meta.copy()
    m["section"] = sect
    metas.append(m)

embs = embed(docs)
col.upsert(ids=ids, documents=docs, embeddings=embs, metadatas=metas)
client.persist()
print(f"Indexed {len(ids)} chunks for {pmcid}")
