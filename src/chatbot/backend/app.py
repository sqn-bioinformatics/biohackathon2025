import os
from typing import Optional, List
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
import chromadb
from chromadb.config import Settings as ChromaSettings

from settings import Settings
from models import Embedder, Generator

cfg = Settings()
app = FastAPI(title="RAG Chat API")
app.add_middleware(
    CORSMiddleware,
    allow_origins=cfg.cors_origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- Model + Chroma init ---
embedder = Embedder(cfg.embed_model)

# NOTE: swap Generator for your hosted endpoint if needed
generator = Generator(cfg.gen_model, cfg.max_new_tokens, cfg.temperature, cfg.top_p)

client = chromadb.Client(ChromaSettings(persist_directory=cfg.chroma_dir, anonymized_telemetry=False))
col = client.get_or_create_collection(cfg.chroma_collection, embedding_function=None)

SYSTEM = (
    "You are a concise biomedical assistant. Cite PMIDs/PMCIDs when present in the context. "
    "Do not provide medical advice."
)


class ChatReq(BaseModel):
    message: str
    top_k: Optional[int] = None
    year_min: Optional[int] = None
    mesh_any: Optional[List[str]] = None


@app.get("/health")
def health():
    return {"status": "ok"}


@app.post("/chat")
def chat(req: ChatReq):
    top_k = req.top_k or cfg.top_k

# 1) retrieve
    q_emb = embedder([req.message])[0]
    where = {}
    if req.year_min:
        where["year"] = {"$gte": req.year_min}
    if req.mesh_any:
        where["mesh_headings"] = {"$in": req.mesh_any}

    res = col.query(query_embeddings=[q_emb], n_results=top_k, where=where)
    docs = res.get("documents", [[]])[0]
    metas = res.get("metadatas", [[]])[0]


# 2) format context with citations
def cite(m):
    pmid = m.get("pmid") or ""
    pmcid = m.get("pmcid") or ""
    j = m.get("journal") or ""
    y = m.get("year") or ""
    s = m.get("section") or ""
    return f"[{j} {y} {s} PMID:{pmid} PMCID:{pmcid}]"

context = "\n\n".join(f"{d}\n{cite(m)}" for d, m in zip(docs, metas))

prompt = (
    f"<s>[SYSTEM]\n{SYSTEM}\n"
    f"[CONTEXT]\n{context or 'No context.'}\n[/CONTEXT]\n"
    f"[USER]\n{req.message}\n[/USER]\n[ASSISTANT]"
)


# 3) stream generation to client
def stream():

    for token in generator.stream(prompt):
        yield token

    return StreamingResponse(stream(), media_type="text/plain")
