# app.py
from fastapi import FastAPI
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
import chromadb
from chromadb.config import Settings
from transformers import AutoTokenizer, AutoModel, AutoModelForCausalLM, TextIteratorStreamer
import torch, threading

EMBED_MODEL = "dmis-lab/biobert-base-cased-v1.1"
tokE = AutoTokenizer.from_pretrained(EMBED_MODEL)
enc = AutoModel.from_pretrained(EMBED_MODEL).eval()

def embed_one(text: str):
    with torch.no_grad():
        inp = tokE([text], padding=True, truncation=True, return_tensors="pt")
        out = enc(**inp).last_hidden_state
        mask = inp["attention_mask"].unsqueeze(-1)
        pooled = (out * mask).sum(1) / mask.sum(1).clamp(min=1)
        return pooled[0].cpu().tolist()

GEN_MODEL = "mistralai/Mistral-7B-Instruct-v0.2"  # swap for your biomedical generator
tokG = AutoTokenizer.from_pretrained(GEN_MODEL)
gen = AutoModelForCausalLM.from_pretrained(GEN_MODEL, torch_dtype=torch.float16, device_map="auto")

client = chromadb.Client(Settings(persist_directory="./chroma", anonymized_telemetry=False))
col = client.get_or_create_collection("epmc_biomed", embedding_function=None)

SYSTEM = "You are a concise biomedical assistant. Cite PMIDs/PMCIDs when relevant. Do not give medical advice."

class ChatReq(BaseModel):
    message: str
    top_k: int = 6
    year_min: int | None = None
    mesh_any: list[str] | None = None

@app.post("/chat")
def chat(req: ChatReq):
    q_emb = embed_one(req.message)
    where = {}
    if req.year_min: where["year"] = {"$gte": req.year_min}
    if req.mesh_any: where["mesh_headings"] = {"$in": req.mesh_any}

    res = col.query(query_embeddings=[q_emb], n_results=req.top_k, where=where)
    docs = res.get("documents", [[]])[0]
    metas = res.get("metadatas", [[]])[0]

    # Build a compact, citable context
    def cite(m):
        pm = m.get("pmid") or ""
        pc = m.get("pmcid") or ""
        j = m.get("journal") or ""
        y = m.get("year") or ""
        return f"[{j} {y}; PMID:{pm} PMCID:{pc}]"
    context = "\n\n".join(f"{d}\n{cite(m)}" for d, m in zip(docs, metas))

    prompt = (
        f"<s>[SYSTEM]\n{SYSTEM}\n"
        f"[CONTEXT]\n{context}\n"
        f"[/CONTEXT]\n[USER]\n{req.message}\n[/USER]\n[ASSISTANT]"
    )

    inputs = tokG(prompt, return_tensors="pt").to(gen.device)
    streamer = TextIteratorStreamer(tokG, skip_special_tokens=True)
    kwargs = dict(**inputs, max_new_tokens=512, temperature=0.2, top_p=0.9, streamer=streamer)
    threading.Thread(target=gen.generate, kwargs=kwargs).start()

    def stream():
        for chunk in streamer:
            yield chunk
    return StreamingResponse(stream(), media_type="text/plain")
