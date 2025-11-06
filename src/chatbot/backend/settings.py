from pydantic import BaseModel
import os


class Settings(BaseModel):
    embed_model: str = os.getenv("EMBED_MODEL", "dmis-lab/biobert-base-cased-v1.1")
    gen_model: str = os.getenv("GEN_MODEL", "mistralai/Mistral-7B-Instruct-v0.2")
    chroma_dir: str = os.getenv("CHROMA_DIR", "./chroma")
    chroma_collection: str = os.getenv("CHROMA_COLLECTION", "epmc_biomed")
    top_k: int = int(os.getenv("TOP_K", 6))
    max_new_tokens: int = int(os.getenv("MAX_NEW_TOKENS", 512))
    temperature: float = float(os.getenv("TEMPERATURE", 0.2))
    top_p: float = float(os.getenv("TOP_P", 0.9))
    cors_origins: list[str] = os.getenv("CORS_ORIGINS", "http://localhost:5173").split(",")