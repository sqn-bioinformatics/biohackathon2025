# vectordb.py
import json
from pprint import pp
from typing import Iterable, Optional

import chromadb
import numpy as np
from embeddings import Embedder
from pydantic import BaseModel


class TextMetadata(BaseModel):
    pubmed_id: int
    pmc_id: str
    doi: str
    title: str
    authors: str
    year: int
    license: str
    mesh_terms: str
    abstract: str


class VectorDB:
    def __init__(self, db_path: str, embedder: Embedder):
        # Persistent on-disk client (0.4+/0.5+ style)
        self.client = chromadb.PersistentClient(path=db_path)
        self.collection = self.client.get_or_create_collection(name="embeddings")
        self.embedder = embedder

    @staticmethod
    def _ensure_2d(vectors: np.ndarray) -> np.ndarray:
        if vectors.ndim == 1:
            return vectors.reshape(1, -1)
        if vectors.ndim != 2:
            raise ValueError(
                f"Expected 2D array (segments, hidden_dim); got shape {vectors.shape}"
            )
        return vectors

    def query(self, text: str):
        embeddings = self.embedder.embed_text(text)[1].tolist()
        return self.collection.query(query_embeddings=embeddings)

    @staticmethod
    def _sanitize_metadata(md: dict) -> dict:
        """Coerce metadata so every value is a primitive allowed by Chroma."""
        out = {}
        for k, v in md.items():
            if isinstance(v, (str, int, float, bool)) or v is None:
                out[k] = v
            elif isinstance(v, (list, tuple, set)):
                # Empty collections become None; otherwise JSON-serialize
                out[k] = None if not v else json.dumps(list(v), ensure_ascii=False)
            elif isinstance(v, dict):
                out[k] = None if not v else json.dumps(v, ensure_ascii=False)
            else:
                # Fallback to string for e.g. enums, dataclasses, pydantic types, numpy scalars, etc.
                out[k] = str(v)
        return out

    def add_article_vectors(self, text: str, metadata: TextMetadata) -> None:
        """
        Store all segment embeddings for one article.

        vectors: np.ndarray of shape (segments, hidden_dim)
        """
        segments, vectors = self.embedder.embed_text(text=text)
        vectors = self._ensure_2d(np.asarray(vectors, dtype=np.float32))
        n_segments = vectors.shape[0]

        base_md = (
            metadata.model_dump()
        )  # self._sanitize_metadata(metadata.model_dump())

        # Build stable string IDs; metadata carries both pubmed_id and segment number
        ids = [f"{metadata.pubmed_id}:{i}" for i in range(n_segments)]
        metadata_ = [{**base_md, "segment": i} for i in range(n_segments)]
        embeddings = vectors.tolist()

        # Use upsert so re-running doesn't raise on existing IDs
        self.collection.upsert(
            ids=ids, embeddings=embeddings, metadatas=metadata_, documents=segments
        )

    def get_article_vector(
        self, pubmed_id: int, segment_number: int
    ) -> Optional[np.ndarray]:
        """
        Retrieve a single segment embedding by (pubmed_id, segment_number).
        """
        # Prefer metadata filter to avoid having to reconstruct the ID outside
        res = self.collection.get(
            where={
                "$and": [{"pubmed_id": pubmed_id}, {"segment": int(segment_number)}]
            },
            include=[
                chromadb.api.types.IncludeEnum.embeddings,
                chromadb.api.types.IncludeEnum.metadatas,
            ],
            limit=1,
        )
        embs = res.get("embeddings")

        if embs.shape[0] == 0:
            return None
        return embs[0]

    def get_article_vectors(self, pubmed_id: int) -> Optional[np.ndarray]:
        """
        Retrieve all embeddings for an article, sorted by segment number.
        Returns an array of shape (segments, hidden_dim), or None if not found.
        """
        res = self.collection.get(
            where={"pubmed_id": pubmed_id},
            include=[
                chromadb.api.types.IncludeEnum.embeddings,
                chromadb.api.types.IncludeEnum.metadatas,
            ],  # "embeddings", "metadata"],
        )
        # print(res.keys())
        embs = res.get("embeddings")
        metas = res.get("metadatas")

        if embs.shape[0] == 0 or metas is None:
            return None

        # Sort by segment to ensure deterministic order
        order = np.argsort([m.get("segment", 0) for m in metas])
        embs_sorted = embs[order]
        return embs_sorted

    def article_exists(self, pubmed_id: int) -> bool:
        key = f"{pubmed_id}:0"
        return bool(self.collection.get(ids=[key])["ids"])
