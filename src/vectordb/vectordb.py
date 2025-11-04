# vectordb.py
from typing import Iterable, Optional

import chromadb
import numpy as np


class VectorDB:
    def __init__(self, db_path: str):
        # Persistent on-disk client (0.4+/0.5+ style)
        self.client = chromadb.PersistentClient(path=db_path)
        self.collection = self.client.get_or_create_collection(name="embeddings")

    @staticmethod
    def _ensure_2d(vectors: np.ndarray) -> np.ndarray:
        if vectors.ndim == 1:
            return vectors.reshape(1, -1)
        if vectors.ndim != 2:
            raise ValueError(
                f"Expected 2D array (segments, hidden_dim); got shape {vectors.shape}"
            )
        return vectors

    def add_article_vectors(self, pubmed_id: int, vectors: np.ndarray) -> None:
        """
        Store all segment embeddings for one article.

        vectors: np.ndarray of shape (segments, hidden_dim)
        """
        vectors = self._ensure_2d(np.asarray(vectors, dtype=np.float32))
        n_segments = vectors.shape[0]

        # Build stable string IDs; metadata carries both pubmed_id and segment number
        ids = [f"{pubmed_id}:{i}" for i in range(n_segments)]
        metadata = [
            {"pubmed_id": pubmed_id, "segment": int(i)} for i in range(n_segments)
        ]
        embeddings = vectors.tolist()

        # Use upsert so re-running doesn't raise on existing IDs
        self.collection.upsert(ids=ids, embeddings=embeddings, metadata=metadata)

    def get_article_vector(
        self, pubmed_id: int, segment_number: int
    ) -> Optional[np.ndarray]:
        """
        Retrieve a single segment embedding by (pubmed_id, segment_number).
        """
        # Prefer metadata filter to avoid having to reconstruct the ID outside
        res = self.collection.get(
            where={"pubmed_id": pubmed_id, "segment": int(segment_number)},
            include=["embeddings", "metadata"],
            limit=1,
        )
        embs = res.get("embeddings") or []
        if not embs:
            return None
        return np.array(embs[0], dtype=np.float32)

    def get_article_vectors(self, pubmed_id: int) -> Optional[np.ndarray]:
        """
        Retrieve all embeddings for an article, sorted by segment number.
        Returns an array of shape (segments, hidden_dim), or None if not found.
        """
        res = self.collection.get(
            where={"pubmed_id": pubmed_id},
            include=["embeddings", "metadata"],
        )
        embs = res.get("embeddings") or []
        metas = res.get("metadata") or []

        if not embs:
            return None

        # Sort by segment to ensure deterministic order
        order = np.argsort([m.get("segment", 0) for m in metas])
        embs_sorted: Iterable[list[float]] = [embs[i] for i in order]
        return np.asarray(list(embs_sorted), dtype=np.float32)

    def delete_article(self, pubmed_id: int) -> int:
        """
        Delete all segments for an article. Returns number of deleted items.
        """
        res = self.collection.get(where={"pubmed_id": pubmed_id}, include=["metadata"])
        ids = res.get("ids") or []
        if ids:
            self.collection.delete(ids=ids)
        return len(ids)
