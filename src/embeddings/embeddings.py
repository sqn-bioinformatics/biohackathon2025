import math
from typing import Dict, Iterable, List, Tuple

import numpy as np
import torch
from transformers import AutoModel, AutoTokenizer


class Embedder:
    def __init__(self, model_name: str, device: str | torch.device = None):
        try:
            # Try local files, and if not available, retry with download
            self.model = AutoModel.from_pretrained(model_name, local_files_only=True)
            self.tokenizer = AutoTokenizer.from_pretrained(
                model_name, use_fast=True, local_file_only=True
            )
        except Exception:
            self.model = AutoModel.from_pretrained(model_name, local_files_only=False)
            self.tokenizer = AutoTokenizer.from_pretrained(
                model_name, use_fast=True, local_file_only=False
            )

        self.device = (
            torch.device(device)
            if device
            else torch.device("cuda" if torch.cuda.is_available() else "cpu")
        )
        print("Running embedder on device:", self.device)
        self.model.to(self.device).eval()

        if hasattr(self.model.config, "max_position_embeddings"):
            self.max_len_default = int(self.model.config.max_position_embeddings)
        else:
            self.max_len_default = 512

        self.max_len_default = min(self.max_len_default, 512)

    def _chunk_with_overlap(
        self,
        text: str,
        max_tokens: int = None,
        overlap: int = 128,
        add_special_tokens: bool = True,
        return_offsets: bool = True,
    ) -> Tuple[List[Dict[str, List[int]]], List[Tuple[int, int]]]:
        """
        Splits `text` into overlapping token windows.
        Returns:
          - a list of dicts like {"input_ids": [...], "attention_mask": [...]}
          - a parallel list of (char_start, char_end) spans for each chunk (if return_offsets)
        """
        max_tokens = max_tokens or self.max_len_default
        assert overlap < max_tokens, "overlap must be smaller than max_tokens"

        enc = self.tokenizer(
            text,
            max_length=max_tokens,
            stride=overlap,
            truncation=True,
            padding=True,
            return_overflowing_tokens=True,
            return_offsets_mapping=return_offsets,
            add_special_tokens=add_special_tokens,
        )

        # Each overflowed window is one "encoding"
        input_id_windows = enc["input_ids"]
        attn_windows = enc["attention_mask"]

        spans = []
        if return_offsets:
            # For each window, offsets map each token to (char_start,char_end) in the original text
            # We take the min start and max end among non-special tokens
            for offsets in enc["offset_mapping"]:
                # Filter out special tokens which often carry (0,0) or (0,0)-like placeholders
                valid = [(s, e) for (s, e) in offsets if (e > s)]
                if valid:
                    spans.append((valid[0][0], valid[-1][1]))
                else:
                    spans.append((0, 0))

        windows = []
        for ids, mask in zip(input_id_windows, attn_windows):
            windows.append({"input_ids": ids, "attention_mask": mask})

        return windows, spans

    @torch.no_grad()
    def _embed_batch(self, batch: Dict[str, List[List[int]]]) -> np.ndarray:
        """Mean-pools token embeddings for each sequence in the batch."""
        input_ids = torch.tensor(
            batch["input_ids"], dtype=torch.long, device=self.device
        )
        attention_mask = torch.tensor(
            batch["attention_mask"], dtype=torch.long, device=self.device
        )

        outputs = self.model(input_ids=input_ids, attention_mask=attention_mask)
        token_embeddings = outputs.last_hidden_state  # [B, T, H]

        # Mean-pool over non-padded tokens
        mask = attention_mask.unsqueeze(-1)  # [B, T, 1]
        summed = (token_embeddings * mask).sum(dim=1)
        counts = mask.sum(dim=1).clamp(min=1)
        sentence_embeddings = summed / counts
        return sentence_embeddings.detach().cpu().numpy()

    def embed_text(
        self,
        text: str,
        max_tokens: int = None,
        overlap: int = 128,
        batch_size: int = 8,
        aggregate: str = "none",  # "none" | "mean" | "max"
    ) -> tuple[list[str], np.ndarray]:
        """
        Embeds long text by sliding-window tokenization with overlap.
        Returns:
          - if aggregate == "none": array of shape [num_chunks, hidden_size]
          - else: single vector [hidden_size]
        """
        windows, spans = self._chunk_with_overlap(
            text, max_tokens=max_tokens, overlap=overlap
        )

        # Simple batching
        all_vecs = []
        segment_strings = []
        for i in range(0, len(windows), batch_size):
            batch = {"input_ids": [], "attention_mask": []}
            batch_windows = windows[i : i + batch_size]
            batch_spans = spans[i : i + batch_size]

            for w, (start, end) in zip(batch_windows, batch_spans):
                # Extract original text chunk using character offsets
                segment_strings.append(text[start:end])
                batch["input_ids"].append(w["input_ids"])
                batch["attention_mask"].append(w["attention_mask"])
            vecs = self._embed_batch(batch)
            all_vecs.append(vecs)

        chunk_embeddings = (
            np.vstack(all_vecs)
            if all_vecs
            else np.zeros((0, self.model.config.hidden_size))
        )

        if aggregate == "none":
            return segment_strings, chunk_embeddings

        if aggregate == "mean":
            return segment_strings, chunk_embeddings.mean(axis=0)
        elif aggregate == "max":
            return segment_strings, chunk_embeddings.max(axis=0)
        else:
            raise ValueError("aggregate must be one of {'none','mean','max'}")
