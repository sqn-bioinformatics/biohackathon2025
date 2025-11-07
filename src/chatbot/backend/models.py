from typing import List

# Optional dependency: torch. We avoid hard-failing on import so the app can start
# and provide a clear error only when GPU/torch functionality is actually used.
try:
    import torch  # type: ignore
    TORCH_AVAILABLE = True
except Exception:  # pragma: no cover - environment without torch
    torch = None  # type: ignore
    TORCH_AVAILABLE = False

# Optional dependency: transformers. Avoid hard-failing at import time so the app can start.
try:
    from transformers import AutoTokenizer, AutoModel, AutoModelForCausalLM, TextIteratorStreamer  # type: ignore
    TRANSFORMERS_AVAILABLE = True
except Exception:  # pragma: no cover - environment without transformers
    AutoTokenizer = AutoModel = AutoModelForCausalLM = TextIteratorStreamer = None  # type: ignore
    TRANSFORMERS_AVAILABLE = False


class Embedder:
    def __init__(self, model_name: str):
        if not TORCH_AVAILABLE:
            raise ImportError(
                "PyTorch is not installed but is required for Embedder. Install it with: pip install 'torch>=2.2.0' --index-url https://download.pytorch.org/whl/cpu (for CPU) or follow https://pytorch.org/get-started/locally/ for CUDA."
            )
        if not TRANSFORMERS_AVAILABLE:
            raise ImportError(
                "Transformers is not installed but is required for Embedder. Install it with: pip install 'transformers>=4.44.2'"
            )
        self.tok = AutoTokenizer.from_pretrained(model_name)
        self.mdl = AutoModel.from_pretrained(model_name).eval()


    def __call__(self, texts: List[str]) -> List[List[float]]:
        if not TORCH_AVAILABLE:
            raise ImportError("PyTorch is required to compute embeddings.")
        with torch.no_grad():
            inp = self.tok(texts, padding=True, truncation=True, return_tensors="pt")
            out = self.mdl(**inp).last_hidden_state
            mask = inp["attention_mask"].unsqueeze(-1)
            pooled = (out * mask).sum(1) / mask.sum(1).clamp(min=1)
            return pooled.cpu().tolist()


class Generator:
    def __init__(self, model_name: str, max_new_tokens: int, temperature: float, top_p: float):
        if not TORCH_AVAILABLE:
            raise ImportError(
                "PyTorch is not installed but is required for text generation. Install it with: pip install 'torch>=2.2.0' --index-url https://download.pytorch.org/whl/cpu (for CPU) or follow https://pytorch.org/get-started/locally/ for CUDA."
            )
        if not TRANSFORMERS_AVAILABLE:
            raise ImportError(
                "Transformers is not installed but is required for text generation. Install it with: pip install 'transformers>=4.44.2'"
            )
        self.tok = AutoTokenizer.from_pretrained(model_name)
        # Use float16 only if supported; otherwise fall back to default
        torch_dtype = torch.float16 if TORCH_AVAILABLE else None
        self.mdl = AutoModelForCausalLM.from_pretrained(
            model_name, torch_dtype=torch_dtype, device_map="auto"
        )
        self.max_new_tokens = max_new_tokens
        self.temperature = temperature
        self.top_p = top_p

    def stream(self, prompt: str):
        inputs = self.tok(prompt, return_tensors="pt").to(self.mdl.device)
        streamer = TextIteratorStreamer(self.tok, skip_special_tokens=True)
        kwargs = dict(
            **inputs,
            max_new_tokens=self.max_new_tokens,
            temperature=self.temperature,
            top_p=self.top_p,
            do_sample=True,
            streamer=streamer,
        )
        import threading

        threading.Thread(target=self.mdl.generate, kwargs=kwargs).start()
        for token in streamer:
            yield token
