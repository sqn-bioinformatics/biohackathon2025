from typing import List
import torch
from transformers import AutoTokenizer, AutoModel, AutoModelForCausalLM, TextIteratorStreamer


class Embedder:
    def __init__(self, model_name: str):
        self.tok = AutoTokenizer.from_pretrained(model_name)
        self.mdl = AutoModel.from_pretrained(model_name).eval()


    def __call__(self, texts: List[str]) -> List[List[float]]:
        with torch.no_grad():
            inp = self.tok(texts, padding=True, truncation=True, return_tensors="pt")
            out = self.mdl(**inp).last_hidden_state
            mask = inp["attention_mask"].unsqueeze(-1)
            pooled = (out * mask).sum(1) / mask.sum(1).clamp(min=1)
            return pooled.cpu().tolist()


class Generator:
    def __init__(self, model_name: str, max_new_tokens: int, temperature: float, top_p: float):
        self.tok = AutoTokenizer.from_pretrained(model_name)
        self.mdl = AutoModelForCausalLM.from_pretrained(
            model_name, torch_dtype=torch.float16, device_map="auto"
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
