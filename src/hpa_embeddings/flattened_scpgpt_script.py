import numpy as np
import scanpy as sc
import sys
import os
import warnings
import torch
from pathlib import Path
import json
import pandas as pd
from scipy.stats import mode

import torch

from hpa_embeddings.load_datasets import load_datasets
from hpa_embeddings.utils import get_batch_cell_embeddings

print(torch.__version__)


# extra dependency for similarity search
try:
    import hnswlib

    hnswlib_imported = True
except ImportError:
    hnswlib_imported = False
    print(
        "hnswlib not installed! We highly recommend installing it for fast similarity search."
    )
    print("To install it, run: pip install hnswlib")

sys.path.insert(0, "../")


from scgpt.loss import masked_mse_loss, masked_relative_error
from scgpt.model import TransformerModel
from scgpt.preprocess import Preprocessor
from scgpt.tokenizer import tokenize_and_pad_batch, random_mask_value
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt import logger
from scgpt import SubsetsBatchSampler
import sklearn


adata_dict, gene_col = load_datasets()



n_hvg = max_seq_len = 1200
DO_HVG = True
HAS_BATCH = batch_key = None

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")



# Load model vocabulary and configuration
model_dir = Path("./model/")  # Update this path to your model directory
vocab_file = model_dir / "vocab.json"
model_config_file = model_dir / "args.json"
model_file = model_dir / "best_model.pt"
pad_token = "<pad>"
special_tokens = [pad_token, "<cls>", "<eoc>"]

# Load vocabulary
vocab = GeneVocab.from_file(vocab_file)
for s in special_tokens:
    if s not in vocab:
        vocab.append_token(s)

vocab.set_default_index(vocab["<pad>"])

# Load model configuration
with open(model_config_file, "r") as f:
    model_configs = json.load(f)



for dset in adata_dict:
    print(dset)
    adata = adata_dict[dset]
    # Match genes to vocabulary
    adata.var["id_in_vocab"] = [
        vocab[gene] if gene in vocab else -1 for gene in adata.var[gene_col]
    ]
    gene_ids_in_vocab = np.array(adata.var["id_in_vocab"])
    logger.info(
        f"match {np.sum(gene_ids_in_vocab >= 0)}/{len(gene_ids_in_vocab)} genes "
        f"in vocabulary of size {len(vocab)}."
    )
    adata = adata[:, adata.var["id_in_vocab"] >= 0] # what does it mean?

    adata_dict[dset] = adata
    print()



model = TransformerModel(
    ntoken=len(vocab),
    d_model=model_configs["embsize"],
    nhead=model_configs["nheads"],
    d_hid=model_configs["d_hid"],
    nlayers=model_configs["nlayers"],
    nlayers_cls=model_configs["n_layers_cls"],
    n_cls=1,
    vocab=vocab,
    dropout=model_configs["dropout"],
    pad_token=model_configs["pad_token"],
    pad_value=model_configs["pad_value"],
    do_mvc=True,
    do_dab=False,
    use_batch_labels=False,
    # num_batch_labels=num_batch_types,
    domain_spec_batchnorm=False,
    explicit_zero_prob=False,
    use_fast_transformer=True,
    fast_transformer_backend="flash",
    pre_norm=False,
)

# Load the checkpoint
checkpoint = torch.load(model_file, map_location=device)

# Handle naming differences between PyTorch 1.x and 2.x Transformer layers
state_dict = checkpoint

# Map "Wqkv" -> "in_proj_weight"/"in_proj_bias" if needed
mapped_state_dict = {}
for k, v in state_dict.items():
    if "Wqkv.weight" in k:
        base = k.replace("Wqkv.weight", "")
        # Split Wqkv.weight into in_proj_weight format if possible
        mapped_state_dict[base + "in_proj_weight"] = v
    elif "Wqkv.bias" in k:
        base = k.replace("Wqkv.bias", "")
        mapped_state_dict[base + "in_proj_bias"] = v
    else:
        mapped_state_dict[k] = v

# Load model parameters, ignoring missing/unexpected keys safely
model_dict = model.state_dict()
pretrained_dict = {
    k: v for k, v in mapped_state_dict.items()
    if k in model_dict and v.shape == model_dict[k].shape
}
model_dict.update(pretrained_dict)
model.load_state_dict(model_dict, strict=False)

model.to(device)
model.eval()

embeddings_dict = {}

for dset in adata_dict:
    print(dset)
    adata = adata_dict[dset]

    data_is_raw = adata.X.max() > 30
    print("data is raw: ", data_is_raw)
    preprocessor = Preprocessor(
        use_key="X",  # the key in adata.layers to use as raw data
        filter_gene_by_counts=False,  # step 1
        filter_cell_by_counts=False,  # step 2
        normalize_total=1e4,  # 3. whether to normalize the raw data and to what sum
        result_normed_key="X_normed",  # the key in adata.layers to store the normalized data
        log1p=data_is_raw,  # 4. whether to log1p the normalized data
        result_log1p_key="X_log1p",
        subset_hvg=n_hvg if DO_HVG else False,
        hvg_flavor="seurat_v3" if data_is_raw else "cell_ranger",
        binning=51,  # 6. whether to bin the raw data and to what number of bins
        result_binned_key="counts",  # the key in adata.layers to store the binned data
    )
    preprocessor(adata, batch_key=batch_key if HAS_BATCH else None)

    # Get gene IDs for the preprocessed genes
    genes = adata.var[gene_col].tolist()
    gene_ids = np.array(vocab(genes), dtype=int)
    print(len(genes), len(gene_ids))

    if HAS_BATCH:  # None
        adata.obs["batch_id"] = adata.obs[batch_key].astype("category").cat.codes.values
    else:
        adata.obs["batch_id"] = adata.obs[cell_type_key].astype("category").cat.codes.values
    batch_ids = adata.obs["batch_id"].tolist()
    num_batch_types = len(set(batch_ids))
    batch_ids = np.array(batch_ids)

    # Tokenize input # What is it?
    if isinstance(adata.layers["counts"], np.ndarray):
        all_counts = adata.layers["counts"]
    else:
        all_counts = adata.layers["counts"].A

    num_of_non_zero_genes = [np.count_nonzero(all_counts[i]) for i in range(all_counts.shape[0])]

    max_length = np.max(num_of_non_zero_genes) + 1  # plus 1 for appending <cls>
    print(max_length, max_seq_len + 1)
    max_length = min(max_length, max_seq_len + 1)
    print(max_length, "\n")

    ref_cell_embeddings = get_batch_cell_embeddings(
        adata,
        cell_embedding_mode="cls",
        model=model,
        vocab=vocab,
        max_length=max_length,
        model_configs=model_configs,
        gene_ids=gene_ids,
        use_batch_labels=False,
    )

    embeddings_dict[dset] = ref_cell_embeddings

    print()



all_embeddings = [embeddings_dict[k] for k in embeddings_dict]
all_embeddings = np.concatenate(all_embeddings, axis=0)
all_embeddings.shape



import anndata as ad



combined = ad.concat(adata_dict.values(), axis=0, join="outer")
# combined



# Optional step to visualize the reference dataset using the embeddings
ref_embed_adata = sc.AnnData(all_embeddings, obs=pd.DataFrame(combined.obs))

sc.pp.neighbors(ref_embed_adata, use_rep="X")
sc.tl.umap(ref_embed_adata)

sc.pl.umap(ref_embed_adata, color=cell_type_key, frameon=False, wspace=0.4, size=100)

for dset in adata_dict:
    adata = adata_dict[dset]
    embeddings = embeddings_dict[dset]

    # Optional step to visualize the reference dataset using the embeddings
    ref_embed_adata = sc.AnnData(embeddings, obs=pd.DataFrame(adata.obs))

    sc.pp.neighbors(ref_embed_adata, use_rep="X")
    sc.tl.umap(ref_embed_adata)

    sc.pl.umap(ref_embed_adata, color=cell_type_key, frameon=False, wspace=0.4, size=100)



# Optional step to visualize the reference dataset using the embeddings
ref_embed_adata = sc.AnnData(ref_cell_embeddings, obs=pd.DataFrame(adata.obs))
ref_embed_adata



# Optional step to visualize the reference dataset using the embeddings
ref_embed_adata = sc.AnnData(ref_cell_embeddings, obs=pd.DataFrame(adata.obs))

sc.pp.neighbors(ref_embed_adata, use_rep="X")
sc.tl.umap(ref_embed_adata)



sc.pl.umap(ref_embed_adata, color=cell_type_key, frameon=False, wspace=0.4, size=100)



sc.pl.umap(ref_embed_adata, color=cell_type_key, frameon=False, wspace=0.4, size=100)



sc.pl.umap(ref_embed_adata, color="lineage", frameon=False, wspace=0.4, size=100)