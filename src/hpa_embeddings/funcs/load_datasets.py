import sys
from importlib import resources

import scanpy as sc
import torch

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


def load_datasets():
    data_folder = resources.files("hpa_embeddings") / "data"
    adata_dict = {}
    for dset in [
        "251020_blood_atlas.h5ad",
        "251020_mann2017.h5ad",
        "251020_monaco.h5ad",
        "251020_schmiedel.h5ad",
    ]:
        adata = sc.read_h5ad(
            # "/Users/vedran/Downloads/251020_blood_atlas.h5ad"
            data_folder / dset
        )
        cell_type_key = "cell_type"
        assert cell_type_key in adata.obs

        gene_col = "index"
        if gene_col == "index":
            adata.var["index"] = adata.var.index
        else:
            assert gene_col in adata.var

        adata_dict[dset] = adata

    adata_dict.keys()

    return adata_dict, gene_col
