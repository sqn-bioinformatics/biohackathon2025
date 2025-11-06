import time

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
from tqdm import tqdm

from hpa_embeddings.load_datasets import load_datasets

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


adata_dict, gene_col = load_datasets()

print("print(adata_dict.keys())", adata_dict.keys())
print('adata_dict["251020_mann2017.h5ad"]', adata_dict["251020_mann2017.h5ad"])
print('adata.obs["cell_type"]', adata_dict["251020_mann2017.h5ad"].obs)
# print('adata.obs["cell_type"]', adata_dict["251020_mann2017.h5ad"].obs["cell_type"])
print('adata_dict["251020_mann2017.h5ad"].var', adata_dict["251020_mann2017.h5ad"].var)
print('adata_dict["251020_mann2017.h5ad"].uns', adata_dict["251020_mann2017.h5ad"].uns.keys())



adata = adata_dict["251020_mann2017.h5ad"]

# count_matrix = (
#     adata.layers["counts"]
#     if isinstance(adata.layers["counts"], np.ndarray)
#     else adata.layers["counts"].A
# )
#
# print("count_matrix", count_matrix)
# print("count_matrix", count_matrix.shape)


# # Ensure count_matrix is a proper numpy array with copy to avoid issues
# count_matrix = np.array(count_matrix, dtype=np.float32)
#
# # gene vocabulary ids
# if gene_ids is None:
#     gene_ids = np.array(adata.var["id_in_vocab"])
#     assert np.all(gene_ids >= 0)
#
# # Ensure gene_ids is a proper numpy array
# gene_ids = np.array(gene_ids, dtype=np.int64)


# print('adata', adata.layers.keys())
print("adata.obs_names", len(adata.obs_names), adata.obs_names)
print("adata.var_names", len(adata.var_names), adata.var_names)

print("adata", adata)

rows = []
lineages = dict()
celltypes = dict()
for dataset_name, adata in adata_dict.items():
    num_genes = adata.X.shape[0] * adata.X.shape[1]
    num_genes_kept = 0
    for row in tqdm(adata):
        row.write('my_results.h5ad', compression="gzip")
        # print("row", row)
        # print("row.obs", list(row.obs.keys()))
        # print("row.obs", row.obs.values)
        # print("row.var", row.X)
        # print("row.var")
        # print(row.var)
        # print("Genes", list(row.var.keys()))
        # print("Genes:", len(row.var['index'].values), row.var['index'].values)
        celltypes[row.obs['cell_type'][0]] = 0
        row_metadata = dict(zip(list(row.obs.keys()), row.obs.values[0])) # celltype and lineage
        row_metadata["type"] = "data"
        row_metadata["dataset name"] = dataset_name
        # row_metadata["meshterms"] =
        # row_metadata["ontology_tree"] =
        row_genes = dict(zip(row.var['index'].values, row.X[0].data))

        for key in list(row_genes.keys()):
            if row_genes[key] == 0:
                del row_genes[key]
        # print("num genes left", len(list(row_genes.keys())))
        num_genes_kept += len(list(row_genes.keys()))
        row = str(row_metadata) + " " + str(row_genes)
        # print("result", row)
        # print("Genes", row.var["index"][:, 0])

        rows.append(row)
    print("kept", num_genes_kept, "out of", num_genes, "so", (num_genes_kept / num_genes) * 100, "%")

print("celltypes", len(list(celltypes.keys())), list(celltypes.keys()))

# for row in adata.to_df():
#     print("row")
#     print(row)

time.sleep(10)


# TODO: drop percentage of genes that make up below a threshold percentage of the total gene counts
# TODO: instead of indexing every singlerow of every database only index the database and then have a separate vector database for retrieving the specific row
#   such that the anndata can remain together and sparse
