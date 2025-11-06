import time

import sys
import json

import numpy as np
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

documents = []
lineages = dict()
celltypes = dict()
for dataset_name, adata in adata_dict.items():
    num_genes = adata.X.shape[0] * adata.X.shape[1]
    num_genes_kept = 0
    for row in tqdm(adata):
        # row.write('my_results.h5ad', compression="gzip")
        # print("row", row)
        # print("row.obs", list(row.obs.keys()))
        # print("row.obs", row.obs.values)
        # print("row.var", row.X)
        # print("row.var")
        # print(row.var)
        # print("Genes", list(row.var.keys()))
        # print("Genes:", len(row.var['index'].values), row.var['index'].values)

        # Obtain all celltypes
        celltype = row.obs['cell_type'].values[0]
        celltypes[celltype] = 0

        row_metadata = dict(zip(list(row.obs.keys()), row.obs.values[0])) # celltype and lineage
        row_metadata["type"] = "dataset"
        row_metadata["dataset_name"] = dataset_name
        row_metadata["dataset_description"] = \
            "This is a row from a dataset about Celltypes, their correpsonding lineages " \
            "and their corresponding gene counts. "
        # TODO: Use Lornas mappings to convert the celltype to synonims and ontologies and add these to the meta data
        # row_metadata["meshterms"] =
        # row_metadata["ontology_tree"] =

        row_genes = dict(zip(row.var['index'].values, np.array(row.X[0].data).astype(np.int64).tolist()))
        for key in list(row_genes.keys()):
            if row_genes[key] == 0:
                del row_genes[key]
        num_genes_kept += len(list(row_genes.keys()))

        document = {
            "metadata": {
                "pmcid": None,
                "author": None,
                "License": None,
                "source_url": None,
                "dataset_meta": json.dumps(row_metadata)},
            "body": json.dumps(row_genes)}
        documents.append(document)

    print("kept", num_genes_kept, "out of", num_genes, "so", (num_genes_kept / num_genes) * 100, "%")

print("celltypes", len(list(celltypes.keys())), list(celltypes.keys()))

print("example document")
print(documents[0])

time.sleep(10)


# TODO: drop percentage of genes that make up below a threshold percentage of the total gene counts
# TODO: instead of indexing every singlerow of every database only index the database and then have a separate vector database for retrieving the specific row
#   such that the anndata can remain together and sparse
