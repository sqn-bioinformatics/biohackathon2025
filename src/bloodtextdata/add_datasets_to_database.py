import json
import sys
import time
from importlib import resources

import numpy as np
import torch
from tqdm import tqdm
from vectordb import VectorDB

from bloodtextdata.funcs.create_jsons import create_document_jsons
from bloodtextdata.funcs.load_datasets import load_datasets

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

vdb = VectorDB()

adata_dict, gene_col = load_datasets()
documents = create_document_jsons(adata_dict)

for doc in tqdm(documents):
    vdb.add_blood_text_data(metadata=doc["metadata"]["dataset_meta"], body=doc["body"])


# TODO: drop percentage of genes that make up below a threshold percentage of the total gene counts
# TODO: instead of indexing every singlerow of every database only index the database and then have a separate vector database for retrieving the specific row
#   such that the anndata can remain together and sparse
