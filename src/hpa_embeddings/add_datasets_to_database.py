import time

import sys
import json

import numpy as np
import torch
from tqdm import tqdm

from hpa_embeddings.funcs.create_jsons import create_document_jsons
from hpa_embeddings.funcs.load_datasets import load_datasets
from hpa_embeddings.funcs.insert_into_db import insert_into_db

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
documents = create_document_jsons(adata_dict)
insert_into_db(documents)


time.sleep(10)


# TODO: drop percentage of genes that make up below a threshold percentage of the total gene counts
# TODO: instead of indexing every singlerow of every database only index the database and then have a separate vector database for retrieving the specific row
#   such that the anndata can remain together and sparse

