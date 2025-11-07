import sys

import torch
import sys
from os.path import dirname
sys.path.append(dirname(dirname(__file__)))

from vectordb import VectorDB

from bloodtextdata.funcs.create_jsons import create_document_jsons
from bloodtextdata.funcs.load_datasets import load_datasets

print(torch.__version__)

adata_dict, gene_col = load_datasets()
documents = create_document_jsons(adata_dict)

vdb = VectorDB()
vdb.add_blood_text_data_bulk(metadatas=[doc["metadata"] for doc in documents],
                             bodies=[doc["body"] for doc in documents])


# TODO: drop percentage of genes that make up below a threshold percentage of the total gene counts
# TODO: instead of indexing every singlerow of every database only index the database and then have a separate vector database for retrieving the specific row
#   such that the anndata can remain together and sparse
