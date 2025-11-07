import json

import numpy as np
from tqdm import tqdm

from bloodtextdata.funcs.dataset_metadatas import metadatas


def create_document_jsons(adata_dict):
    documents = []
    lineages = dict()
    celltypes = dict()
    num_genes_kept = 0
    num_genes = 0
    for dataset_name, adata in adata_dict.items():
        num_genes += adata.X.shape[0] * adata.X.shape[1]
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
            celltype = row.obs["cell_type"].values[0]
            celltypes[celltype] = 0

            row_metadata = dict(
                zip(list(row.obs.keys()), row.obs.values[0])
            )  # celltype and lineage
            row_metadata["type"] = "dataset"
            row_metadata["dataset_name"] = dataset_name
            row_metadata["dataset_description"] = (
                "This is a row from a dataset about Celltypes, their corresponding lineages "
                "and their corresponding gene counts. "
            )

            # TODO: Use Lornas mappings to convert the celltype to synonims and ontologies and add these to the meta data
            # row_metadata["meshterms"] =
            # row_metadata["ontology_tree"] =

            row_genes = dict(
                zip(
                    row.var["index"].values,
                    np.array(row.X[0].data).astype(np.int64).tolist(),
                )
            )
            for key in list(row_genes.keys()):
                if row_genes[key] == 0:
                    del row_genes[key]
            num_genes_kept += len(list(row_genes.keys()))

            # pubmed_id = json_data["pmid"],
            # pmc_id = json_data.get("pmcid", ""),
            # doi = json_data.get("doi") or "",
            # title = json_data["title"],
            # authors = json_data["authorString"],
            # year = json_data["pubYear"],
            # license = json_data.get("license") or "",
            # mesh_terms = ",".join(mesh_terms),

            document = {
                "metadata": metadatas[dataset_name] | row_metadata,
                "body": json.dumps(row_genes),
            }
            documents.append(document)

    print(
        "kept",
        num_genes_kept,
        "out of",
        num_genes,
        "so",
        (num_genes_kept / num_genes) * 100,
        "%",
    )

    print("celltypes", len(list(celltypes.keys())), list(celltypes.keys()))

    print("example document")
    print(documents[0])

    return documents
