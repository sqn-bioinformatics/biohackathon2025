metadatas = \
    {
        "251020_blood_atlas.h5ad": {
            "pubmed_id": None,
            "pmc_id": None,
            "doi": None,
            "title": None,
            "authors": None,
            "year": None,
            "license": None
        },
        "251020_mann2017.h5ad": {
            "pubmed_id": None,
            "pmc_id": None,
            "doi": None,
            "title": None,
            "authors": None,
            "year": None,
            "license": None
        },
        "251020_monaco.h5ad": {
            "pmc_id": None,
            "doi": None,
            "title": None,
            "authors": None,
            "year": None,
            "license": None
        },
        "251020_schmiedel.h5ad": {
            "pubmed_id": None,
            "pmc_id": None,
            "doi": None,
            "title": None,
            "authors": None,
            "year": None,
            "license": None
        }
    }

for key in metadatas.keys():
    for key2 in metadatas[key]:
        if metadatas[key][key2] is None:
            metadatas[key][key2] = "None"