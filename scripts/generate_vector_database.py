import argparse
import json
import pathlib
from os.path import join, dirname
from pprint import pp

from embeddings import Embedder
from tqdm import tqdm
from vectordb import TextMetadata, VectorDB


def main():
    parser = argparse.ArgumentParser(description="Create vector database.")
    parser.add_argument(
        "--articles-dir", default=None, help="Path to europepmc_articles directory"
    )
    parser.add_argument(
        "--db_path",
        default="vectors.chromadb",
        help="Path to ChromaDB store.",
    )
    args = parser.parse_args()

    articles_dir = args.articles_dir if args.articles_dir is not None else dirname(join("data", "europepmc_articles"))
    db_path = args.db_path if args.db_path is not None else dirname(join("data", "database"))

    # json_files = tuple(pathlib.Path(articles_dir).glob("*.json"))
    embedder = Embedder("michiyasunaga/BioLinkBERT-base", device="cuda")
    vectordb = VectorDB(db_path=db_path, embedder=embedder)
    print(tuple(pathlib.Path(articles_dir).glob("*.json")))
    for file in tqdm(pathlib.Path(articles_dir).glob("*.json")):
        with file.open() as fp:
            json_data = json.load(fp)
            if "body" in json_data and "pmid" in json_data and json_data["body"]:
                mesh_terms = []
                try:
                    mesh_terms = [
                        mesh_heading["descriptorName"]
                        for mesh_heading in json_data["meshHeadingList"]["meshHeading"]
                    ]
                except KeyError:
                    pass

                try:
                    metadata = TextMetadata(
                        pubmed_id=json_data["pmid"],
                        pmc_id=json_data.get("pmcid", ""),
                        doi=json_data.get("doi") or "",
                        title=json_data["title"],
                        authors=json_data["authorString"],
                        year=json_data["pubYear"],
                        license=json_data.get("license") or "",
                        mesh_terms=",".join(mesh_terms),
                    )
                except KeyError:
                    print(json_data["pmid"])
                else:
                    vectordb.add_article_vectors(
                        text=json_data["body"]["text"], metadata=metadata
                    )


if __name__ == "__main__":
    main()
