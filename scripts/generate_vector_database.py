import sys
from os.path import join, dirname
import argparse
import json
import pathlib

from pprint import pp

sys.path.append(join(dirname(dirname(__file__)), "src"))

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

    articles_dir = (
        args.articles_dir
        if args.articles_dir is not None
        else join(dirname(dirname(__file__)), "data", "europepmc_articles")
    )
    db_path = (
        args.db_path
        if args.db_path is not None
        else join(dirname(dirname(__file__)), "data", "database")
    )

    json_files = tuple(pathlib.Path(articles_dir).glob("*.json"))
    embedder = Embedder("michiyasunaga/BioLinkBERT-large", device="cuda")
    vectordb = VectorDB(db_path=db_path, embedder=embedder)

    for file in tqdm(json_files):
        with file.open(encoding="utf-8") as fp:
            try:
                json_data = json.load(fp)
            except UnicodeError:
                print("UnicodeError:", fp)
                continue

            if (
                "body" in json_data
                and "pmid" in json_data
                and json_data["body"]
                and vectordb.get_article_vector(
                    pubmed_id=int(json_data["pmid"]), segment_number=0
                )
                is None
            ):  # Check if article already indexed for resuming
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
                    metadata_prefix = f"""
                    Title: {metadata.title}
                    Authors: {metadata.authors}
                    Year: {metadata.year}
                    MeSH Terms: {metadata.mesh_terms}
                    
                    """
                    full_text = metadata_prefix + json_data["body"]["text"]
                    vectordb.add_article_vectors(text=full_text, metadata=metadata)


if __name__ == "__main__":
    main()
