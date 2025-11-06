import argparse
import json
import pathlib
import sys
from os.path import dirname, join
from pprint import pp

sys.path.append(join(dirname(dirname(__file__)), "src"))

from tqdm import tqdm
from vectordb import TextMetadata, VectorDB


def main():
    parser = argparse.ArgumentParser(description="Create vector database.")
    parser.add_argument(
        "--articles-dir", default=None, help="Path to europepmc_articles directory"
    )
    args = parser.parse_args()

    articles_dir = (
        args.articles_dir
        if args.articles_dir is not None
        else join(dirname(dirname(__file__)), "data", "europepmc_articles")
    )

    json_files = tuple(pathlib.Path(articles_dir).glob("*.json"))
    vectordb = VectorDB()

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
                and not vectordb.article_exists(json_data["pmid"])
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
                        abstract=json_data.get("abstract", ""),
                    )
                except KeyError:
                    print(json_data["pmid"])
                else:
                    metadata_prefix = f"""
                    Title: {metadata.title}
                    Authors: {metadata.authors}
                    Year: {metadata.year}
                    MeSH Terms: {metadata.mesh_terms}
                    Abstract: {metadata.abstract}
                    
                    """
                    full_text = metadata_prefix + json_data["body"]["text"]
                    vectordb.add_article_vectors(text=full_text, metadata=metadata)


if __name__ == "__main__":
    main()
