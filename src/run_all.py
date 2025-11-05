import os

from download_europepmc import main as download_data
from embeddings import Embedder
from update_body_from_xml import main as xml_patching, find_repo_root
from vectordb.vectordb import VectorDB

#
# def fill_database(articles_directory):
#     embedder = Embedder()
#
#     for article_name in os.listdir(articles_directory):
#

if __name__ == "__main__":
    download_data()
    xml_patching()
    #
    # repo_root = find_repo_root(os.getcwd())
    # database = VectorDB(os.path.join(repo_root, 'data', 'database'))
    # fill_database(os.path.join(repo_root, 'data', 'europepmc_articles'))
