import os

from download_europepmc import main as download_data
from update_body_from_xml import main as xml_patching

if __name__ == "__main__":
    download_data()
    xml_patching()
