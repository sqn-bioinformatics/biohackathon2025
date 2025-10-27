# Project 7: Embedding Multimodal Data into a Vector Database to Make a Blood Atlas (BioHackathon 2025)

ELIXIR Biohackathon 2025 project nr. 7

## Project Description

This project aims to develop a multimodal vector embedding dataset that integrates blood cell data across multiple modalities, including text corpora, proteomics profiles, transcriptomic signatures, and cellular imaging, to be made available for future Retrieval-Augmented Generation (RAG) applications. This will enable researchers to use natural language queries to access knowledge about various blood cell types, thus improving quality control of lab-grown cellular blood cultures and a general understanding of cellular types and their dynamics in health and disease.

We want to establish this database because there is a shift in the field towards the usage of multimodal data for the training of deep learning and generative models. At the same time, while many datasets are available for such a well-studied system as blood, a united way to access knowledge is still lacking. The generation of high-quality pre-computed vector embeddings that preserve the inherent relationships between modalities while enabling efficient similarity search and cross-modal queries is a non-trivial problem. This project aims to integrate omics data — such as RNA-Seq and mass spectrometry — with tokens of text obtained from open-source research paper abstracts and with microscopy images generated using segmented blood smear images or imaging flow cytometry. We expect to deliver a curated, validated dataset with standardized formats, comprehensive metadata, and pre-computed vector embeddings that enable seamless integration across modalities.
