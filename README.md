# Project 7: Embedding Multimodal Data into a Vector Database to Make a Blood Atlas (BioHackathon 2025)

ELIXIR Biohackathon 2025 project nr. 7

## Project Description

This project aims to develop a multimodal vector embedding dataset that integrates blood cell data across multiple modalities, including text corpora, proteomics profiles, transcriptomic signatures, and cellular imaging, to be made available for future Retrieval-Augmented Generation (RAG) applications. This will enable researchers to use natural language queries to access knowledge about various blood cell types, thus improving quality control of lab-grown cellular blood cultures and a general understanding of cellular types and their dynamics in health and disease.

We want to establish this database because there is a shift in the field towards the usage of multimodal data for the training of deep learning and generative models. At the same time, while many datasets are available for such a well-studied system as blood, a united way to access knowledge is still lacking. The generation of high-quality pre-computed vector embeddings that preserve the inherent relationships between modalities while enabling efficient similarity search and cross-modal queries is a non-trivial problem. This project aims to integrate omics data — such as RNA-Seq and mass spectrometry — with tokens of text obtained from open-source research paper abstracts and with microscopy images generated using segmented blood smear images or imaging flow cytometry. We expect to deliver a curated, validated dataset with standardized formats, comprehensive metadata, and pre-computed vector embeddings that enable seamless integration across modalities.

## Team

- Tess Afanasyeva
- Marina Pominova
- Cecilia Lindskog
- Mahfouz Shehu
- Benjamin Gerritsen
- Emanuel Quadros

## Installation

This repository uses **conda + conda-lock** for a fully reproducible environment across macOS, Linux, and Windows.


### 1. Choose your tool

You can use either **conda** or the faster alternative **micromamba**.

### 2. Create the environment


#### Option 0 — (Windows only) Barebones install
##### Create python env in project directory
```bash
python -m venv venv
.\venv\Scripts\activate
```
##### Run install and reproduce results
```bash 
.\install_and_run_tests.ps1
```

#### Option A — using conda
```bash
# Install conda-lock

# Create environment from lockfile
conda create --name hackathon --file conda-lock.yml

# Or directly from the environment.yml (w/o lockfile)
conda env create -f environment.yml

# Activate
conda activate hackathon
```
#### Option B — using micromamba (recommended)
##### Install micromamba (if not installed)
```bash
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
export PATH=$PWD/bin:$PATH
```

##### Create environment
```bash
micromamba create -n hackathon -f conda-lock.yml
micromamba activate hackathon
```

If you are getting `'micromamba' is running as a subprocess and can't modify the parent shell.` error you should run the following commands before activating the environment:
```bash
micromamba shell init --shell bash --root-prefix=~/.local/share/mamba
source ~/.bashrc
```

### 3. Install PyTorch

We don’t include PyTorch in the lockfile (to avoid GPU/CPU conflicts).
Please install one of the following fixed versions after activating your environment:

#### CPU version
##### Conda
```bash
conda install pytorch=2.4.1 torchvision=0.19.1 cpuonly -c pytorch -c conda-forge
```

##### Micromamba
```bash
micromamba install pytorch=2.4.1 torchvision=0.19.1 cpuonly -c pytorch -c conda-forge
```

#### GPU version (CUDA 12.4)
##### Conda
```bash
conda install pytorch=2.4.1 torchvision=0.19.1 pytorch-cuda=12.4 -c pytorch -c nvidia -c conda-forge
```

##### Micromamba
```bash
micromamba install pytorch=2.4.1 torchvision=0.19.1 pytorch-cuda=12.4 -c pytorch -c nvidia -c conda-forge
```

## Project Structure

```
biohack2025/
├── data/               # Data files (raw and processed)
├── models/             # Trained models and model artifacts
├── notebooks/          # Jupyter notebooks for analysis
├── src/                # Source code modules
├── tests/              # Test files
└── environment.yml     # Project dependencies

