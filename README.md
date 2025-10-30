# Project 7: Embedding Multimodal Data into a Vector Database to Make a Blood Atlas (BioHackathon 2025)

ELIXIR Biohackathon 2025 project nr. 7

## Project Description

This project aims to develop a multimodal vector embedding dataset that integrates blood cell data across multiple modalities, including text corpora, proteomics profiles, transcriptomic signatures, and cellular imaging, to be made available for future Retrieval-Augmented Generation (RAG) applications. This will enable researchers to use natural language queries to access knowledge about various blood cell types, thus improving quality control of lab-grown cellular blood cultures and a general understanding of cellular types and their dynamics in health and disease.

We want to establish this database because there is a shift in the field towards the usage of multimodal data for the training of deep learning and generative models. At the same time, while many datasets are available for such a well-studied system as blood, a united way to access knowledge is still lacking. The generation of high-quality pre-computed vector embeddings that preserve the inherent relationships between modalities while enabling efficient similarity search and cross-modal queries is a non-trivial problem. This project aims to integrate omics data — such as RNA-Seq and mass spectrometry — with tokens of text obtained from open-source research paper abstracts and with microscopy images generated using segmented blood smear images or imaging flow cytometry. We expect to deliver a curated, validated dataset with standardized formats, comprehensive metadata, and pre-computed vector embeddings that enable seamless integration across modalities.

## Team

- Tess Afanasyeva
- Marina Pominova
- Cecilia Lindskog

## Installation

### Prerequisites

- Docker and Docker Compose (recommended)
- OR Python 3.10+ and Poetry 1.7+ (for local installation)

### Option 1: Docker Installation (Recommended)

This is the easiest way to get started with a fully configured environment including Jupyter Lab.

1. Clone the repository:
```bash
git clone <repository-url>
cd biohack2025
```

2. Build and start the Docker container:
```bash
docker-compose up --build
```

3. Access Jupyter Lab in your browser at:
```
http://localhost:8888
```

The container includes all dependencies and mounts your local directory, so changes are immediately reflected.

#### Docker Commands

- Start the container: `docker-compose up`
- Start in detached mode: `docker-compose up -d`
- Stop the container: `docker-compose down`
- Rebuild after dependency changes: `docker-compose up --build`
- Access shell in running container: `docker exec -it blood-atlas-dev bash`

### Option 2: Local Installation with Poetry

1. Install Poetry if not already installed:
```bash
curl -sSL https://install.python-poetry.org | python3 -
```

2. Clone the repository:
```bash
git clone <repository-url>
cd biohack2025
```

3. Install dependencies:
```bash
poetry install
```

4. Activate the virtual environment:
```bash
poetry shell
```

5. Start Jupyter Lab:
```bash
poetry run jupyter lab
```

### Option 3: Local Installation with pip

If you prefer not to use Poetry, you can export the dependencies:

```bash
poetry export -f requirements.txt --output requirements.txt --without-hashes
pip install -r requirements.txt
```

## Project Structure

```
biohack2025/
├── data/               # Data files (raw and processed)
├── models/             # Trained models and model artifacts
├── notebooks/          # Jupyter notebooks for analysis
├── src/                # Source code modules
├── tests/              # Test files
├── download_europepmc.py  # Script to download papers from Europe PMC
├── pyproject.toml      # Poetry dependencies and project config
├── Dockerfile          # Docker container definition
└── docker-compose.yml  # Docker Compose configuration
```

## Usage

### Running Jupyter Notebooks

After installation, you can create and run Jupyter notebooks in the `notebooks/` directory. The environment includes all necessary libraries for:

- Data processing and analysis (pandas, numpy, scipy)
- Machine learning (scikit-learn, PyTorch)
- Vector embeddings (sentence-transformers, ChromaDB, FAISS)
- Image processing (OpenCV, scikit-image)
- Bioinformatics (Biopython, Scanpy, AnnData)
- Visualization (matplotlib, seaborn, plotly)

### Running Scripts

Example - downloading papers from Europe PMC:
```bash
# With Docker
docker exec -it blood-atlas-dev poetry run python download_europepmc.py

# With Poetry (local)
poetry run python download_europepmc.py

# With activated virtual environment
python download_europepmc.py
```

## Development

### Installing Development Dependencies

Development dependencies (pytest, black, ruff, mypy) are included by default with Poetry installation.

### Code Formatting and Linting

```bash
# Format code with Black
poetry run black .

# Lint with Ruff
poetry run ruff check .

# Type checking with MyPy
poetry run mypy src/
```

### Running Tests

```bash
poetry run pytest
```

## Adding New Dependencies

```bash
# Add a production dependency
poetry add package-name

# Add a development dependency
poetry add --group dev package-name

# After adding dependencies in Docker, rebuild the container
docker-compose up --build
```

## Troubleshooting

### Docker Issues

- If port 8888 is already in use, modify the port mapping in `docker-compose.yml`
- If you encounter permission issues, ensure Docker has access to your project directory

### Poetry Issues

- Clear Poetry cache: `poetry cache clear pypi --all`
- Update Poetry: `poetry self update`
- Recreate virtual environment: `poetry env remove python && poetry install`

## License

[Add license information]

## Acknowledgments

This project is part of ELIXIR Biohackathon 2025.
