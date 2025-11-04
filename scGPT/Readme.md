# scGPT Embeddings

A specialized environment for scGPT reference mapping and cell embedding generation.

## Setup

This subproject uses Poetry for dependency management:

```bash
# Install dependencies
poetry install

# Activate the environment
poetry shell

# Run the tutorial notebook
poetry run jupyter notebook Tutorial_Reference_Mapping_dataset.ipynb
```

## Requirements

- Python 3.10+
- Optimized for Mac M1 Pro (ARM64)
- Uses MPS backend for GPU acceleration
