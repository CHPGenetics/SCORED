# SCORED <img align="right" style="margin-left: 20px; margin-bottom: 10px;" src="./pictures/sticker.png" width="180" height="195">

![GitHub](https://img.shields.io/github/license/CHPGenetics/SCORED) [![pypiDownloads](https://static.pepy.tech/badge/scored)](https://pepy.tech/project/scored)[![Pytest](https://github.com/CHPGenetics/SCORED/workflows/py308/badge.svg)](https://github.com/CHPGenetics/SCORED)

Long-read single-cell sequencing provides a foundational tool opportunity to resolve full-length isoform expression, but the resulting data is often characterized by profound sparsity and high levels of technical noise (dropouts). This sparsity obscures important biological signals, such as correlations between isoforms and dynamic expression trends along cellular trajectories.

`SCORED` is a computational method tailored for the refinement of sparse single-cell RNA isoform expression data generated from long-read sequencing technologies. It addresses the data sparsity challenge by implementing a graph-based diffusion algorithm that borrows information from functionally similar cells to infer a more complete and accurate isoform expression profile for each cell. 

This repository contains the Python implementation, example usage, and experiment results of the `SCORED` algorithm.

## Installation

```bash
pip install SCORED
```

## Usage

The `SCORED` function can be applied directly to raw Scanpy objects, with no prior normalization or preprocessing required.

```python
import scanpy as sc
import torch
from scored import SCORED

# Load your AnnData object
adata = sc.read_h5ad("your_data.h5ad")

# Run SCORED
imputed_matrix = SCORED(
    adata_tr=adata,
    condition_key="condition",
    device="cuda" if torch.cuda.is_available() else "cpu"
)
```

## Requirements

- Python >= 3.8
- scikit-learn >= 1.0
- scipy >= 1.7
- matplotlib >= 3.4
- seaborn >= 0.11
- tqdm >= 4.0
- networkx >= 2.6
- numpy >= 1.20
- torch >= 1.9
- scanpy >= 1.8
- pandas >= 1.3


