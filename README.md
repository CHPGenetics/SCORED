# SCORED <img align="right" style="margin-left: 20px; margin-bottom: 10px;" src="./pictures/sticker.png" width="215" height="215">

SCORED is a Python package for long read single-cell RNA sequencing data refinement. The SCORED method evaluates cell-cell similarity by combining SimRank similarity metrics with Gaussian kernel weights derived from cell-cell distances in the reduced-dimensional space and then leverage a Markov process to incorporate information from similar cells, enabling the inference of the true transcriptomic profile for each cell.

## Installation

```bash
pip install SCORED
```

## Usage

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


