# SCORED

SCORED is a Python package for single-cell RNA sequencing data imputation using SimRank and Random Walk with Restart algorithms. It processes AnnData objects, performs quality control, builds graphs, and imputes missing data.

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

## License

MIT License
