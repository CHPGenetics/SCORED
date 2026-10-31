# SCORED <img align="right" style="margin-left: 20px; margin-bottom: 10px;" src="./pictures/sticker_SCORED.png" width="195" height="215">

![GitHub](https://img.shields.io/github/license/CHPGenetics/SCORED) [![pypiDownloads](https://static.pepy.tech/badge/scored)](https://pepy.tech/project/scored)

Long-read single-cell sequencing offers a powerful opportunity to resolve full-length isoform expression, but the resulting data are often characterized by profound sparsity and high levels of technical noise (dropouts). This sparsity obscures important biological signals, such as correlations between isoforms and dynamic expression trends along cellular trajectories.

`SCORED` (Single-Cell isOform data REfinement via graph-based Diffusion) is a computational method tailored for the refinement of sparse single-cell RNA isoform expression data generated from long-read sequencing technologies. It addresses the data sparsity challenge by implementing a graph-based diffusion algorithm that borrows information from functionally similar cells to infer a more complete and accurate isoform expression profile for each cell. 

This repository accompanies the manuscript “Time-resolved single-cell multi-omics reveals CREM isoform dynamics underlying PGE2-mediated inhibition of CD4⁺ memory T cell activation,” providing all associated data processing and analytical workflows. The reproducible results notebook is available [here](https://dogmaisoformpaperresults.netlify.app).

![Overview of the SCORED algorithm](./pictures/overview.png)


## Installation

```bash
pip install scanpy
pip install SCORED
```

## Running SCORED

### Quick start

The `SCORED` function can be applied directly to raw AnnData objects, with no prior normalization or preprocessing required.

```python
import scanpy as sc
import torch
from scored import SCORED

# Load your AnnData object
adata_transcript = sc.read_h5ad("isoform_level_data.h5ad") # Isoform-level UMI matrix
adata_gene = sc.read_h5ad("gene_level_data.h5ad") # Gene-level UMI matrix that accompanies isoform-level data 
# (e.g., from preprocessing pipelines such as wf-single-cell).
# By default, gene-level data is used to infer cell similarity.

# Run SCORED
refined_matrix = SCORED(
    adata_tr = adata_transcript,
    adata_g = adata_gene
)

# To replace the original data with the refined data
adata_transcript_refined = adata_transcript.copy()
adata_transcript_refined.X = refined_matrix
# Data is normalized. Use log transformation to adjust the scale if needed for your analysis.
```

### Advanced tutorial

#### 1. Isoform-level cell similarity (optional)

`SCORED` can also use isoform-level expression to infer cell similarity (recommended only when sequencing depth is high and isoform-level data is of good quality):

```python
refined_matrix = SCORED(
    adata_tr = adata_transcript
)
```
#### 2. Condition-specific refinement

If the AnnData object contains multiple conditions and you want to preserve subtle, condition-specific expression signatures, you can refine within each condition:

```python
refined_matrix = SCORED(
    adata_tr = adata_transcript,
    adata_g = adata_gene,
    condition_key="condition"
)
```
#### 3. Tuning parameters

The key hyperparameter for `SCORED` is the restart probability parameter `rwr_alpha`. Default is `rwr_alpha = 0.3`,  which is robust in most cases. 
* Smaller values: cells borrow more information from neighbors.
* Larger values: more local information is preserved.

```python
refined_matrix = SCORED(
    adata_tr = adata_transcript,
    adata_g = adata_gene,
    rwr_alpha = 0.3
) 
```

#### 4. GPU acceleration

`SCORED` can leverage GPU for faster runtime. If a GPU is available, set `device="cuda"`:

```python
refined_matrix = SCORED(
    adata_tr = adata_transcript,
    adata_g  = adata_gene,
    device   = "cuda"
)
```


## System requirements

### Hardware
SCORED requires only a standard computer with enough RAM for the
in-memory operations. A CUDA-capable GPU is optional and
speeds up computation on larger datasets.

### Operating system
Tested on Linux (Ubuntu 22.04).

### Software dependencies
Python (tested on 3.13). Minimum required versions, with the versions
SCORED has been tested on:

| Package      | Required | Tested  |
|--------------|----------|---------|
| numpy        | ≥1.20    | 2.1.3   |
| scipy        | ≥1.7     | 1.16.3  |
| pandas       | ≥1.3     | 3.1.0   |
| scikit-learn | ≥1.0     | 1.6.1   |
| torch        | ≥1.9     | 2.11.0  |
| scanpy       | ≥1.8     | 1.12.5  |
| anndata      | —        | 0.13.4  |
| networkx     | ≥2.6     | 3.7     |
| matplotlib   | ≥3.4     | 3.10.0  |
| seaborn      | ≥0.11    | 0.13.2  |
| tqdm         | ≥4.0     | 4.67.3  |

## Installation time
Typical install time: ~1 minute on a standard desktop computer

