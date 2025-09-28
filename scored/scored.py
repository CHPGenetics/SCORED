import sklearn
import scipy
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from tqdm.notebook import tqdm
import networkx as nx
from scipy.stats import spearmanr
import numpy as np
import torch
import scanpy as sc
import pandas as pd
from scipy.sparse import csr_matrix
from anndata import AnnData
from typing import Optional, List

def SCORED(
    adata_tr: AnnData,
    adata_g: Optional[AnnData]   = None,
    condition_key: Optional[str] = None,
    # --- Scanpy‐preprocessing params ---
    filter_min_cells: int        = 3,
    target_sum: float            = 1e4,
    regress_vars: List[str]      = ["total_counts", "pct_counts_mt"],
    n_neighbors: int             = 60,
    n_pcs: int                   = 40,
    # --- SimRank params ---
    k: int                       = 30,
    simrank_decay: float         = 0.8,
    simrank_max_iter: int        = 100,
    simrank_eps: float           = 1e-4,
    # --- RWR params ---
    rwr_alpha: float             = 0.3,
    rwr_max_iter: int            = 100,
    rwr_tol: float               = 1e-4,
    # --- device ---
    device: Optional[str]        = None
):
    """
      1) Performs comprehensive QC on adata_tr.
      2) Optionally processes adata_g (or clones adata_tr).
      3) Builds a masked graph.
      4) Runs an optimized SimRank algorithm.
      5) Computes refined Gaussian–kernel weights.
      6) Conducts a Random Walk With Restart for imputation.
      7) Returns the refined matrix.
    """

    # === STEP 0: Device initialization ===
    # Choose computing device intelligently
    if device is None:
        cuda_available = torch.cuda.is_available()
        if cuda_available:
            selected_device = "cuda"
        else:
            selected_device = "cpu"
    else:
        selected_device = device
    # Create torch.device
    torch_device = torch.device(selected_device)


    # === STEP 1: Transcriptome preprocessing ===
    # 1.1: Filter out genes with very low counts
    sc.pp.filter_genes(adata_tr, min_cells=filter_min_cells)
    # 1.2: Annotate mitochondrial genes for QC
    adata_tr.var["mt"] = adata_tr.var_names.str.startswith("MT-")
    # 1.3: Compute QC metrics (percent mitochondrial, etc.)
    sc.pp.calculate_qc_metrics(
        adata_tr,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True
    )
    # 1.4: Normalize total counts per cell
    sc.pp.normalize_total(adata_tr, target_sum=target_sum)
    # 1.5: Keep a raw copy for later imputation
    adata_tr_raw_copy = adata_tr.copy()
    # 1.6: Log‐transform the data
    sc.pp.log1p(adata_tr)
    # 1.7: Identify highly variable genes
    sc.pp.highly_variable_genes(adata_tr)
    # 1.8
    _ = adata_tr.var.highly_variable.sum()



    # === STEP 2: Graph AnnData preparation ===
    if adata_g is not None:
        # 2.1: Apply the same QC pipeline to adata_g
        sc.pp.filter_genes(adata_g, min_cells=filter_min_cells)
        adata_g.var["mt"] = adata_g.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(
            adata_g,
            qc_vars=["mt"],
            percent_top=None,
            log1p=False,
            inplace=True
        )
        sc.pp.normalize_total(adata_g, target_sum=target_sum)
        sc.pp.log1p(adata_g)
        sc.pp.highly_variable_genes(adata_g)
        # Subset to HVGs
        hvg_mask = adata_g.var.highly_variable
        adata_g = adata_g[:, hvg_mask]
    else:
        # 2.2: clone the transcriptome
        adata_g = adata_tr.copy()

    # 2.3: Regress out unwanted sources of variation
    sc.pp.regress_out(adata_g, regress_vars)
    # 2.4: Compute neighborhood graph
    sc.pp.neighbors(adata_g, n_neighbors=n_neighbors, n_pcs=n_pcs)
    # 2.5: Embed with UMAP
    sc.tl.umap(adata_g)
    _umap_coords = adata_g.obsm["X_umap"]

    # === STEP 3: Extract adjacency & distances ===
    adjacency_sparse = adata_g.obsp["connectivities"]
    distance_sparse  = adata_g.obsp["distances"]
    num_cells       = adjacency_sparse.shape[0]

    # === STEP 4: Build condition mask (optional) ===
    if condition_key is None:
        # No condition filtering: allow all edges
        adjacency_mask = np.ones((num_cells, num_cells), dtype=bool)
    else:
        # Prefer transcriptome metadata
        if condition_key in adata_tr.obs:
            condition_array = adata_tr.obs[condition_key].astype(str).values
        elif condition_key in adata_g.obs:
            condition_array = adata_g.obs[condition_key].astype(str).values
        else:
            # Fallback: no mask
            condition_array = None

        if condition_array is not None:
            adjacency_mask = np.equal.outer(condition_array, condition_array)
        else:
            adjacency_mask = np.ones((num_cells, num_cells), dtype=bool)


    # Convert sparse→dense, apply mask
    adjacency_dense = adjacency_sparse.toarray()
    adjacency_dense[~adjacency_mask] = 0
    A_tensor = torch.tensor(adjacency_dense, dtype=torch.float32, device=torch_device)


    # === STEP 5: SimRank computation ===
    def _simrank_internal(adj_tensor, decay_factor, max_iterations, epsilon):
        # initialize similarity matrix
        n_nodes = adj_tensor.size(0)
        sim = torch.eye(n_nodes, device=adj_tensor.device)
        # normalize adjacency by out‐degree
        outdeg = adj_tensor.sum(dim=1, keepdim=True)
        M = torch.where(outdeg > 0, adj_tensor / outdeg, adj_tensor)
        # iterative update
        for iteration in range(max_iterations):
            update_term = M @ sim @ M.T
            new_sim   = decay_factor * update_term
            new_sim.fill_diagonal_(1.0)
            diff_norm = torch.norm(new_sim - sim)
            if diff_norm < epsilon:
                break
            sim = new_sim
        return sim


    sim_matrix = _simrank_internal(
        A_tensor,
        simrank_decay,
        simrank_max_iter,
        simrank_eps
    )
    sim_numpy = (sim_matrix + sim_matrix.T).cpu().numpy()


    # === STEP 6: Gaussian‐kernel weights ∩ top-k SimRank neighbors ===
    def _gaussian_kernel_intersection(dist_sp: csr_matrix, sim_arr: np.ndarray, k_nn: int):
        D = dist_sp.toarray()
        Dt = torch.tensor(D, dtype=torch.float32, device=torch_device)
        n   = Dt.size(0)
        # compute per-cell sigma^2
        sigma2 = torch.zeros(n, device=torch_device)
        for idx in range(n):
            nonzero_vals = Dt[idx][Dt[idx] > 0]
            if nonzero_vals.numel() > 0:
                sigma2[idx] = torch.median(nonzero_vals ** 2)
        sigma = torch.sqrt(sigma2)
        # allocate weight matrix
        W = torch.zeros_like(Dt)
        # fill via intersection
        for i in range(n):
            knn_neighbors = torch.where(Dt[i] > 0)[0].cpu().numpy()
            topk_by_sim    = np.argsort(-sim_arr[i])[:k_nn]
            common_inds    = np.intersect1d(knn_neighbors, topk_by_sim)
            for j in common_inds:
                dij2 = Dt[i, j] ** 2
                numerator   = 2 * sigma[i] * sigma[j]
                denominator = sigma2[i] + sigma2[j]
                W[i, j] = torch.sqrt(numerator / denominator) * torch.exp(-dij2 / denominator)
        return W + W.T


    W_tensor = _gaussian_kernel_intersection(distance_sparse, sim_numpy, k)

    # === STEP 7: Combine with SimRank & normalize rows ===
    combined_W = W_tensor * torch.tensor(sim_numpy, dtype=torch.float32, device=torch_device)
    combined_np = combined_W.cpu().numpy()
    W_df = pd.DataFrame(combined_np)
    row_sums = W_df.sum(axis=1)
    W_normalized = W_df.div(row_sums, axis=0).fillna(0)


    # === STEP 8: Random Walk With Restart (RWR) ===
    def _random_walk_restart(P_df: pd.DataFrame, raw_adata: AnnData,
                             alpha, max_iters, tolerance):
        P_mat = torch.tensor(P_df.values, dtype=torch.float32, device=torch_device)
        n_dim = P_mat.size(0)
        P0    = torch.eye(n_dim, device=torch_device)
        Pcurr = P0.clone()
        for it in range(max_iters):
            Pnext = (1 - alpha) * (Pcurr @ P_mat) + alpha * P0
            if torch.norm(Pnext - Pcurr) < tolerance:
                break
            Pcurr = Pnext
        P_final = Pcurr.cpu().numpy()
        # refine expression
        X_raw = raw_adata.X.toarray() if hasattr(raw_adata.X, "toarray") else np.array(raw_adata.X)
        X_t   = torch.tensor(X_raw, dtype=torch.float32, device=torch_device)
        imputed = (Pcurr.T @ X_t).cpu().numpy()


        return imputed, P_final


    new_matrix, P_final = _random_walk_restart(
        W_normalized,
        adata_tr_raw_copy,
        rwr_alpha,
        rwr_max_iter,
        rwr_tol
    )

    # return final result
    return new_matrix
