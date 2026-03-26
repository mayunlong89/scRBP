#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scRBP getSketch — Stratified GeoSketch cell downsampling

Applies geometric sketching (GeoSketch) to select a transcriptomically
representative subset of cells from large single-cell datasets.  Compared
with random downsampling, GeoSketch preserves rare cell populations and
transcriptional diversity in PCA space, yielding consistently improved
accuracy for downstream GRN inference.

For .h5ad input (with cell-type metadata), cells are sampled per cell type
so that each type contributes at least a user-specified minimum number of
cells, while keeping the total close to a global target (stratified mode).

For .feather input (no cell-type metadata), a global GeoSketch is applied
across all cells without per-type constraints.
"""

import argparse
import logging
import time
import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix, issparse
from sklearn.decomposition import TruncatedSVD
from geosketch import gs


# ------------------------------
# Helper: stratified GeoSketch for h5ad
# ------------------------------
def stratified_geosketch(adata, X_pca, celltype_col, n_total_target, min_cells_per_type, seed=42):
    """
    Perform stratified GeoSketch downsampling in the global PCA space.

    Each cell type contributes at least *min_cells_per_type* cells
    (or all available cells if fewer exist).  Any remaining quota up to
    *n_total_target* is distributed proportionally across cell types
    based on their remaining capacity.

    Parameters
    ----------
    adata : AnnData
        Input single-cell dataset.
    X_pca : np.ndarray
        PCA embedding of all cells, shape (n_cells, n_components).
    celltype_col : str
        Column in ``adata.obs`` containing cell-type labels.
    n_total_target : int
        Target total number of cells to retain.
    min_cells_per_type : int
        Minimum cells to retain per cell type.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    np.ndarray
        Array of selected global cell indices (sorted).
    """
    rng = np.random.RandomState(seed)

    celltypes = adata.obs[celltype_col].astype(str)
    unique_ct = celltypes.unique()

    logging.info(f"Number of cell types in '{celltype_col}': {len(unique_ct)}")
    total_cells = adata.n_obs
    n_total_target = min(n_total_target, total_cells)  # cannot exceed total cell count

    # Count cells per cell type
    ct_counts = {ct: int(np.sum(celltypes == ct)) for ct in unique_ct}
    for ct in unique_ct:
        logging.info(f"  Cell type '{ct}': total cells = {ct_counts[ct]}")

    # Step 1: guarantee at least min_cells_per_type per cell type
    base_n = {}
    for ct in unique_ct:
        n_avail = ct_counts[ct]
        base_n[ct] = min(n_avail, min_cells_per_type)

    base_total = sum(base_n.values())
    logging.info(f"Base allocation (min_cells_per_type={min_cells_per_type}) gives {base_total} cells.")

    # If the per-type minimum already exceeds the global target, accept the
    # larger total (per-type lower bound takes precedence over global target).
    if base_total >= n_total_target:
        logging.warning(
            f"Base_total = {base_total} >= target {n_total_target}. "
            f"Cannot satisfy both per-type minimum and target total. "
            f"Using base allocation only."
        )
        final_n = base_n
    else:
        # Step 2: distribute remaining quota proportionally by remaining capacity
        remaining = n_total_target - base_total
        logging.info(f"Remaining cells to distribute above base: {remaining}")

        extra_capacity = {ct: ct_counts[ct] - base_n[ct] for ct in unique_ct}
        total_extra_capacity = sum(max(v, 0) for v in extra_capacity.values())

        final_n = base_n.copy()
        if total_extra_capacity <= 0:
            logging.warning("No extra capacity in any cell type; using base allocation only.")
        else:
            # Allocate extra cells proportionally; cap at available capacity per type
            for ct in unique_ct:
                cap = max(extra_capacity[ct], 0)
                if cap == 0:
                    continue
                extra = int(round(remaining * cap / total_extra_capacity))
                extra = min(extra, cap)
                final_n[ct] += extra

            # Rounding may cause slight deviation from target; log for transparency
            allocated_total = sum(final_n.values())
            logging.info(
                f"After distributing extra cells, allocated_total = {allocated_total}, "
                f"target = {n_total_target} (may differ slightly due to rounding)."
            )

    # Step 3: run GeoSketch within each cell type
    selected_idx = []

    for ct in unique_ct:
        mask = (celltypes == ct).values
        idx_ct = np.where(mask)[0]
        X_ct = X_pca[mask, :]

        n_ct_target = final_n[ct]
        n_ct_target = min(n_ct_target, X_ct.shape[0])

        if n_ct_target <= 0:
            logging.info(f"Cell type '{ct}': n_ct_target=0, skipping.")
            continue

        logging.info(
            f"Cell type '{ct}': total={X_ct.shape[0]}, "
            f"sampling n={n_ct_target} cells via GeoSketch."
        )

        if X_ct.shape[0] == n_ct_target:
            # Retain all cells; no sketching needed
            selected_idx.extend(idx_ct.tolist())
        else:
            sketch_local = gs(X_ct, n_ct_target, replace=False)
            sketch_global = idx_ct[sketch_local]
            selected_idx.extend(sketch_global.tolist())

    selected_idx = np.array(sorted(set(selected_idx)))
    logging.info(f"Total sampled cells after stratified GeoSketch: {len(selected_idx)}")

    # Log final per-cell-type counts
    ct_sampled_counts = {}
    celltypes_sampled = celltypes.iloc[selected_idx]
    for ct in unique_ct:
        ct_sampled_counts[ct] = int(np.sum(celltypes_sampled == ct))
        logging.info(
            f"  FINAL: cell type '{ct}': sampled = {ct_sampled_counts[ct]} "
            f"(requested base+extra = {final_n[ct]}, total avail = {ct_counts[ct]})"
        )

    return selected_idx


def register_subcommand(subparsers):
    """Register this subcommand with the main parser."""
    parser = subparsers.add_parser(
        'getSketch',
        help='Stratified GeoSketch cell downsampling from h5ad / feather input',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )

    parser.add_argument("--input", type=str, required=True,
                        help="Input file: .h5ad or .feather (Seurat .rds not directly supported)")
    parser.add_argument("--output", type=str, required=True,
                        help="Output file: .h5ad / .csv / .feather / .npz")

    parser.add_argument("--n_cells", type=int, default=250000,
                        help="Target total number of cells to sketch (approximate)")
    parser.add_argument("--n_pca", type=int, default=100,
                        help="Number of PCA components for GeoSketch")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed")

    # h5ad-specific arguments (for stratified sampling)
    parser.add_argument("--celltype_col", type=str, default="celltype",
                        help="Column in adata.obs containing cell-type labels (for .h5ad input)")
    parser.add_argument("--min_cells_per_type", type=int, default=50,
                        help="Minimum number of cells per cell type (for .h5ad input)")

    parser.set_defaults(func=main)
    return parser


def main(args):
    """Main execution function."""
    # ------------------------------
    # Logging setup
    # ------------------------------
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler("geosketch_log.txt"), logging.StreamHandler()]
    )

    np.random.seed(args.seed)
    start_time = time.time()

    # Flag indicating whether stratified sampling was performed (h5ad input only)
    did_stratified_h5ad = False

    # ------------------------------
    # Detect input file type
    # ------------------------------
    input_ext = os.path.splitext(args.input)[-1].lower()
    logging.info(f"Detected file type: {input_ext}")

    # ------------------------------
    # Case 1: .h5ad input
    # ------------------------------
    if input_ext == ".h5ad":
        logging.info("Loading .h5ad file with scanpy...")
        adata = sc.read_h5ad(args.input)

        logging.info(f"Original shape (cells, genes): {adata.shape}")

        if args.celltype_col not in adata.obs.columns:
            raise ValueError(
                f"celltype_col '{args.celltype_col}' not found in adata.obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )

        logging.info("Performing PCA using TruncatedSVD on adata.X ...")
        svd = TruncatedSVD(n_components=args.n_pca, random_state=args.seed)
        X_pca = svd.fit_transform(adata.X)

        logging.info(f"Running stratified GeoSketch with target n_cells={args.n_cells}, "
                     f"min_cells_per_type={args.min_cells_per_type} ...")

        sketch_idx = stratified_geosketch(
            adata=adata,
            X_pca=X_pca,
            celltype_col=args.celltype_col,
            n_total_target=args.n_cells,
            min_cells_per_type=args.min_cells_per_type,
            seed=args.seed
        )

        # subset
        adata = adata[sketch_idx, :].copy()
        logging.info(f"Final subset shape (cells, genes): {adata.shape}")

        did_stratified_h5ad = True  # used below to export cell-to-celltype mapping

    # ------------------------------
    # Case 2: .feather input (no celltype info, global GeoSketch only)
    # ------------------------------
    elif input_ext == ".feather":
        logging.info("Loading feather expression matrix...")
        df = pd.read_feather(args.input)

        # Assume rows = genes, columns = cells (with an optional 'Gene' column for gene IDs)
        expr = df.set_index('Gene') if 'Gene' in df.columns else df
        X_sparse = csr_matrix(expr.values.T)  # shape: cells x genes

        logging.info("Performing PCA using TruncatedSVD...")
        svd = TruncatedSVD(n_components=args.n_pca, random_state=args.seed)
        X_pca = svd.fit_transform(X_sparse)

        logging.warning(
            "Feather input has no cell-type metadata; using global GeoSketch only. "
            "Cell-type-wise minimum sampling cannot be enforced."
        )

        logging.info(f"Running GeoSketch on {X_pca.shape[0]} cells...")
        n_cells = min(args.n_cells, X_pca.shape[0])
        sketch_idx = gs(X_pca, n_cells, replace=False)

        # Subset and reformat (genes x selected_cells)
        sketch_df = expr.iloc[:, sketch_idx]
        sketch_df.reset_index(inplace=True)
        sketch_df.rename(columns={"index": "Gene"}, inplace=True)

        # Convert to AnnData for unified saving logic (cells x genes inside AnnData)
        adata = sc.AnnData(X=sketch_df.set_index("Gene").T.values)
        adata.var_names = sketch_df["Gene"].values
        adata.obs_names = sketch_df.columns[1:]  # skip "Gene"
        logging.info(f"Converted feather to AnnData object; shape = {adata.shape}")

    # ------------------------------
    # Case 3: .rds input (unsupported)
    # ------------------------------
    elif input_ext == ".rds":
        logging.error("Direct .rds (Seurat) input is not supported in Python.")
        raise ValueError("Please convert .rds to .h5ad using SeuratDisk in R.")

    # ------------------------------
    # Unsupported format
    # ------------------------------
    else:
        logging.error("Unsupported input file format. Use .h5ad or .feather.")
        raise ValueError("Unsupported input file format.")

    # ------------------------------
    # For h5ad input with stratified sampling, export a cell-to-celltype mapping CSV
    # ------------------------------
    if did_stratified_h5ad:
        celltype_col = args.celltype_col
        mapping_df = (
            adata.obs[[celltype_col]]
            .copy()
            .reset_index()
            .rename(columns={"index": "cell", celltype_col: "cell_type"})
        )
        mapping_out = os.path.splitext(args.output)[0] + "_cell_to_celltype.csv"
        mapping_df.to_csv(mapping_out, index=False)
        logging.info(f"Saved cell-to-celltype mapping to {mapping_out}")

    # ------------------------------
    # Save in format determined by output suffix
    # ------------------------------
    output_ext = os.path.splitext(args.output)[-1].lower()

    if output_ext == ".h5ad":
        logging.info(f"Saving sketch to {args.output} (AnnData: cells x genes) ...")
        adata.write(args.output)

    elif output_ext == ".csv":
        logging.info(f"Saving sketch to CSV (gene-by-cell): {args.output}")
        if isinstance(adata.X, np.ndarray):
            dense_matrix = adata.X
        else:
            dense_matrix = adata.X.toarray()
        df_out = pd.DataFrame(dense_matrix.T, index=adata.var_names, columns=adata.obs_names)
        df_out = df_out.reset_index().rename(columns={"index": "Gene"})
        df_out.to_csv(args.output, index=False)

    elif output_ext == ".feather":
        logging.info(f"Saving sketch to Feather (gene-by-cell): {args.output}")
        if isinstance(adata.X, np.ndarray):
            dense_matrix = adata.X
        else:
            dense_matrix = adata.X.toarray()
        df_out = pd.DataFrame(dense_matrix.T, index=adata.var_names, columns=adata.obs_names)
        df_out = df_out.reset_index().rename(columns={"index": "Gene"})
        df_out.to_feather(args.output)

    elif output_ext == ".npz":
        logging.info(f"Saving sketch matrix to sparse .npz (gene-by-cell): {args.output}")
        from scipy.sparse import save_npz
        if issparse(adata.X):
            mat = adata.X.T
        else:
            mat = csr_matrix(adata.X.T)
        save_npz(args.output, mat)

    else:
        logging.error("Unsupported output format. Use .h5ad, .csv, .feather, or .npz")
        raise ValueError("Unsupported output format")

    logging.info(f"Sketching complete. Total time: {time.time() - start_time:.2f} seconds.")


if __name__ == "__main__":
    import sys
    parser = argparse.ArgumentParser(description="Stratified GeoSketch downsampling of cells from h5ad / feather input.")
    register_subcommand(parser.add_subparsers()).set_defaults(func=main)
    args = parser.parse_args()
    main(args)
