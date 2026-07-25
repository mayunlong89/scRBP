#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scRBP getMetacell — Aggregate single cells into mini-metacells for robust GRN inference

Single-cell RNA-seq is sparse (high dropout), which weakens the co-expression
signal that GRNBoost2 relies on — an effect that is *more* pronounced for RBP
regulators than for transcription factors, because post-transcriptional
co-regulation typically has a smaller dynamic range.  Following the
mini-metacell strategy used in large regulon atlases, ``getMetacell`` pools
groups of ~10-15 transcriptomically similar cells into denser pseudo-profiles,
markedly improving signal-to-noise for downstream ``getGRN``.

Relationship to ``getSketch`` (they are COMPLEMENTARY, not alternatives)
------------------------------------------------------------------------
  - ``getSketch``  *selects* a diversity-preserving subset of REAL single cells
                   (solves the "too many cells" / compute problem; keeps cells
                   sparse).  Use it to feed the ACTIVITY branch (``ras``), which
                   needs real single cells to keep cell-type resolution.
  - ``getMetacell`` *aggregates* similar cells into DENSE units (solves the
                   "too sparse" / signal problem).  Use it to feed the GRN
                   branch (``getGRN`` -> ``getModule`` -> ``getPrune``).

Recommended composition for a large multi-tissue atlas:

    per lineage/organ:
      [GRN branch]  (optional getSketch pre-thin) -> getMetacell --within-celltype -> getGRN ...
      [RAS branch]  getSketch (single cells, keep cell types)      -> ras -> rss -> trs

IMPORTANT — metacells are built WITHIN each cell type by default, so that the
diversity ``getSketch`` protects (and cell-type boundaries) are never blurred
by pooling across heterogeneous cells.

Normalisation contract
-----------------------
Cell *similarity* (the neighbour graph / clustering used to decide which cells
to pool) is computed on a library-size-normalised, log1p-transformed PCA
embedding, while the metacell *profile* is formed by aggregating the ORIGINAL
input counts (sum by default).  Therefore feed ``getMetacell`` raw/linear
counts, exactly as you would feed ``getGRN`` afterwards.

Grouping methods (``--method``)
-------------------------------
  knn     Greedy non-overlapping nearest-neighbour pooling in PCA space.
          Produces size-controlled, transcriptomically coherent metacells.
          (default; most faithful to the metacell concept.)
  kmeans  MiniBatchKMeans with k = round(n_cells / metacell_size) clusters.
          Scales well; cluster sizes vary around the target.
  random  Random partition into fixed-size bins (no embedding needed).
          Fastest; ignores similarity — use only as a cheap baseline / for
          already-homogeneous inputs.

Inputs / outputs
----------------
  Input  : .h5ad (with a cell-type column for stratified pooling) or
           .feather (gene-by-cell; no cell-type metadata -> global pooling).
  Output : metacell expression matrix (.h5ad / .csv / .feather / .npz),
           gene-by-metacell, directly consumable by ``getGRN``.
           Plus  <out>_metacell_to_celltype.csv  (cell, cell_type, n_cells)
           and   <out>_metacell_summary.csv      (per-cell-type size stats).
"""

import argparse
import logging
import os
import time
from math import ceil

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix, diags, issparse, save_npz
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors


# ──────────────────────────────────────────────────────────────
# Input loading
# ──────────────────────────────────────────────────────────────
def load_input(input_path, celltype_col):
    """
    Load expression as a (cells x genes) CSR matrix plus identifiers.

    Returns
    -------
    X_raw      : scipy.sparse.csr_matrix, shape (n_cells, n_genes) — ORIGINAL counts
    gene_names : np.ndarray[str]
    cell_names : np.ndarray[str]
    celltypes  : pd.Series indexed by cell_names, or None if unavailable
    """
    ext = os.path.splitext(input_path)[-1].lower()

    if ext == ".h5ad":
        logging.info("Loading .h5ad with scanpy ...")
        adata = sc.read_h5ad(input_path)
        logging.info(f"Original shape (cells, genes): {adata.shape}")
        X_raw = adata.X
        X_raw = X_raw.tocsr() if issparse(X_raw) else csr_matrix(X_raw)
        gene_names = np.asarray(adata.var_names, dtype=object)
        cell_names = np.asarray(adata.obs_names, dtype=object)
        celltypes = None
        if celltype_col in adata.obs.columns:
            celltypes = pd.Series(
                adata.obs[celltype_col].astype(str).values, index=cell_names, name="cell_type"
            )
        else:
            logging.warning(
                f"celltype_col '{celltype_col}' not in adata.obs; falling back to GLOBAL pooling. "
                f"Available: {list(adata.obs.columns)}"
            )
        return X_raw, gene_names, cell_names, celltypes

    elif ext == ".feather":
        logging.info("Loading .feather (gene-by-cell) ...")
        df = pd.read_feather(input_path)
        idx_col = "Gene" if "Gene" in df.columns else df.columns[0]
        expr = df.set_index(idx_col)            # genes x cells
        gene_names = np.asarray(expr.index, dtype=object)
        cell_names = np.asarray(expr.columns, dtype=object)
        X_raw = csr_matrix(expr.values.T)       # cells x genes
        logging.warning(
            "Feather input has no cell-type metadata; using GLOBAL pooling "
            "(metacells may mix cell types). Provide an .h5ad with --celltype_col "
            "for within-cell-type pooling."
        )
        return X_raw, gene_names, cell_names, None

    elif ext == ".rds":
        raise ValueError("Direct .rds (Seurat) input is not supported. Convert to .h5ad with SeuratDisk.")
    else:
        raise ValueError(f"Unsupported input format: {ext}. Use .h5ad or .feather.")


# ──────────────────────────────────────────────────────────────
# Embedding: normalise for similarity, then PCA (similarity only;
# aggregation always uses the ORIGINAL counts)
# ──────────────────────────────────────────────────────────────
def compute_embedding(X_raw, n_pca, seed):
    """
    Library-size normalise (to median total) + log1p, then TruncatedSVD.

    The embedding is used ONLY to decide which cells are similar; metacell
    profiles are aggregated from the original counts in `X_raw`.
    """
    n_cells, n_genes = X_raw.shape
    lib = np.asarray(X_raw.sum(axis=1)).ravel().astype(np.float64)
    med = np.median(lib[lib > 0]) if np.any(lib > 0) else 1.0
    with np.errstate(divide="ignore", invalid="ignore"):
        scale = np.where(lib > 0, med / lib, 0.0)

    Xn = diags(scale) @ X_raw            # row-scale (sparse-safe)
    Xn = Xn.tocsr()
    Xn.data = np.log1p(Xn.data)

    n_comp = int(min(n_pca, n_genes - 1, n_cells - 1))
    n_comp = max(2, n_comp)
    logging.info(f"Computing PCA embedding via TruncatedSVD (n_components={n_comp}) ...")
    svd = TruncatedSVD(n_components=n_comp, random_state=seed)
    Z = svd.fit_transform(Xn)
    return np.ascontiguousarray(Z, dtype=np.float32)


# ──────────────────────────────────────────────────────────────
# Grouping strategies (operate on one cell-type block; return local labels 0..g-1)
# ──────────────────────────────────────────────────────────────
def _group_random(n, size, rng):
    order = rng.permutation(n)
    n_groups = max(1, ceil(n / size))
    labels = np.empty(n, dtype=np.int64)
    for g in range(n_groups):
        labels[order[g * size:(g + 1) * size]] = g
    return labels


def _group_kmeans(Z, size, seed):
    n = Z.shape[0]
    k = max(1, int(round(n / size)))
    if k == 1:
        return np.zeros(n, dtype=np.int64)
    km = MiniBatchKMeans(n_clusters=k, random_state=seed, n_init=3, batch_size=max(256, 10 * size))
    return km.fit_predict(Z).astype(np.int64)


def _group_knn_greedy(Z, size, seed):
    """Greedy non-overlapping nearest-neighbour pooling. Metacells of up to `size` cells."""
    n = Z.shape[0]
    if n <= size:
        return np.zeros(n, dtype=np.int64)
    n_neighbors = int(min(size, n))
    nn = NearestNeighbors(n_neighbors=n_neighbors).fit(Z)
    _, idx = nn.kneighbors(Z)                 # (n, n_neighbors), nearest-first
    assigned = np.full(n, -1, dtype=np.int64)
    rng = np.random.default_rng(seed)
    mc = 0
    for i in rng.permutation(n):
        if assigned[i] != -1:
            continue
        assigned[i] = mc
        filled = 1
        for j in idx[i]:
            if filled >= size:
                break
            if assigned[j] == -1:
                assigned[j] = mc
                filled += 1
        mc += 1
    return assigned


def _merge_small(labels, Z, min_size):
    """Merge metacells smaller than `min_size` into the nearest remaining metacell.

    Distance = Euclidean between metacell centroids in Z. If Z is None
    (random method) small metacells are merged into the largest one. Labels
    are relabelled to a contiguous 0..M-1 range. Returns new labels.
    """
    if min_size <= 1:
        return labels
    labels = labels.copy()
    counts = np.bincount(labels)
    small = np.where(counts < min_size)[0]
    large = np.where(counts >= min_size)[0]
    if small.size == 0:
        return labels
    if large.size == 0:
        # everything is small (tiny cell type): collapse to a single metacell
        return np.zeros_like(labels)

    if Z is not None:
        cent = np.vstack([Z[labels == g].mean(axis=0) for g in range(counts.size)])
        large_cent = cent[large]
        for g in small:
            d = np.sum((large_cent - cent[g]) ** 2, axis=1)
            labels[labels == g] = large[int(np.argmin(d))]
    else:
        target = large[int(np.argmax(counts[large]))]
        for g in small:
            labels[labels == g] = target

    # relabel to contiguous ids
    uniq = np.unique(labels)
    remap = {old: new for new, old in enumerate(uniq)}
    return np.array([remap[v] for v in labels], dtype=np.int64)


# ──────────────────────────────────────────────────────────────
# Aggregation
# ──────────────────────────────────────────────────────────────
def aggregate(X_raw, global_labels, n_metacells, agg):
    """Aggregate ORIGINAL counts into (n_metacells x genes) using a sparse pooling matrix."""
    n_cells = X_raw.shape[0]
    sizes = np.bincount(global_labels, minlength=n_metacells).astype(np.float64)
    if agg == "mean":
        with np.errstate(divide="ignore"):
            w = np.where(sizes[global_labels] > 0, 1.0 / sizes[global_labels], 0.0)
    else:  # sum
        w = np.ones(n_cells, dtype=np.float64)
    G = csr_matrix((w, (global_labels, np.arange(n_cells))), shape=(n_metacells, n_cells))
    MC = G @ X_raw                       # (n_metacells x genes), sparse
    return MC.tocsr(), sizes.astype(int)


# ──────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────
def register_subcommand(subparsers):
    """Register this subcommand with the main parser."""
    parser = subparsers.add_parser(
        "getMetacell",
        help="Aggregate single cells into mini-metacells for robust GRN inference",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )

    parser.add_argument("--input", required=True,
                        help="Input file: .h5ad (with cell-type column) or .feather (gene-by-cell).")
    parser.add_argument("--output", required=True,
                        help="Output file: .h5ad / .csv / .feather / .npz (gene-by-metacell).")

    parser.add_argument("--metacell_size", type=int, default=10,
                        help="Target number of cells pooled per metacell (default: 10; atlas typical 10-15).")
    parser.add_argument("--method", choices=["knn", "kmeans", "random"], default="knn",
                        help="Cell-grouping strategy (default: knn).")
    parser.add_argument("--agg", choices=["sum", "mean"], default="sum",
                        help="Aggregate pooled cells by sum (default) or mean of ORIGINAL counts.")

    parser.add_argument("--within_celltype", dest="within_celltype", action="store_true", default=True,
                        help="Pool only within each cell type (default: on; requires --celltype_col in .h5ad).")
    parser.add_argument("--global_pooling", dest="within_celltype", action="store_false",
                        help="Disable within-cell-type pooling; pool across all cells.")
    parser.add_argument("--celltype_col", type=str, default="celltype",
                        help="Column in adata.obs holding cell-type labels (.h5ad input).")

    parser.add_argument("--min_metacell_size", type=int, default=1,
                        help="Merge metacells smaller than this into the nearest one "
                             "(default: 1 = keep all). Recommended ~metacell_size//2 for cleaner pooling.")
    parser.add_argument("--n_pca", type=int, default=50,
                        help="PCA components for the similarity embedding (knn/kmeans).")
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument("--save_members", action="store_true",
                        help="Also write <out>_metacell_members.csv mapping each metacell to its source cells.")

    parser.set_defaults(func=main)
    return parser


def main(args):
    """Main execution function."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler("getMetacell_log.txt"), logging.StreamHandler()],
    )
    t0 = time.time()
    rng = np.random.default_rng(args.seed)

    if args.metacell_size < 2:
        raise ValueError("--metacell_size must be >= 2 (a metacell pools multiple cells).")

    # 1) Load original counts
    X_raw, gene_names, cell_names, celltypes = load_input(args.input, args.celltype_col)
    n_cells, n_genes = X_raw.shape
    logging.info(f"Loaded {n_cells} cells x {n_genes} genes.")

    # 2) Decide grouping blocks (cell types or one global block)
    use_ct = (celltypes is not None) and args.within_celltype
    if use_ct:
        ct_values = celltypes.reindex(cell_names).astype(str).values
        blocks = [(ct, np.where(ct_values == ct)[0]) for ct in pd.unique(ct_values)]
        logging.info(f"Within-cell-type pooling across {len(blocks)} cell types.")
    else:
        ct_values = np.array(["all"] * n_cells, dtype=object)
        blocks = [("all", np.arange(n_cells))]
        logging.info("Global pooling (single block).")

    # 3) Embedding (only for knn/kmeans)
    Z = None
    if args.method in ("knn", "kmeans"):
        Z = compute_embedding(X_raw, args.n_pca, args.seed)

    # 4) Assign global metacell labels block by block
    global_labels = np.full(n_cells, -1, dtype=np.int64)
    metacell_celltype = []
    next_id = 0
    for ct, idx in blocks:
        if idx.size == 0:
            continue
        if args.method == "random":
            local = _group_random(idx.size, args.metacell_size, rng)
        elif args.method == "kmeans":
            local = _group_kmeans(Z[idx], args.metacell_size, args.seed)
        else:  # knn
            local = _group_knn_greedy(Z[idx], args.metacell_size, args.seed)

        if args.min_metacell_size > 1:
            local = _merge_small(local, None if Z is None else Z[idx], args.min_metacell_size)

        n_local = int(local.max()) + 1
        global_labels[idx] = local + next_id
        metacell_celltype.extend([ct] * n_local)
        next_id += n_local
        logging.info(f"  cell type '{ct}': {idx.size} cells -> {n_local} metacells")

    n_metacells = next_id
    if (global_labels < 0).any():
        raise RuntimeError("Some cells were not assigned to a metacell (internal error).")
    logging.info(f"Total: {n_cells} cells -> {n_metacells} metacells "
                 f"(median {np.median(np.bincount(global_labels)):.1f} cells/metacell).")

    # 5) Aggregate ORIGINAL counts
    MC, sizes = aggregate(X_raw, global_labels, n_metacells, args.agg)
    metacell_ids = np.array([f"MC{i:07d}" for i in range(n_metacells)], dtype=object)
    metacell_celltype = np.array(metacell_celltype, dtype=object)

    # 6) Write metacell matrix (gene-by-metacell), reusing getSketch's output conventions
    _write_matrix(args.output, MC, gene_names, metacell_ids)

    # 7) Write metacell -> cell-type mapping (drop-in for `ras --celltypes-csv`)
    base = os.path.splitext(args.output)[0]
    map_out = f"{base}_metacell_to_celltype.csv"
    pd.DataFrame({"cell": metacell_ids, "cell_type": metacell_celltype, "n_cells": sizes}).to_csv(
        map_out, index=False
    )
    logging.info(f"Saved metacell->cell-type mapping: {map_out}")

    # 8) Per-cell-type size summary
    summ = (pd.DataFrame({"cell_type": metacell_celltype, "n_cells": sizes})
            .groupby("cell_type")["n_cells"]
            .agg(n_metacells="size", median_cells="median",
                 min_cells="min", max_cells="max", total_cells="sum")
            .reset_index())
    summ_out = f"{base}_metacell_summary.csv"
    summ.to_csv(summ_out, index=False)
    logging.info(f"Saved size summary: {summ_out}")

    # 9) Optional: full membership table
    if args.save_members:
        mem_out = f"{base}_metacell_members.csv"
        pd.DataFrame({
            "cell": cell_names,
            "metacell": metacell_ids[global_labels],
            "cell_type": ct_values,
        }).to_csv(mem_out, index=False)
        logging.info(f"Saved membership table: {mem_out}")

    logging.info(f"getMetacell complete in {time.time() - t0:.2f}s. "
                 f"Feed '{args.output}' to: scRBP getGRN --matrix {args.output} ...")


def _write_matrix(output, MC_csr, gene_names, metacell_ids):
    """Write a (metacells x genes) CSR as gene-by-metacell, matching getSketch output formats."""
    ext = os.path.splitext(output)[-1].lower()
    parent = os.path.dirname(output)
    if parent:
        os.makedirs(parent, exist_ok=True)

    if ext == ".npz":
        logging.info(f"Saving sparse .npz (gene-by-metacell): {output}")
        save_npz(output, MC_csr.T.tocsr())     # genes x metacells
        return

    if ext == ".h5ad":
        logging.info(f"Saving AnnData (metacells x genes): {output}")
        adata = sc.AnnData(X=MC_csr)            # cells(metacells) x genes
        adata.var_names = gene_names.astype(str)
        adata.obs_names = metacell_ids.astype(str)
        adata.write(output)
        return

    # csv / feather -> dense gene-by-metacell with leading "Gene" column
    dense = MC_csr.toarray().T                  # genes x metacells
    df_out = pd.DataFrame(dense, index=gene_names, columns=metacell_ids)
    df_out = df_out.reset_index().rename(columns={"index": "Gene"})
    if ext == ".csv":
        logging.info(f"Saving CSV (gene-by-metacell): {output}")
        df_out.to_csv(output, index=False)
    elif ext == ".feather":
        logging.info(f"Saving Feather (gene-by-metacell): {output}")
        df_out.columns = df_out.columns.astype(str)
        df_out.to_feather(output)
    else:
        raise ValueError(f"Unsupported output format: {ext}. Use .h5ad, .csv, .feather, or .npz.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Aggregate single cells into mini-metacells for robust GRN inference."
    )
    register_subcommand(parser.add_subparsers())
    args = parser.parse_args()
    main(args)
