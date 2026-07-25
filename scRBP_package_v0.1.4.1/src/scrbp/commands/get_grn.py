#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scRBP getGRN - Gene/Isoform regulatory network inference with GRNBoost2/GENIE3

Supports two inference modes:
  --mode gene     : RBP → Gene network (classic mode)
  --mode isoform  : RBP → Isoform network (requires --isoform_annotation)

In gene mode, RBPs present in the expression matrix are used as regulators and
all genes are treated as targets.

In isoform mode, multiple isoforms of the same RBP gene are aggregated into a
single regulator signal (sum/mean/max), and individual isoforms serve as targets
after expression-level filtering.

Correlation is computed as Spearman (vectorised rank-based Pearson on ranks).
"""

from __future__ import annotations

import os
import sys
import time
import gzip
import logging
import argparse
from typing import Optional, Set, Dict, List

import numpy as np
import pandas as pd
import polars as pl
import tqdm

from multiprocessing import Pool, cpu_count
from concurrent.futures import ProcessPoolExecutor, as_completed

from arboreto.utils import load_tf_names
from arboreto.core import (
    SGBM_KWARGS,
    RF_KWARGS,
    EARLY_STOP_WINDOW_LENGTH,
    infer_partial_network,
)

from scrbp.utils.stats import assign_mode

# ──────────────────────────────────────────────────────────────
# Logging
# ──────────────────────────────────────────────────────────────
LOGGER = logging.getLogger("scRBP_Infer_GRN_Pipeline")


def create_logging_handler(debug: bool, log_file: Optional[str]) -> logging.Handler:
    handler = (
        logging.FileHandler(log_file)
        if log_file
        else logging.StreamHandler(stream=sys.stderr)
    )
    handler.setLevel(logging.DEBUG if debug else logging.INFO)
    handler.setFormatter(
        logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    )
    return handler


# ──────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────

def str2bool(v: str) -> bool:
    if isinstance(v, bool):
        return v
    val = str(v).strip().lower()
    if val in {"true", "1", "yes", "y"}:
        return True
    if val in {"false", "0", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Boolean value expected, got: {v}")


def check_file_exists(path: str, file_desc: str):
    if path is None:
        raise ValueError(f"{file_desc} is required but got None.")
    if not os.path.exists(path):
        raise FileNotFoundError(f"{file_desc} not found: {path}")


def strip_tsv_suffix(path: str) -> str:
    return path[:-4] if path.endswith(".tsv") else path


def build_main_output_path(output_arg: str, mode: str) -> str:
    base = strip_tsv_suffix(output_arg)
    if mode == "gene":
        return f"{base}_scRBP_gene_GRNs.tsv"
    return f"{base}_scRBP_isoform_GRNs.tsv"


def build_aux_output_path(output_arg: str, suffix: str) -> str:
    base = strip_tsv_suffix(output_arg)
    return f"{base}_{suffix}"


# ──────────────────────────────────────────────────────────────
# Expression matrix loader  (gene mode: keep_features=RBP set)
# ──────────────────────────────────────────────────────────────

def load_expression_matrix(
    matrix_path: str,
    keep_features: Optional[Set[str]] = None,
    dtype: str = "float32",
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Load expression matrix and return cells × features DataFrame.

    Input files are assumed to have rows = features (genes or isoforms)
    and columns = cells; the function transposes automatically.

    Supported formats: .csv  .csv.gz  .feather  .loom
    """
    def _log(msg):
        if verbose:
            print(msg)

    if not os.path.exists(matrix_path):
        raise FileNotFoundError(f"Matrix file not found: {matrix_path}")
    if not matrix_path.lower().endswith((".csv", ".csv.gz", ".feather", ".loom")):
        raise ValueError(f"Matrix must be .csv/.csv.gz, .feather, or .loom, got: {matrix_path}")

    LOGGER.info(f"Loading expression matrix from {matrix_path} ...")
    lower = matrix_path.lower()

    if lower.endswith((".csv", ".csv.gz")):
        _log(f"[Matrix] Loading CSV via Polars: {matrix_path}")
        if lower.endswith(".gz"):
            with gzip.open(matrix_path, "rt") as fh:
                df_pl = pl.read_csv(fh, has_header=True)
        else:
            df_pl = pl.read_csv(matrix_path, has_header=True)
        index_col = df_pl.columns[0]
        if keep_features is not None:
            df_pl = df_pl.filter(pl.col(index_col).is_in(list(keep_features)))
        df = df_pl.to_pandas()
        df.set_index(index_col, inplace=True)
        df = df.T

    elif lower.endswith(".feather"):
        _log(f"[Matrix] Loading Feather via pandas: {matrix_path}")
        df = pd.read_feather(matrix_path)
        index_col = df.columns[0]
        df.set_index(index_col, inplace=True)
        if keep_features is not None:
            df = df.loc[df.index.intersection(keep_features)]
        df = df.T

    elif lower.endswith(".loom"):
        _log(f"[Matrix] Loading Loom via loompy: {matrix_path}")
        try:
            import loompy
        except ImportError as exc:
            raise ImportError(
                "Please install loompy to read .loom files: pip install loompy"
            ) from exc
        ds = loompy.connect(matrix_path, mode="r")
        feature_keys = ("Gene", "gene", "gene_symbol", "gene_names", "var_names", "Genes")
        feature_names = next(
            (ds.row_attrs[k] for k in feature_keys if k in ds.row_attrs),
            [f"feature_{i}" for i in range(ds.shape[0])],
        )
        cell_keys = ("CellID", "cell", "cell_names", "obs_names", "Cells")
        cell_names = next(
            (ds.col_attrs[k] for k in cell_keys if k in ds.col_attrs),
            [f"cell_{i}" for i in range(ds.shape[1])],
        )
        if keep_features is not None:
            f2i = {f: i for i, f in enumerate(feature_names)}
            sel = sorted(f2i[f] for f in keep_features if f in f2i)
            if not sel:
                ds.close()
                raise ValueError("None of keep_features found in loom row attrs.")
            matrix = ds[sel, :]
            cols = [feature_names[i] for i in sel]
        else:
            matrix = ds[:, :]
            cols = list(feature_names)
        ds.close()
        df = pd.DataFrame(matrix.T, index=cell_names, columns=cols)
    else:
        raise ValueError(f"Unsupported matrix format: {matrix_path}")

    df = df.loc[~df.index.duplicated(keep="first")]
    df = df.loc[:, ~pd.Index(df.columns).duplicated(keep="first")]
    if keep_features is not None:
        df = df.loc[:, df.columns.intersection(keep_features)]
    if dtype:
        df = df.astype(dtype, copy=False)

    LOGGER.info(f"Expression matrix loaded: {df.shape[0]} cells × {df.shape[1]} features")
    return df


# ──────────────────────────────────────────────────────────────
# Isoform annotation loader
# ──────────────────────────────────────────────────────────────

def load_isoform_annotation(annotation_path: str) -> pd.DataFrame:
    """
    Load isoform→gene annotation table.

    Accepted column aliases:
      isoform side : isoform, isoform_id, transcript, transcript_id, tx, tx_id
      gene side    : gene, gene_symbol, symbol, gene_name

    Returns a DataFrame with columns ["Isoform", "Gene"].
    Raises ValueError if any isoform maps to more than one gene.
    """
    check_file_exists(annotation_path, "Isoform annotation file")
    LOGGER.info(f"Loading isoform annotation from {annotation_path} ...")

    sep = "\t" if annotation_path.lower().endswith((".tsv", ".txt")) else ","
    anno = pd.read_csv(annotation_path, sep=sep)
    colmap = {c.lower(): c for c in anno.columns}

    isoform_col = next(
        (colmap[c] for c in ["isoform", "isoform_id", "transcript", "transcript_id", "tx", "tx_id"] if c in colmap),
        None,
    )
    gene_col = next(
        (colmap[c] for c in ["gene", "gene_symbol", "symbol", "gene_name"] if c in colmap),
        None,
    )
    if isoform_col is None or gene_col is None:
        raise ValueError(
            "Could not detect required columns in isoform annotation file. "
            "Need one isoform/transcript column and one gene/gene_symbol column. "
            f"Available columns: {list(anno.columns)}"
        )

    anno = anno[[isoform_col, gene_col]].copy()
    anno.columns = ["Isoform", "Gene"]
    anno = anno.astype(str).dropna().drop_duplicates()

    dup_check = anno.groupby("Isoform")["Gene"].nunique()
    bad = dup_check[dup_check > 1]
    if len(bad) > 0:
        raise ValueError(
            f"Found isoforms mapping to multiple genes. Example isoforms: {list(bad.index[:10])}"
        )

    LOGGER.info(f"Isoform annotation loaded: {anno.shape[0]} isoform-gene pairs")
    return anno


# ──────────────────────────────────────────────────────────────
# RBP aggregation (isoform mode)
# ──────────────────────────────────────────────────────────────

def aggregate_rbp_expression_from_isoforms(
    expr_cells_x_isoforms: pd.DataFrame,
    rbp_genes: List[str],
    isoform_annotation: pd.DataFrame,
    agg_method: str = "sum",
):
    """
    Aggregate isoform-level expression into gene-level RBP regulators.

    Returns
    -------
    rbp_expr_df          : cells × matched_RBPs aggregated matrix
    rbp_isoform_map_df   : mapping table {RBP, Isoform, n_isoforms_for_RBP, agg_method}
    annotation_summary_df: summary statistics DataFrame
    gene_to_isoforms     : dict {RBP_gene: [isoform_ids]}
    """
    LOGGER.info(f"Aggregating RBP expression from isoforms using method: {agg_method}")
    anno = isoform_annotation.copy()

    n_matrix_isoforms = expr_cells_x_isoforms.shape[1]
    n_input_rbps = len(rbp_genes)

    anno_in_matrix = anno[anno["Isoform"].isin(expr_cells_x_isoforms.columns)].copy()
    anno_rbp = anno_in_matrix[anno_in_matrix["Gene"].isin(rbp_genes)].copy()

    if anno_rbp.empty:
        raise ValueError("No annotated RBP isoforms found in expression matrix.")

    gene_to_isoforms: Dict[str, List[str]] = (
        anno_rbp.groupby("Gene")["Isoform"]
        .apply(lambda x: list(pd.unique(x)))
        .to_dict()
    )

    matched_rbps = sorted(gene_to_isoforms.keys())
    if not matched_rbps:
        raise ValueError("No RBPs could be matched to isoforms via annotation.")

    agg_frames, map_rows = [], []
    for rbp in matched_rbps:
        isoforms = [x for x in gene_to_isoforms[rbp] if x in expr_cells_x_isoforms.columns]
        if not isoforms:
            continue
        sub = expr_cells_x_isoforms[isoforms]
        if agg_method == "sum":
            agg_vec = sub.sum(axis=1)
        elif agg_method == "mean":
            agg_vec = sub.mean(axis=1)
        elif agg_method == "max":
            agg_vec = sub.max(axis=1)
        else:
            raise ValueError(f"Unsupported agg_method: {agg_method}")
        agg_frames.append(agg_vec.rename(rbp))
        for iso in isoforms:
            map_rows.append({"RBP": rbp, "Isoform": iso,
                             "n_isoforms_for_RBP": len(isoforms), "agg_method": agg_method})

    if not agg_frames:
        raise ValueError("No RBP expression vectors could be aggregated from isoforms.")

    rbp_expr_df = pd.concat(agg_frames, axis=1)
    rbp_expr_df = rbp_expr_df.loc[:, ~rbp_expr_df.columns.duplicated(keep="first")]

    rbp_isoform_map_df = (
        pd.DataFrame(map_rows).sort_values(["RBP", "Isoform"]).reset_index(drop=True)
    )

    isoforms_per_rbp = rbp_isoform_map_df.groupby("RBP")["Isoform"].nunique()
    annotation_summary_df = pd.DataFrame([
        {"metric": "input_isoforms_in_matrix",            "value": int(n_matrix_isoforms)},
        {"metric": "input_rbps_in_list",                  "value": int(n_input_rbps)},
        {"metric": "annotation_pairs_total",              "value": int(isoform_annotation.shape[0])},
        {"metric": "annotation_isoforms_matched_to_matrix","value": int(anno_in_matrix["Isoform"].nunique())},
        {"metric": "matched_rbp_isoform_pairs",           "value": int(anno_rbp.shape[0])},
        {"metric": "matched_rbps",                        "value": int(len(matched_rbps))},
        {"metric": "aggregated_rbp_matrix_n_cells",       "value": int(rbp_expr_df.shape[0])},
        {"metric": "aggregated_rbp_matrix_n_rbps",        "value": int(rbp_expr_df.shape[1])},
        {"metric": "median_isoforms_per_rbp",             "value": float(np.median(isoforms_per_rbp.values))},
        {"metric": "max_isoforms_per_rbp",                "value": int(np.max(isoforms_per_rbp.values))},
        {"metric": "agg_method",                          "value": agg_method},
    ])

    LOGGER.info(f"Aggregated RBP matrix: {rbp_expr_df.shape[0]} cells × {rbp_expr_df.shape[1]} RBPs")
    return rbp_expr_df, rbp_isoform_map_df, annotation_summary_df, gene_to_isoforms


# ──────────────────────────────────────────────────────────────
# Target isoform filtering
# ──────────────────────────────────────────────────────────────

def filter_target_isoforms(
    expr_cells_x_isoforms: pd.DataFrame,
    min_cells_expressed: int = 10,
    min_mean_expr: float = 0.01,
):
    """Remove low-expressed target isoforms before inference."""
    expr_nonzero_cells = (expr_cells_x_isoforms > 0).sum(axis=0)
    expr_mean = expr_cells_x_isoforms.mean(axis=0)
    keep_mask = (expr_nonzero_cells >= min_cells_expressed) & (expr_mean >= min_mean_expr)
    kept = list(expr_cells_x_isoforms.columns[keep_mask])

    summary_df = pd.DataFrame({
        "Isoform":           expr_cells_x_isoforms.columns,
        "n_cells_expressed": expr_nonzero_cells.values,
        "mean_expr":         expr_mean.values,
        "kept_as_target":    keep_mask.values,
    })
    LOGGER.info(
        f"Target isoform filtering: kept {len(kept)} / {expr_cells_x_isoforms.shape[1]} "
        f"(min_cells_expressed>={min_cells_expressed}, min_mean_expr>={min_mean_expr})"
    )
    return kept, summary_df


# ──────────────────────────────────────────────────────────────
# Network inference worker
# ──────────────────────────────────────────────────────────────

def run_network_for_target(task):
    target_name, target_expr, regulator_matrix, regulator_names, method, seed = task
    try:
        regressor_type = "GBM" if method == "grnboost2" else "RF"
        regressor_kwargs = SGBM_KWARGS if method == "grnboost2" else RF_KWARGS
        result = infer_partial_network(
            regressor_type=regressor_type,
            regressor_kwargs=regressor_kwargs,
            tf_matrix=regulator_matrix,
            tf_matrix_gene_names=regulator_names,
            target_gene_name=target_name,
            target_gene_expression=target_expr,
            include_meta=False,
            early_stop_window_length=EARLY_STOP_WINDOW_LENGTH,
            seed=seed,
        )
        return result
    except Exception as exc:
        LOGGER.error(f"Error on target {target_name}: {exc}")
        return pd.DataFrame()


def process_batch(batch_id, batch_tasks, args):
    LOGGER.info(f"Processing batch {batch_id + 1}/{args.batch_size} ({len(batch_tasks)} targets)")
    with Pool(args.n_workers) as pool:
        results = list(
            tqdm.tqdm(
                pool.imap(run_network_for_target, batch_tasks),
                total=len(batch_tasks),
                desc=f"Batch {batch_id + 1}",
            )
        )
    valid = [x for x in results if x is not None and not x.empty]
    return pd.concat(valid, ignore_index=True) if valid else pd.DataFrame()


# ──────────────────────────────────────────────────────────────
# Output column standardisation
# ──────────────────────────────────────────────────────────────

def standardize_result_columns(df: pd.DataFrame) -> pd.DataFrame:
    rename_map = {}
    cols = set(df.columns)
    if "regulator" not in cols:
        if "TF" in cols:
            rename_map["TF"] = "regulator"
        elif "feature" in cols:
            rename_map["feature"] = "regulator"
    if "importance" not in cols and "weight" in cols:
        rename_map["weight"] = "importance"
    if rename_map:
        df = df.rename(columns=rename_map)
    missing = {"regulator", "target", "importance"} - set(df.columns)
    if missing:
        raise ValueError(f"Missing expected columns from network result: {sorted(missing)}")
    return df


# ──────────────────────────────────────────────────────────────
# Spearman correlation (vectorised)
# ──────────────────────────────────────────────────────────────

def compute_pairwise_spearman_vectorized(
    regulator_expr_df: pd.DataFrame,
    target_expr_df: pd.DataFrame,
    pairs_df: pd.DataFrame,
) -> np.ndarray:
    """Vectorised Spearman correlation for all (regulator, target) pairs."""
    LOGGER.info("Computing regulator-target Spearman correlations ...")
    reg_rank = regulator_expr_df.rank(axis=0, method="average")
    tar_rank = target_expr_df.rank(axis=0, method="average")
    reg_rank = (reg_rank - reg_rank.mean(axis=0)) / reg_rank.std(axis=0, ddof=0).replace(0, np.nan)
    tar_rank = (tar_rank - tar_rank.mean(axis=0)) / tar_rank.std(axis=0, ddof=0).replace(0, np.nan)
    corr_mat = (reg_rank.T @ tar_rank) / reg_rank.shape[0]
    corr_df = pd.DataFrame(corr_mat.values, index=reg_rank.columns, columns=tar_rank.columns)

    r_idx = corr_df.index.get_indexer(pairs_df["regulator"])
    c_idx = corr_df.columns.get_indexer(pairs_df["target"])
    vals = corr_df.values
    corr = np.full(len(pairs_df), np.nan, dtype=np.float32)
    mask = (r_idx >= 0) & (c_idx >= 0)
    corr[mask] = vals[r_idx[mask], c_idx[mask]]
    return corr


# ──────────────────────────────────────────────────────────────
# CLI registration
# ──────────────────────────────────────────────────────────────

def register_subcommand(subparsers):
    """Register getGRN with the main scRBP CLI."""
    parser = subparsers.add_parser(
        "getGRN",
        help="Gene/Isoform regulatory network inference (GRNBoost2/GENIE3)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )

    # ── Required ──────────────────────────────────────────────
    parser.add_argument(
        "--matrix", required=True,
        help="Expression matrix file (.csv/.csv.gz/.feather/.loom). "
             "Rows=features (genes or isoforms), columns=cells."
    )
    parser.add_argument(
        "--rbp_list", required=True,
        help="Plain-text file of RBP gene symbols (one per line)."
    )
    parser.add_argument(
        "--output", required=True,
        help="Output prefix or TSV path. "
             "The actual filename gains a mode suffix automatically: "
             "<prefix>_scRBP_gene_GRNs.tsv  or  <prefix>_scRBP_isoform_GRNs.tsv."
    )

    # ── Mode ──────────────────────────────────────────────────
    parser.add_argument(
        "--mode", choices=["gene", "isoform"], default="gene",
        help="Inference mode (default: gene).\n"
             "  gene    : RBP → Gene network; RBPs are matched directly from the gene matrix.\n"
             "  isoform : RBP → Isoform network; multiple isoforms of the same RBP are\n"
             "            aggregated into one regulator signal (see --rbp_agg_method).\n"
             "            Requires --isoform_annotation."
    )

    # ── Isoform-mode options ───────────────────────────────────
    iso = parser.add_argument_group("isoform mode options")
    iso.add_argument(
        "--isoform_annotation", default=None,
        help="[isoform mode] TSV/CSV file mapping isoform IDs → gene symbols. "
             "Required columns (auto-detected): isoform/transcript_id  and  gene/gene_symbol."
    )
    iso.add_argument(
        "--rbp_agg_method", choices=["sum", "mean", "max"], default="sum",
        help="[isoform mode] How to aggregate multiple isoforms of the same RBP gene "
             "into one regulator signal (default: sum)."
    )
    iso.add_argument(
        "--remove_self_targets", type=str2bool, default=True,
        help="[isoform mode] Remove edges where the target isoform belongs to the same "
             "gene as the regulator RBP (default: True)."
    )
    iso.add_argument(
        "--save_rbp_agg_matrix", type=str2bool, default=True,
        help="[isoform mode] Save aggregated RBP expression matrix to TSV (default: True)."
    )
    iso.add_argument(
        "--save_rbp_isoform_map", type=str2bool, default=True,
        help="[isoform mode] Save RBP→isoform mapping table to TSV (default: True)."
    )
    iso.add_argument(
        "--min_target_cells_expressed", type=int, default=10,
        help="[isoform mode] Minimum number of cells with expression > 0 to keep a target "
             "isoform (default: 10)."
    )
    iso.add_argument(
        "--min_target_mean_expr", type=float, default=0.01,
        help="[isoform mode] Minimum mean expression across all cells to keep a target "
             "isoform (default: 0.01)."
    )

    # ── Algorithm & performance ────────────────────────────────
    parser.add_argument(
        "--method", choices=["genie3", "grnboost2"], default="grnboost2",
        help="Tree-based algorithm to use (default: grnboost2)."
    )
    parser.add_argument(
        "--n_workers", type=int, default=cpu_count(),
        help="Number of parallel worker processes inside each batch (default: all CPUs)."
    )
    parser.add_argument(
        "--batch_size", type=int, default=10,
        help="Number of outer batches for target genes (default: 10)."
    )

    # ── Correlation ────────────────────────────────────────────
    parser.add_argument(
        "--correlation", type=str2bool, default=True,
        help="Compute Spearman correlation and assign regulatory Mode "
             "(activating/repressing) for each edge (default: True). "
             "Set False to skip (faster)."
    )
    parser.add_argument(
        "--threshold", type=float, default=0.03,
        help="Absolute Spearman correlation threshold. Edges with |r| <= threshold are "
             "removed (default: 0.03). Only used when --correlation True."
    )

    # ── Misc ──────────────────────────────────────────────────
    parser.add_argument(
        "--seed", type=int, default=1234,
        help="Random seed for the tree ensemble (default: 1234)."
    )
    parser.add_argument(
        "--log", type=str, default=None,
        help="Path to a log file. If not set, logs go to stderr."
    )
    parser.add_argument(
        "--debug", action="store_true",
        help="Enable verbose DEBUG logging."
    )

    parser.set_defaults(func=main)
    return parser


# ──────────────────────────────────────────────────────────────
# Main execution
# ──────────────────────────────────────────────────────────────

def main(args):
    start_time = time.time()

    LOGGER.handlers.clear()
    LOGGER.addHandler(create_logging_handler(args.debug, args.log))
    LOGGER.setLevel(logging.DEBUG if args.debug else logging.INFO)

    LOGGER.info("Starting scRBP network inference pipeline")
    LOGGER.info(f"Mode   : {args.mode}")
    LOGGER.info(f"Method : {args.method}")

    # Step 1: Load expression matrix
    expr_cells_x_features = load_expression_matrix(args.matrix)

    # Step 2: Load RBP list
    rbps = list(pd.unique(pd.Series(load_tf_names(args.rbp_list)).astype(str)))
    LOGGER.info(f"Input RBP list size: {len(rbps)}")

    # Step 3: Build regulator matrix + target list (mode-dependent)
    isoform_to_gene_map = None

    if args.mode == "gene":
        rbps_in_expr = sorted(set(expr_cells_x_features.columns) & set(rbps))
        if not rbps_in_expr:
            raise ValueError("No RBP found in expression matrix for gene mode.")
        regulator_expr_df = expr_cells_x_features[rbps_in_expr].copy()
        regulator_names   = rbps_in_expr
        target_names      = list(expr_cells_x_features.columns)
        target_expr_df    = expr_cells_x_features.copy()
        LOGGER.info(f"Matched RBPs: {len(regulator_names)}  |  Target genes: {len(target_names)}")

    elif args.mode == "isoform":
        if args.isoform_annotation is None:
            raise ValueError("--isoform_annotation is required when --mode isoform")

        isoform_anno = load_isoform_annotation(args.isoform_annotation)
        isoform_to_gene_map = dict(zip(isoform_anno["Isoform"], isoform_anno["Gene"]))

        (regulator_expr_df, rbp_isoform_map_df,
         annotation_summary_df, gene_to_isoforms) = aggregate_rbp_expression_from_isoforms(
            expr_cells_x_isoforms=expr_cells_x_features,
            rbp_genes=rbps,
            isoform_annotation=isoform_anno,
            agg_method=args.rbp_agg_method,
        )
        regulator_names = list(regulator_expr_df.columns)

        # Save isoform-mode auxiliary files
        if args.save_rbp_agg_matrix:
            agg_out = build_aux_output_path(args.output, "scRBP_isoform_RBP_aggregated_expression.tsv")
            regulator_expr_df.to_csv(agg_out, sep="\t", index=True)
            LOGGER.info(f"Saved aggregated RBP matrix → {agg_out}")

        if args.save_rbp_isoform_map:
            map_out = build_aux_output_path(args.output, "scRBP_isoform_RBP_isoform_map.tsv")
            rbp_isoform_map_df.to_csv(map_out, sep="\t", index=False)
            LOGGER.info(f"Saved RBP-isoform map → {map_out}")

        summary_out = build_aux_output_path(args.output, "scRBP_isoform_annotation_summary.tsv")
        annotation_summary_df.to_csv(summary_out, sep="\t", index=False)
        LOGGER.info(f"Saved annotation summary → {summary_out}")

        # Filter low-expressed target isoforms
        target_names, target_filter_df = filter_target_isoforms(
            expr_cells_x_isoforms=expr_cells_x_features,
            min_cells_expressed=args.min_target_cells_expressed,
            min_mean_expr=args.min_target_mean_expr,
        )
        target_filter_out = build_aux_output_path(args.output, "scRBP_isoform_target_filter_summary.tsv")
        target_filter_df.to_csv(target_filter_out, sep="\t", index=False)
        LOGGER.info(f"Saved target filter summary → {target_filter_out}")

        if not target_names:
            raise ValueError("No isoform targets left after filtering.")

        target_expr_df = expr_cells_x_features[target_names].copy()
        LOGGER.info(f"Aggregated RBPs: {len(regulator_names)}  |  Filtered target isoforms: {len(target_names)}")

    else:
        raise ValueError(f"Unsupported mode: {args.mode}")

    regulator_matrix = regulator_expr_df.values

    # Step 4: Build tasks
    tasks = [
        (tgt, target_expr_df[tgt].values, regulator_matrix, regulator_names, args.method, args.seed)
        for tgt in target_expr_df.columns
    ]
    if not tasks:
        raise ValueError("No targets available for network inference.")

    batch_chunk = max(1, len(tasks) // args.batch_size)
    batches     = [tasks[i : i + batch_chunk] for i in range(0, len(tasks), batch_chunk)]
    LOGGER.info(f"Total targets: {len(tasks)}  |  Batches: {len(batches)}")

    # Step 5: Run inference in parallel
    final_result = pd.DataFrame()
    with ProcessPoolExecutor(max_workers=args.batch_size) as executor:
        futures = {executor.submit(process_batch, i, b, args): i for i, b in enumerate(batches)}
        for future in as_completed(futures):
            bid = futures[future]
            try:
                res = future.result()
                if res is not None and not res.empty:
                    final_result = pd.concat([final_result, res], ignore_index=True)
                else:
                    LOGGER.warning(f"Batch {bid + 1} returned no results.")
            except Exception as exc:
                LOGGER.error(f"Batch {bid + 1} failed: {exc}")

    if final_result.empty:
        LOGGER.error("No network inference results produced.")
        sys.exit(1)

    final_result = standardize_result_columns(final_result)

    # Step 6: Remove self-targets (isoform mode only)
    if args.mode == "isoform" and args.remove_self_targets:
        regulator_gene = final_result["regulator"].astype(str)
        target_gene    = final_result["target"].map(isoform_to_gene_map).astype(str)
        self_mask      = regulator_gene == target_gene
        n_self         = int(self_mask.sum())

        if n_self > 0:
            self_out = build_aux_output_path(args.output, "scRBP_isoform_self_targets_removed.tsv")
            final_result.loc[self_mask].to_csv(self_out, sep="\t", index=False)
            LOGGER.info(f"Saved removed self-target edges → {self_out}")

        final_result = final_result.loc[~self_mask].copy()
        LOGGER.info(f"Self-target removal: {n_self} edges removed ({len(final_result)} remaining)")

        pd.DataFrame([
            {"metric": "remove_self_targets",               "value": str(args.remove_self_targets)},
            {"metric": "n_self_targets_removed",            "value": int(n_self)},
            {"metric": "n_edges_after_self_target_filter",  "value": int(len(final_result))},
        ]).to_csv(
            build_aux_output_path(args.output, "scRBP_isoform_self_target_summary.tsv"),
            sep="\t", index=False,
        )

    # Step 7: Spearman correlation & Mode
    if args.correlation:
        LOGGER.info("Computing Spearman correlation and assigning Mode ...")
        reg_log2 = np.log2(regulator_expr_df + 1.0)
        unique_tgts = pd.unique(final_result["target"])
        tgt_log2 = np.log2(target_expr_df[unique_tgts] + 1.0)

        corr = compute_pairwise_spearman_vectorized(
            regulator_expr_df=reg_log2,
            target_expr_df=tgt_log2,
            pairs_df=final_result[["regulator", "target"]],
        )
        final_result["Correlation"] = corr

        before_n = len(final_result)
        final_result = final_result[final_result["Correlation"].abs() > args.threshold].copy()
        LOGGER.info(f"Correlation filter (|r|>{args.threshold}): {before_n} → {len(final_result)}")

        final_result["Mode"] = np.where(
            final_result["Correlation"] > 0, "activating", "repressing"
        )
    else:
        LOGGER.info("Skipping correlation & mode assignment (--correlation=False).")
        final_result["Correlation"] = np.nan
        final_result["Mode"] = pd.NA

    # Step 8: Rename columns & save
    if args.mode == "gene":
        final_result = final_result.rename(
            columns={"regulator": "RBP", "target": "Gene", "importance": "Importance"}
        )[["RBP", "Gene", "Importance", "Correlation", "Mode"]]
    else:
        final_result = final_result.rename(
            columns={"regulator": "RBP", "target": "Isoform", "importance": "Importance"}
        )[["RBP", "Isoform", "Importance", "Correlation", "Mode"]]

    final_result = final_result.sort_values("Importance", ascending=False)
    output_path  = build_main_output_path(args.output, args.mode)
    final_result.to_csv(output_path, sep="\t", index=False)

    LOGGER.info(f"Final result saved → {output_path}")
    LOGGER.info(f"Pipeline completed in {time.time() - start_time:.2f} s")
