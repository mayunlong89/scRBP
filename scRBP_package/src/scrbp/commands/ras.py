#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scRBP ras — Regulon Activity Score (RAS) computation via AUCell

Quantifies RBP regulon activity at single-cell or cell-type resolution by
evaluating whether the target genes of each regulon are preferentially
enriched among the most highly expressed genes in each cell.  The RAS is
computed as the normalised area under the cumulative recovery curve within
the top-ranked gene window (default: top 5%), using the AUCell algorithm
from pyscenic.  Because AUCell uses within-cell expression ranks rather
than absolute values, the RAS is robust to differences in expression scale
and normalisation procedures.

In cell-type mode (``--mode ct``), single-cell RAS values are aggregated
per cell type using the Jensen–Shannon divergence-based Regulon Specificity
Score (RSS), which measures how specifically each regulon is active in a
given cell type relative to all others.

Outputs
-------
- AUCell score matrix (regulons × cells): CSV and/or LOOM format
- RSS matrix (cell types × regulons): ``<base>_rss.csv``
- Per-gene expression statistics for MAGMA null matching: ``<base>_expr_stats.tsv``

Supported expression input formats: .csv, .csv.gz, .feather, .loom
  (rows = genes, columns = cells; transposed automatically)
Supported regulon input formats: .gmt, .gmt.gz, .pkl/.pickle
"""

from __future__ import annotations

# ---- prevent BLAS over-parallelism (before numpy/pandas import) ----
import os as _os
for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_MAX_THREADS"):
    _os.environ.setdefault(_k, "1")

import io
import gzip
import time
import argparse
import threading
import multiprocessing
import pickle
from typing import Dict, Iterable, List, Set, Optional, Tuple

import numpy as np
import pandas as pd
import polars as pl

from pyscenic.aucell import aucell
from ctxcore.genesig import Regulon

from scrbp.utils.io import (
    check_matrix_path,
    load_expression_matrix,
    compute_expr_stats_from_df,
    read_gmt,
)
from scrbp.utils.stats import regulon_specificity_scores


# =============================
# Progress (time bar)
# =============================
class TimeBar:
    def __init__(self, desc: str = "Working"):
        self.desc = desc
        self._stop = threading.Event()
        self._thread = None
        self._use_tqdm = False
        self._t0 = None

    def __enter__(self):
        self._t0 = time.time()
        try:
            from tqdm import tqdm  # noqa: F401
            self._use_tqdm = True
        except Exception:
            self._use_tqdm = False

        if self._use_tqdm:
            from tqdm import tqdm
            def runner():
                with tqdm(total=0,
                          bar_format="{desc} | elapsed: {elapsed}",
                          desc=self.desc,
                          dynamic_ncols=True,
                          leave=True) as pbar:
                    while not self._stop.is_set():
                        time.sleep(0.2)
                        pbar.update(0)
            self._thread = threading.Thread(target=runner, daemon=True)
            self._thread.start()
        else:
            print(f"[{self.desc}] started ...")
        return self

    def __exit__(self, exc_type, exc, tb):
        self._stop.set()
        if self._thread is not None:
            self._thread.join()
        dt = time.time() - self._t0 if self._t0 is not None else 0.0
        print(f"[{self.desc}] completed in {dt:.2f}s")


# =============================
# Path checks
# =============================
def _check_exists(path: str, msg: str):
    if not path or not _os.path.exists(path):
        raise FileNotFoundError(msg + f": {path}")


# check_matrix_path moved to scrbp.utils.io


def check_regulon_path(path: str):
    _check_exists(path, "Regulons file not found")
    if not path.lower().endswith((".pkl", ".pickle", ".gmt", ".gmt.gz")):
        raise ValueError(f"Regulons must be .pkl/.pickle or .gmt/.gmt.gz, got: {path}")


def check_aucell_in_path(path: str):
    _check_exists(path, "AUCell file not found")
    if not path.lower().endswith((".csv", ".csv.gz", ".loom")):
        raise ValueError(f"AUCell must be .csv/.csv.gz or .loom, got: {path}")


# load_expression_matrix moved to scrbp.utils.io


# =============================
# Regulon loader (PKL / GMT)
# =============================
def _iterable_genes(x: Iterable) -> List[str]:
    return [str(g).strip() for g in x if str(g).strip()]


def _ensure_regulon_dict(obj) -> Dict[str, Set[str]]:
    if isinstance(obj, dict):
        return {str(k): set(_iterable_genes(v)) for k, v in obj.items()}
    if isinstance(obj, list) and all(isinstance(x, (tuple, list)) and len(x) >= 2 for x in obj):
        return {str(name): set(_iterable_genes(genes)) for (name, genes) in obj}
    try:
        from ctxcore.genesig import Regulon as _R
        if isinstance(obj, list) and all(isinstance(x, _R) for x in obj):
            return {x.name: set(_iterable_genes(x.genes)) for x in obj}
    except Exception:
        pass
    raise ValueError("Unsupported PKL content. Expect dict{name->genes}, list[(name, genes)], or list[Regulon].")


# _load_from_gmt replaced by scrbp.utils.io.read_gmt(path, as_sets=True)


def load_regulons(file_path: str, *,
                  to_upper: bool = False,
                  allowed_genes: Optional[Set[str]] = None,
                  min_genes: int = 1) -> List[Regulon]:
    print(f"[Regulons] Loading from {file_path} ...")
    ext = file_path.lower()

    if ext.endswith((".pkl", ".pickle")):
        with open(file_path, "rb") as f:
            obj = pickle.load(f)
        reg_dict = _ensure_regulon_dict(obj)
    elif ext.endswith((".gmt", ".gmt.gz")):
        reg_dict = read_gmt(file_path, as_sets=True)
    else:
        raise ValueError(f"Unsupported file type: {file_path}")

    if to_upper:
        reg_dict = {k: {g.upper() for g in v} for k, v in reg_dict.items()}
    if allowed_genes is not None:
        reg_dict = {k: (set(v) & allowed_genes) for k, v in reg_dict.items()}
    if min_genes > 1:
        reg_dict = {k: v for k, v in reg_dict.items() if len(v) >= min_genes}

    print(f"[Regulons] Kept: {len(reg_dict)} sets")

    out: List[Regulon] = []
    for name, genes in reg_dict.items():
        genes = sorted(set(genes))
        out.append(Regulon(
            name=name,
            gene2weight={g: 1.0 for g in genes},
            gene2occurrence={g: 1 for g in genes},
            transcription_factor=name,
        ))
    return out


# =============================
# AUCell input helpers (cells x regulons)
# =============================
def load_aucell_csv_sets_by_cells(path: str) -> pd.DataFrame:
    # 第一列=Regulon，其余=Cell  -> (regs x cells)，转成 cells x regs
    df = pd.read_csv(path, low_memory=False)
    if df.shape[1] < 2:
        raise ValueError(f"[AUCell CSV] Expect first col=Regulon, others=Cells. Got shape={df.shape}")
    set_col = df.columns[0]
    m = (df.set_index(set_col).apply(pd.to_numeric, errors="coerce").fillna(0.0))
    auc = m.T
    auc = auc.loc[~auc.index.duplicated(keep="first")]
    auc = auc.loc[:, ~pd.Index(auc.columns).duplicated(keep="first")]
    print(f"[AUCell CSV sets-by-cells] Loaded: {auc.shape} (cells x regulons)")
    return auc


def load_aucell_csv_cells_by_sets(path: str) -> pd.DataFrame:
    # 第一列=Cell，其余=Regulon -> (cells x regs)
    df = pd.read_csv(path, low_memory=False)
    if df.shape[1] < 2:
        raise ValueError(f"[AUCell CSV] Expect first col=Cell, others=Regulons. Got shape={df.shape}")
    cell_col = df.columns[0]
    auc = (df.set_index(cell_col).apply(pd.to_numeric, errors="coerce").fillna(0.0))
    auc = auc.loc[~auc.index.duplicated(keep="first")]
    auc = auc.loc[:, ~pd.Index(auc.columns).duplicated(keep="first")]
    print(f"[AUCell CSV cells-by-sets] Loaded: {auc.shape} (cells x regulons)")
    return auc


def load_aucell_loom(path: str, dtype="float32") -> pd.DataFrame:
    # 约定 LOOM: rows=regulons, cols=cells -> 转 cells x regs
    try:
        import loompy
    except ImportError as e:
        raise ImportError("Please install loompy: pip install loompy") from e
    with loompy.connect(path, mode="r") as ds:
        regs  = np.array(ds.ra["Regulon"], dtype=object) if "Regulon" in ds.ra else np.arange(ds.shape[0], dtype=object)
        cells = np.array(ds.ca["CellID"], dtype=object)  if "CellID"  in ds.ca else np.arange(ds.shape[1], dtype=object)
        mat = ds[:, :]
    auc = pd.DataFrame(mat.T, index=cells, columns=regs)
    auc = auc.loc[~auc.index.duplicated(keep="first")]
    auc = auc.loc[:, ~pd.Index(auc.columns).duplicated(keep="first")]
    print(f"[AUCell LOOM] Loaded: {auc.shape} (cells x regulons)")
    return auc


def load_aucell_auto(path: str, labels: Optional[pd.Series] = None) -> pd.DataFrame:
    lower = path.lower()
    if lower.endswith(".loom"):
        return load_aucell_loom(path)

    try_a = load_aucell_csv_sets_by_cells(path)
    try:
        try_b = load_aucell_csv_cells_by_sets(path)
    except Exception:
        try_b = None

    if labels is None or try_b is None:
        return try_a if try_b is None else (try_b if try_b.shape[0] >= try_a.shape[0] else try_a)

    a_match = pd.Index(try_a.index).intersection(labels.index).size
    b_match = pd.Index(try_b.index).intersection(labels.index).size
    if b_match >= a_match:
        print(f"[AUCell CSV auto] Chose cells-by-sets (matches {b_match} cells).")
        return try_b
    else:
        print(f"[AUCell CSV auto] Chose sets-by-cells (matches {a_match} cells).")
        return try_a


# =============================
# Orientation helper (computed AUCell)
# =============================
def orient_aucell_to_cells_by_regulons(auc_df: pd.DataFrame,
                                       cell_ids: Iterable[str],
                                       regulon_names: Iterable[str]) -> pd.DataFrame:
    cells = set(map(str, cell_ids))
    regs  = set(map(str, regulon_names))
    r_rows = len(set(map(str, auc_df.index))   & regs)
    c_rows = len(set(map(str, auc_df.index))   & cells)
    r_cols = len(set(map(str, auc_df.columns)) & regs)
    c_cols = len(set(map(str, auc_df.columns)) & cells)
    if c_rows >= r_rows and r_cols >= c_cols:  # already cells x regs
        return auc_df
    if r_rows >= c_rows and c_cols >= r_cols:  # regs x cells
        return auc_df.T
    n_cells, n_regs = len(cells), len(regs)
    if auc_df.shape == (n_cells, n_regs): return auc_df
    if auc_df.shape == (n_regs, n_cells): return auc_df.T
    raise ValueError(f"Cannot infer AUCell orientation: shape={auc_df.shape}")


# regulon_specificity_scores_scrbp moved to scrbp.utils.stats as regulon_specificity_scores


# =============================
# Cell-type CSV loader
# =============================
def load_celltypes_table(path: str,
                         cell_col: Optional[str],
                         ctype_col: Optional[str]) -> pd.DataFrame:
    if path is None or not _os.path.exists(path):
        raise FileNotFoundError("--celltypes-csv is required for --mode ct and must exist.")
    df = pd.read_csv(path)
    if cell_col is None:
        cand = [c for c in df.columns if c.lower() in {"cell", "cell_id", "barcode", "barcodes", "obs_names", "cellid"}]
        if not cand: raise ValueError("Cannot auto-detect cell ID column; specify --cell-col.")
        cell_col = cand[0]
    if ctype_col is None:
        cand = [c for c in df.columns if c.lower() in {"cell_type", "celltype", "celltype_major", "cell_type_major", "cluster", "annotation", "ctype", "celltypes"}]
        if not cand: raise ValueError("Cannot auto-detect cell type column; specify --ctype-col.")
        ctype_col = cand[0]
    sub = df[[cell_col, ctype_col]].copy()
    sub.columns = ["cell", "cell_type"]
    sub["cell"] = sub["cell"].astype(str)
    sub["cell_type"] = sub["cell_type"].astype(str)
    return sub


# =============================
# Output helpers
# =============================
def _strip_ext(path: str) -> str:
    if path.lower().endswith(".csv.gz"): return path[:-7]
    if path.lower().endswith(".csv"):    return path[:-4]
    if path.lower().endswith(".loom"):   return path[:-5]
    return path


def resolve_paths(out_arg: str, out_format: str, csv_layout: str):
    """
    Returns:
      csv_paths: List[(path, layout)]  layout in {'regulons_by_cells','cells_by_regulons'}
      loom_path: Optional[str]
      base_for_rss: str  (used for *_rss.csv / *_expr_stats.tsv)
    """
    base = _strip_ext(out_arg)
    csv_paths: List[Tuple[str, str]] = []
    loom_path = None

    if out_format in ("csv", "both"):
        if csv_layout == "both":
            csv_paths.append((f"{base}_regulons_by_cells.csv", "regulons_by_cells"))
            csv_paths.append((f"{base}_cells_by_regulons.csv", "cells_by_regulons"))
        elif out_arg.lower().endswith((".csv", ".csv.gz")):
            csv_paths.append((out_arg, csv_layout))
        else:
            csv_paths.append((f"{base}.csv", csv_layout))

    if out_format in ("loom", "both"):
        loom_path = f"{base}.loom" if not out_arg.lower().endswith(".loom") else out_arg

    return csv_paths, loom_path, base


def save_csv_oriented(auc_cells_by_regs: pd.DataFrame, path: str, layout: str):
    _os.makedirs(_os.path.dirname(path) or ".", exist_ok=True)
    if layout == "regulons_by_cells":
        df = auc_cells_by_regs.T
        df.index.name = "Regulon"
    else:  # cells_by_regulons
        df = auc_cells_by_regs.copy()
        df.index.name = "cell"
    df.to_csv(path)


def save_loom_from_sc(auc_df_cell_by_reg: pd.DataFrame, path: str, dtype="float32"):
    try:
        import loompy
    except ImportError as e:
        raise ImportError("Please install loompy: pip install loompy") from e
    _os.makedirs(_os.path.dirname(path) or ".", exist_ok=True)
    mat = auc_df_cell_by_reg.to_numpy(dtype=dtype).T  # rows=regulons, cols=cells
    row_attrs = {"Regulon": np.array(auc_df_cell_by_reg.columns.values, dtype=object)}
    col_attrs = {"CellID": np.array(auc_df_cell_by_reg.index.values, dtype=object)}
    loompy.create(path, mat, row_attrs, col_attrs)


# compute_expr_stats_from_df moved to scrbp.utils.io


def save_expr_stats(expr_df_cells_by_genes: pd.DataFrame, out_tsv: str):
    _os.makedirs(_os.path.dirname(out_tsv) or ".", exist_ok=True)
    df = compute_expr_stats_from_df(expr_df_cells_by_genes)
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[ExprStats] Wrote: {out_tsv} | genes={df.shape[0]} | units={expr_df_cells_by_genes.shape[0]}")


def register_subcommand(subparsers):
    """Register this subcommand with the main parser."""
    p = subparsers.add_parser(
        'ras',
        help='Compute Regulon Activity Scores (AUCell + RSS)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )
    # Modes
    p.add_argument("--mode", choices=["sc", "ct"], default="ct",
                   help="sc: only AUCell; ct: AUCell -> RSS (default)")

    # Compute AUCell from expression
    p.add_argument("--matrix",
                   help="Expression (.csv/.csv.gz/.feather/.loom). Rows=genes, Cols=cells (auto-transpose).")
    p.add_argument("--regulons",
                   help="Regulons (.pkl/.pickle or .gmt/.gmt.gz).")

    # Or: read precomputed AUCell
    p.add_argument("--aucell-in",
                   help="Precomputed AUCell (.csv/.csv.gz/.loom). CSV layout auto-detected.")

    # Outputs
    p.add_argument("--out", required=True,
                   help="SC output path/prefix. If endswith .csv/.csv.gz or .loom it's respected; otherwise used as prefix.")
    p.add_argument("--out_format", choices=["csv", "loom", "both"], default="csv",
                   help="Which SC outputs to write (RSS always CSV).")
    p.add_argument("--csv_layout",
                   choices=["regulons_by_cells", "cells_by_regulons", "both"],
                   default="regulons_by_cells",
                   help="CSV orientation. Default 'regulons_by_cells' matches scRBP_trs.")

    # Performance / preprocessing
    p.add_argument("--n_workers", type=int, default=4, help="Workers for AUCell.")
    p.add_argument("--to_upper", action="store_true", help="Uppercase gene symbols in regulons.")
    p.add_argument("--min_genes", type=int, default=1, help="Drop regulons with < min_genes after filtering.")
    p.add_argument("--dtype", default="float32", help="Numeric dtype for expression (default: float32).")
    p.add_argument("--log", default="run_time.log", help="Runtime log file.")

    # CT mode
    p.add_argument("--celltypes-csv", help="CSV with at least (cell, cell_type) columns. Required for --mode ct.")
    p.add_argument("--cell-col", default=None, help="Column name for cell IDs in --celltypes-csv (auto-detect if omitted).")
    p.add_argument("--ctype-col", default=None, help="Column name for cell types in --celltypes-csv (auto-detect if omitted).")
    p.add_argument("--rss-eps", type=float, default=1e-12, help="Epsilon for JSD logs.")

    # Expr-stats controls (always available)
    p.add_argument("--emit-expr-stats", action="store_true", default=True,
                   help="Emit <out>_expr_stats.tsv from the FULL input matrix (default: True).")
    p.add_argument("--expr-stats-out", default=None,
                   help="Custom path for expr-stats TSV. Default: <out>_expr_stats.tsv")
    p.set_defaults(func=main)
    return p


def main(args):
    """Main execution function."""
    t0 = time.time()

    # ---- read full matrix (for AUCell and expr-stats) ----
    expr_df_full: Optional[pd.DataFrame] = None
    if args.matrix:
        check_matrix_path(args.matrix)
        with TimeBar("Load FULL expression matrix"):
            expr_df_full = load_expression_matrix(args.matrix, dtype=args.dtype)

    # ---- AUCell: compute or load (internal: cells x regulons) ----
    if args.aucell_in:
        check_aucell_in_path(args.aucell_in)
        labels_for_auto = None
        if args.mode == "ct" and args.celltypes_csv:
            ct_df_tmp = load_celltypes_table(args.celltypes_csv, args.cell_col, args.ctype_col)
            labels_for_auto = pd.Series(ct_df_tmp["cell_type"].values, index=ct_df_tmp["cell"].values, name="cell_type")
        with TimeBar("Load AUCell (precomputed)"):
            auc_df = load_aucell_auto(args.aucell_in, labels_for_auto)
    else:
        if expr_df_full is None or not args.regulons:
            raise ValueError("Provide --matrix and --regulons to compute AUCell (or use --aucell-in).")
        check_regulon_path(args.regulons)

        with TimeBar("Load regulons"):
            regulons_list = load_regulons(args.regulons, to_upper=args.to_upper,
                                          allowed_genes=set(expr_df_full.columns), min_genes=args.min_genes)
        n_sets = len(regulons_list)
        if n_sets == 0:
            raise ValueError("No regulons kept after filtering/intersection.")
        cpu_cap = multiprocessing.cpu_count()
        num_workers = max(1, min(cpu_cap, args.n_workers, n_sets))
        print(f"[scRBP_ras] Regulons: {n_sets} | CPU: {cpu_cap} | requested: {args.n_workers} -> using: {num_workers}")

        print("[scRBP_ras] Computing AUCell ...")
        with TimeBar("AUCell"):
            auc_raw = aucell(expr_df_full, regulons_list, num_workers=num_workers)

        reg_names = [r.name for r in regulons_list]
        auc_df = orient_aucell_to_cells_by_regulons(auc_raw, expr_df_full.index, reg_names)
        print(f"[AUCell] Oriented to cells x regulons: {auc_df.shape}")

    # ---- write SC output ----
    csv_paths, loom_path, base = resolve_paths(args.out, args.out_format, args.csv_layout)
    for pth, lay in csv_paths:
        with TimeBar(f"Write SC CSV ({lay}): {pth}"):
            save_csv_oriented(auc_df, pth, lay)
    if loom_path:
        with TimeBar(f"Write SC LOOM: {loom_path}"):
            save_loom_from_sc(auc_df, loom_path, dtype=args.dtype)

    # ---- Expr-stats (full matrix) ----
    if args.emit_expr_stats:
        if expr_df_full is None:
            raise ValueError("--emit-expr-stats needs --matrix to compute full-matrix stats.")
        expr_out = args.expr_stats_out if args.expr_stats_out else f"{base}_expr_stats.tsv"
        with TimeBar(f"Write expr-stats: {expr_out}"):
            save_expr_stats(expr_df_full, expr_out)

    # ---- CT: compute RSS ----
    if args.mode == "ct":
        if not args.celltypes_csv:
            raise ValueError("--mode ct requires --celltypes-csv")
        with TimeBar("Load cell types"):
            ct_df = load_celltypes_table(args.celltypes_csv, args.cell_col, args.ctype_col)
            ct_ser = pd.Series(ct_df["cell_type"].values, index=ct_df["cell"].values, name="cell_type")
            ct_ser = ct_ser.reindex(auc_df.index)

        matched = ct_ser.notna().sum()
        total   = len(ct_ser)
        print(f"[RSS] cell-type matches: {matched}/{total} ({(matched/total if total else 0):.1%}). "
              f"Unlabeled cells will be dropped.")

        with TimeBar("Compute RSS"):
            rss_df = regulon_specificity_scores(auc_df, ct_ser, eps=args.rss_eps)
        rss_out = f"{base}_rss.csv"
        print(f"[RSS] Shape: {rss_df.shape} (cell_types x regulons)")
        with TimeBar(f"Write RSS CSV: {rss_out}"):
            rss_df.to_csv(rss_out)

    # ---- Logging ----
    elapsed = time.time() - t0
    msg = f"scRBP_ras v4 (no-prefilter) completed in {elapsed:.2f}s | CSV files: {len(csv_paths)} | LOOM: {bool(loom_path)}"
    if args.mode == "ct": msg += " | RSS: True"
    if args.emit_expr_stats: msg += " | ExprStats: True"
    msg += "\n"
    print(msg)
    with open(args.log, "a") as f:
        f.write(msg)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    subp = p.add_subparsers()
    register_subcommand(subp)
    args = p.parse_args()
    main(args)
