#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRBP shared I/O utilities
==========================
Functions shared across multiple scRBP commands:

  - check_matrix_path         (get_grn, ras)
  - load_expression_matrix    (get_grn, ras, rgs)
  - compute_expr_stats_from_df(ras, rgs)
  - read_gmt                  (ras, rgs)
  - write_gmt                 (get_regulon, rgs)
"""

from __future__ import annotations

import gzip
import io
import os
from typing import Dict, List, Optional, Set, Union

import numpy as np
import pandas as pd


# ──────────────────────────────────────────────────────────────
# 1. Path validation
# ──────────────────────────────────────────────────────────────

_MATRIX_EXTS = (".csv", ".csv.gz", ".feather", ".loom")
_GMT_EXTS    = (".gmt", ".gmt.gz")


def check_matrix_path(path: str) -> None:
    """
    Validate that *path* exists and has a supported expression-matrix extension.

    Supported extensions: .csv  .csv.gz  .feather  .loom

    Parameters
    ----------
    path : str
        File path to validate.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the extension is not supported.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Matrix file not found: {path}")
    if not path.lower().endswith(_MATRIX_EXTS):
        raise ValueError(
            f"Matrix must be one of {_MATRIX_EXTS}, got: {path}"
        )


# ──────────────────────────────────────────────────────────────
# 2. Expression matrix loader
# ──────────────────────────────────────────────────────────────

def load_expression_matrix(
    matrix_path: str,
    keep_genes: Optional[Set[str]] = None,
    dtype: str = "float32",
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Load a gene-by-cell expression matrix and return a **cells × genes** DataFrame.

    The function handles all formats used across scRBP steps (getGRN, ras, rgs).
    Input files are assumed to have **rows = genes** and **columns = cells**;
    the function transposes automatically.

    Parameters
    ----------
    matrix_path : str
        Path to the expression matrix file.
        Supported formats: ``.csv``, ``.csv.gz``, ``.feather``, ``.loom``.
    keep_genes : set of str, optional
        If provided, only these gene names are retained after loading.
        Genes not present in the file are silently ignored.
        When ``None`` (default), all genes are returned.
    dtype : str, default ``"float32"``
        Numeric dtype to cast the matrix to after loading.
        Use ``"float32"`` (default) to reduce memory usage.
        Pass ``None`` to keep the file's native dtype.
    verbose : bool, default ``True``
        Print progress messages.

    Returns
    -------
    pd.DataFrame
        Shape: ``(n_cells, n_genes)``.  Index = cell barcodes, columns = gene names.
        Duplicate cell/gene names are deduplicated (first occurrence kept).

    Raises
    ------
    ValueError
        For unsupported file formats.
    ImportError
        If ``loompy`` or ``polars`` is not installed when the respective
        format is requested.

    Notes
    -----
    * Originally implemented independently in ``get_grn.py``, ``ras.py``,
      and ``rgs.py``.  Consolidated here in v0.1.0.
    * The ``keep_genes`` parameter is used by ``getGRN`` to restrict the
      matrix to RBP-list genes before inference.

    Examples
    --------
    >>> df = load_expression_matrix("expr.feather")
    >>> df.shape          # (n_cells, n_genes)
    (15000, 33694)

    >>> df = load_expression_matrix("expr.feather", keep_genes={"HNRNPA1","PTBP1"})
    >>> df.shape
    (15000, 2)
    """
    def _log(msg: str) -> None:
        if verbose:
            print(msg)

    lower = matrix_path.lower()

    # ── CSV / CSV.GZ ──────────────────────────────────────────
    if lower.endswith((".csv", ".csv.gz")):
        try:
            import polars as pl
        except ImportError as exc:
            raise ImportError(
                "polars is required to load CSV matrices: pip install polars"
            ) from exc

        _log(f"[Matrix] Loading CSV via Polars: {matrix_path}")
        if lower.endswith(".gz"):
            with gzip.open(matrix_path, "rt") as fh:
                df_pl = pl.read_csv(fh, has_header=True)
        else:
            df_pl = pl.read_csv(matrix_path, has_header=True)

        index_col = df_pl.columns[0]
        if keep_genes is not None:
            df_pl = df_pl.filter(pl.col(index_col).is_in(list(keep_genes)))
        df = df_pl.to_pandas()
        df.set_index(index_col, inplace=True)
        df = df.T  # → cells × genes

    # ── FEATHER ───────────────────────────────────────────────
    elif lower.endswith(".feather"):
        _log(f"[Matrix] Loading Feather via pandas: {matrix_path}")
        df = pd.read_feather(matrix_path)
        idx_col = df.columns[0]
        df.set_index(idx_col, inplace=True)
        if keep_genes is not None:
            df = df.loc[df.index.intersection(keep_genes)]
        df = df.T  # → cells × genes

    # ── LOOM ─────────────────────────────────────────────────
    elif lower.endswith(".loom"):
        _log(f"[Matrix] Loading Loom via loompy: {matrix_path}")
        try:
            import loompy
        except ImportError as exc:
            raise ImportError(
                "loompy is required to load .loom files: pip install loompy"
            ) from exc

        ds = loompy.connect(matrix_path, mode="r")

        # resolve gene names
        _gene_keys = ("Gene", "gene", "gene_symbol", "gene_names", "var_names", "Genes")
        gene_names = next(
            (ds.row_attrs[k] for k in _gene_keys if k in ds.row_attrs),
            [f"gene_{i}" for i in range(ds.shape[0])],
        )

        # resolve cell IDs
        _cell_keys = ("CellID", "cell", "cell_names", "obs_names", "Cells")
        cell_names = next(
            (ds.col_attrs[k] for k in _cell_keys if k in ds.col_attrs),
            [f"cell_{i}" for i in range(ds.shape[1])],
        )

        if keep_genes is not None:
            g2i = {g: i for i, g in enumerate(gene_names)}
            sel = sorted(g2i[g] for g in keep_genes if g in g2i)
            if not sel:
                ds.close()
                raise ValueError("None of keep_genes found in loom row attributes.")
            matrix = ds[sel, :]
            cols   = [gene_names[i] for i in sel]
        else:
            matrix = ds[:, :]
            cols   = list(gene_names)

        ds.close()
        df = pd.DataFrame(matrix.T, index=cell_names, columns=cols)

    else:
        raise ValueError(
            f"Unsupported matrix format: {matrix_path!r}. "
            f"Accepted: {_MATRIX_EXTS}"
        )

    # ── Post-processing ───────────────────────────────────────
    df = df.loc[~df.index.duplicated(keep="first")]
    df = df.loc[:, ~pd.Index(df.columns).duplicated(keep="first")]
    if keep_genes is not None:
        df = df.loc[:, df.columns.intersection(keep_genes)]
    if dtype:
        df = df.astype(dtype, copy=False)

    _log(f"[Matrix] Loaded: {df.shape} (cells × genes)")
    return df


# ──────────────────────────────────────────────────────────────
# 3. Expression statistics
# ──────────────────────────────────────────────────────────────

def compute_expr_stats_from_df(
    expr_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Compute per-gene expression statistics from a **cells × genes** DataFrame.

    Used by ``ras`` (to export statistics for ``rgs`` 4-D null matching)
    and by ``rgs`` (to compute statistics on-the-fly).

    Parameters
    ----------
    expr_df : pd.DataFrame
        Shape ``(n_cells, n_genes)``.  Values should be raw or normalised counts.

    Returns
    -------
    pd.DataFrame
        One row per gene with columns:

        * ``symbol``       – gene name (str)
        * ``mean_expr``    – mean expression across all cells (float)
        * ``pct_detected`` – fraction of cells with expression > 0 (float, 0–1)
        * ``n_units``      – number of cells used (int)

    Notes
    -----
    Originally duplicated verbatim in ``ras.py`` and ``rgs.py``.
    Consolidated here in v0.1.0.

    Examples
    --------
    >>> stats = compute_expr_stats_from_df(df)
    >>> stats.columns.tolist()
    ['symbol', 'mean_expr', 'pct_detected', 'n_units']
    """
    X = expr_df
    mean_expr    = X.mean(axis=0).astype(float)
    pct_detected = (X > 0).mean(axis=0).astype(float)
    n_units      = np.full(X.shape[1], X.shape[0], dtype=int)

    return pd.DataFrame({
        "symbol":       X.columns.astype(str).values,
        "mean_expr":    mean_expr.values,
        "pct_detected": pct_detected.values,
        "n_units":      n_units,
    })


# ──────────────────────────────────────────────────────────────
# 4. GMT file I/O
# ──────────────────────────────────────────────────────────────

def _open_maybe_gz(path: str, mode: str = "r"):
    """Open a plain or gzip-compressed text file transparently."""
    if path.lower().endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def read_gmt(
    path: str,
    as_sets: bool = False,
) -> Union[Dict[str, List[str]], Dict[str, Set[str]]]:
    """
    Read a GMT file and return a ``{name: genes}`` mapping.

    Handles both plain ``.gmt`` and gzip-compressed ``.gmt.gz`` files.
    Comment lines (``#``) and lines with fewer than 3 fields are skipped.
    Duplicate gene entries within a set are removed (first occurrence kept).

    Parameters
    ----------
    path : str
        Path to the GMT file (plain or ``.gz``).
    as_sets : bool, default ``False``
        If ``False`` (default), return ``Dict[str, List[str]]``
        (ordered, duplicate-free gene list — needed by ``rgs`` null matching).
        If ``True``, return ``Dict[str, Set[str]]``
        (unordered — used by ``ras`` AUCell loading).

    Returns
    -------
    dict
        ``{set_name: [gene1, gene2, ...]}`` or ``{set_name: {gene1, gene2, ...}}``
        depending on *as_sets*.

    Notes
    -----
    Previously three separate implementations existed:
    * ``_load_from_gmt`` in ``ras.py``  → Set output, .gz support
    * ``read_gmt``        in ``rgs.py``  → List output, no .gz support
    * ``read_gmt_lines``  in ``merge_regulons.py`` → raw lines (kept separate)

    Examples
    --------
    >>> regulons = read_gmt("regulons_symbol.gmt")
    >>> len(regulons)
    342
    >>> regulons_set = read_gmt("regulons_symbol.gmt", as_sets=True)
    """
    result: Dict[str, List[str]] = {}
    with _open_maybe_gz(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            name  = parts[0].strip()
            # deduplicate while preserving order
            seen: Set[str] = set()
            genes: List[str] = []
            for g in parts[2:]:
                g = g.strip()
                if g and g not in seen:
                    seen.add(g)
                    genes.append(g)
            if genes:
                result[name] = genes

    if as_sets:
        return {k: set(v) for k, v in result.items()}
    return result


def write_gmt(
    path: str,
    gene_sets: Dict[str, Union[List[str], Set[str]]],
    desc_suffix: str = "_regulon",
    sort_genes: bool = True,
    makedirs: bool = True,
) -> int:
    """
    Write a ``{name: genes}`` mapping to a GMT file.

    Parameters
    ----------
    path : str
        Output file path.
    gene_sets : dict
        ``{set_name: [gene, ...]}`` or ``{set_name: {gene, ...}}``.
        Sets with zero genes are skipped.
    desc_suffix : str, default ``"_regulon"``
        Appended to the set name to form the GMT description field
        (column 2).  For example, ``"HNRNPA1"`` → ``"HNRNPA1_regulon"``.
    sort_genes : bool, default ``True``
        Sort genes alphabetically within each set.
    makedirs : bool, default ``True``
        Create parent directories if they do not exist.

    Returns
    -------
    int
        Number of non-empty sets written.

    Notes
    -----
    Originally duplicated in ``get_regulon.py`` and ``rgs.py`` with
    minor differences (dir creation, sorting, return value).
    Consolidated here in v0.1.0.

    Examples
    --------
    >>> write_gmt("out/regulons.gmt", {"HNRNPA1": ["GAPDH", "ACTB"]})
    1
    """
    if makedirs:
        parent = os.path.dirname(path)
        if parent:
            os.makedirs(parent, exist_ok=True)

    kept = 0
    with open(path, "w", encoding="utf-8") as fo:
        for name, genes in sorted(gene_sets.items()):
            genes_list = sorted(genes) if sort_genes else list(genes)
            if not genes_list:
                continue
            line = [name, f"{name}{desc_suffix}"] + genes_list
            fo.write("\t".join(line) + "\n")
            kept += 1
    return kept
