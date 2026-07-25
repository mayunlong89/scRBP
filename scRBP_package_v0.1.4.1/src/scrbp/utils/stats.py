#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRBP shared statistical utilities
====================================
Shared statistical functions used across multiple scRBP pipeline steps:

  - assign_mode                  Regulatory mode labelling from correlation (getGRN, getMerge_GRN)
  - regulon_specificity_scores   Jensen-Shannon divergence-based RSS matrix (ras, trs)
"""

from __future__ import annotations

from typing import Union

import numpy as np
import pandas as pd


# ──────────────────────────────────────────────────────────────
# 1. Regulatory mode assignment
# ──────────────────────────────────────────────────────────────

def assign_mode(
    correlation: float,
    threshold: float = 0.0,
) -> str:
    """
    Assign a regulatory mode label based on an RBP–target correlation value.

    In ``getGRN``, the correlation is Spearman rank correlation (vectorised
    rank-based computation).  In ``getMerge_GRN``, it is the averaged
    Spearman correlation across multi-seed runs.

    Parameters
    ----------
    correlation : float
        Spearman correlation between an RBP and its target gene/isoform.
    threshold : float, default ``0.0``
        Minimum absolute correlation required to assign a directional label.
        Values with ``|correlation| <= threshold`` (or NaN) are labelled
        ``"unknown"``.

    Returns
    -------
    str
        One of ``"activating"``, ``"repressing"``, or ``"unknown"``.

    Notes
    -----
    Edges with Spearman |r| > threshold and positive correlation are labelled
    ``"activating"`` (RBP positively co-expressed with target); negative
    correlation edges are ``"repressing"``; edges at or below the threshold
    (including NaN) are ``"unknown"``.

    Examples
    --------
    >>> assign_mode(0.45)
    'activating'
    >>> assign_mode(-0.3)
    'repressing'
    >>> assign_mode(0.0)
    'unknown'
    >>> assign_mode(0.02, threshold=0.03)
    'unknown'
    >>> assign_mode(float('nan'))
    'unknown'
    """
    if pd.isna(correlation) or abs(correlation) <= threshold:
        return "unknown"
    return "activating" if correlation > 0 else "repressing"


# ──────────────────────────────────────────────────────────────
# 2. Regulon Specificity Score (RSS)
# ──────────────────────────────────────────────────────────────

def regulon_specificity_scores(
    auc_mtx: pd.DataFrame,
    cell_type_series: pd.Series,
    eps: float = 1e-12,
) -> pd.DataFrame:
    """
    Compute Regulon Specificity Scores (RSS) using Jensen-Shannon divergence.

    For each (regulon, cell-type) pair, the RSS quantifies how specifically a
    regulon is active in that cell type relative to all others.  It is defined
    as RSS = 1 − JSD(P_regulon ∥ Q_celltype), where:

    - P_regulon  is the AUCell score distribution normalised across all cells.
    - Q_celltype is the indicator (uniform) distribution over cells of that type.

    RSS = 1 indicates perfect cell-type specificity; RSS = 0 indicates no
    enrichment in the target cell type relative to the background.

    Parameters
    ----------
    auc_mtx : pd.DataFrame
        AUCell score matrix, shape ``(n_cells, n_regulons)``.
        Index must match *cell_type_series*.
    cell_type_series : pd.Series
        Cell-type label for every cell.  Index must align with *auc_mtx*.
        NaN cells are excluded from the computation.
    eps : float, default ``1e-12``
        Small constant added inside log to prevent log(0).

    Returns
    -------
    pd.DataFrame
        Shape ``(n_cell_types, n_regulons)``.
        Index = unique cell types, columns = regulon names.
        Regulons with all-zero AUCell scores receive NaN.

    Notes
    -----
    Uses vectorised NumPy operations throughout.  Regulons with all-zero
    AUCell scores across all cells receive NaN in the output matrix.

    Examples
    --------
    >>> rss = regulon_specificity_scores(aucell_df, obs["cell_type"])
    >>> rss.shape
    (12, 342)   # 12 cell types, 342 regulons
    >>> rss.idxmax(axis=0).head(3)
    HNRNPA1_regulon    Excitatory neurons
    PTBP1_regulon      Oligodendrocytes
    FUS_regulon        Microglia
    dtype: object
    """
    # Align indices
    if not auc_mtx.index.equals(cell_type_series.index):
        cell_type_series = cell_type_series.reindex(auc_mtx.index)

    # Remove NaN cells
    mask = ~cell_type_series.isna()
    if not mask.all():
        auc_mtx          = auc_mtx.loc[mask]
        cell_type_series = cell_type_series.loc[mask]

    A          = auc_mtx.to_numpy(dtype=np.float64)
    regulons   = auc_mtx.columns
    cell_types = pd.Index(pd.unique(cell_type_series), name="cell_type")
    labs       = cell_type_series.to_numpy()

    # Column-normalise: P[cell, reg] = AUC[cell, reg] / sum_cells(AUC[:, reg])
    col_sums  = A.sum(axis=0, keepdims=True)
    zero_cols = (col_sums.ravel() == 0.0)
    col_sums[col_sums == 0.0] = 1.0
    P = A / col_sums

    rss = np.empty((len(cell_types), len(regulons)), dtype=np.float32)

    for i, ct in enumerate(cell_types):
        in_ct = (labs == ct).astype(np.float64)
        n_ct  = in_ct.sum()
        if n_ct == 0:
            rss[i, :] = np.nan
            continue
        q      = (in_ct / n_ct)[:, np.newaxis]   # (n_cells, 1)
        Mmid   = 0.5 * (P + q)
        kl_pm  = (P * (np.log2(P + eps) - np.log2(Mmid + eps))).sum(axis=0)
        kl_qm  = (q * (np.log2(q + eps) - np.log2(Mmid + eps))).sum(axis=0)
        rss[i, :] = 1.0 - np.sqrt(0.5 * (kl_pm + kl_qm))

    # Mark all-zero regulons as NaN
    if zero_cols.any():
        rss[:, zero_cols] = np.nan

    return pd.DataFrame(rss, index=cell_types, columns=regulons)
