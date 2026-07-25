#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scRBP rgs_rare v8 — rare-variant regulon-level genetic association score

v8 change vs v7: the one-sided P used to derive RGS_z = norm.isf(P) is now
clipped on both tails (1e-300 .. 1-1e-15). This prevents RGS_z = -inf for
strongly depleted null regulons (P rounds to exactly 1.0), which previously
propagated -inf/NaN into null_mean, null_sd, z_empirical and P_param in the
ct-mode empirical audit. Empirical P (P_empirical/FDR_empirical) is unchanged.


This module mirrors the interface and output convention of ``scRBP rgs`` but
uses externally generated gene-level rare-variant evidence instead of MAGMA
common-variant gene results.

Main idea
---------
For each RBP regulon, scRBP rgs_rare fits a competitive gene-set regression:

    rare_score_g = beta0 + beta_s * I(g in regulon_s) + C_g^T gamma + error_g

where rare_score_g is a winsorized and z-standardized gene-level rare-variant
score.  For frequentist gene-level rare summaries, the score is -log10(P).  For
TADA-like Bayesian summaries, the score is logBF.  The regression covariate C_g
is kept purely genetic and captures coding rare-variant opportunity.  If
``--cds-scale raw`` is used, the raw union CDS length is log-transformed and
z-standardized within the analyzed gene universe.  If ``--cds-scale z_log`` is
used, the provided column is interpreted as a precomputed z-standardized
log-CDS covariate and used directly without additional transformation

Expression-derived covariates (mean_expr and pct_detected) are NOT included in
the competitive regression.  They are used only for covariate-matched null
regulon construction in ct mode, mirroring scRBP common-variant RGS where
expression covariates are used for null matching rather than MAGMA regression.
The regulon-level genetic association score is defined analogously to common
RGS as beta_s / SE(beta_s), with an optional -log10(P) score.

Modes
-----
- sc: test real regulons only and write a MAGMA-like RGS CSV.
- ct: construct covariate-matched null regulons, run the same rare competitive
      regression for REAL and NULL gene sets, and write RBP__REAL / RBP__NULL
      rows that are directly compatible with scRBP trs.py.

Outputs
-------
- <out>.gsa_RGS.csv
    Columns match common scRBP rgs output:
    sc: Regulon,GeneSet,NGENES,BETA,BETA_STD,SE,P,RGS_z,RGS_mlog10P
    ct: same + RBP,SET_KIND,NULL_ID

- <out>_REAL_PLUS_NULLS.<idtype>.gmt   [ct only]
    REAL+NULL regulons for recomputing RAS with ras.py.

- <out>_empirical.csv                  [ct only]
    Optional audit table summarizing REAL vs NULL RGS distributions.

Notes
-----
This script does not perform primary rare-variant association testing.  It
accepts gene-level rare-variant evidence from external frameworks such as
TADA/extTADA, SCHEMA-style burden tests, SAIGE-GENE+, REGENIE, or STAAR-O.
"""

from __future__ import annotations

import argparse
import gzip
import io
import os
import re
import sys
import tempfile
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import norm

try:
    from scrbp.utils.io import (
        read_gmt,
        write_gmt,
        load_expression_matrix,
        compute_expr_stats_from_df,
    )
except Exception:
    # Standalone fallback for development/testing outside the installed package.
    def read_gmt(path: str, as_sets: bool = False):
        opener = gzip.open if str(path).endswith(".gz") else open
        out = {}
        with opener(path, "rt") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                name = parts[0]
                genes = [g for g in parts[2:] if str(g).strip()]
                out[name] = set(genes) if as_sets else genes
        return out

    def write_gmt(path: str, sets: Dict[str, Iterable[str]]):
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "wt") as f:
            for name, genes in sets.items():
                f.write(str(name) + "\tna\t" + "\t".join(map(str, genes)) + "\n")

    def load_expression_matrix(path: str, dtype="float32"):
        df = pd.read_csv(path, index_col=0)
        # Match scRBP convention after loading: cells x genes.
        if df.shape[0] < df.shape[1]:
            df = df.T
        return df.astype(dtype)

    def compute_expr_stats_from_df(expr_df_cells_by_genes: pd.DataFrame):
        X = expr_df_cells_by_genes
        return pd.DataFrame({
            "symbol": X.columns.astype(str),
            "mean_expr": X.mean(axis=0).values,
            "pct_detected": (X > 0).mean(axis=0).values,
        })


# ------------------------- small utils -------------------------

def str2bool(v):
    if isinstance(v, bool):
        return v
    v = str(v).strip().lower()
    if v in ("y", "yes", "t", "true", "1"):
        return True
    if v in ("n", "no", "f", "false", "0"):
        return False
    raise argparse.ArgumentTypeError("Boolean value expected (True/False).")


def read_table_auto(path: str) -> pd.DataFrame:
    if not path or not os.path.exists(path):
        raise FileNotFoundError(path)
    low = path.lower()
    if low.endswith(".csv") or low.endswith(".csv.gz"):
        return pd.read_csv(path)
    if low.endswith(".tsv") or low.endswith(".txt") or low.endswith(".tsv.gz") or low.endswith(".txt.gz"):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path, sep=None, engine="python")


def qbin(x, q=10):
    x = pd.to_numeric(x, errors="coerce")
    if x.notna().sum() == 0:
        return pd.Series(np.zeros(len(x), dtype=int), index=x.index)
    try:
        qq = pd.qcut(x, q=q, labels=False, duplicates="drop")
        return qq.fillna(0).astype(int)
    except Exception:
        return pd.Series(np.zeros(len(x), dtype=int), index=x.index)


def zscore(x: pd.Series, eps: float = 1e-12) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    mu = x.mean(skipna=True)
    sd = x.std(skipna=True, ddof=1)
    if not np.isfinite(sd) or sd < eps:
        return pd.Series(np.nan, index=x.index)
    return (x - mu) / sd


def winsorize_top(x: pd.Series, top: float = 0.01) -> Tuple[pd.Series, float]:
    x = pd.to_numeric(x, errors="coerce")
    cutoff = float(x.quantile(1.0 - top))
    return pd.Series(np.minimum(x, cutoff), index=x.index), cutoff


def bh_fdr(pvals: Sequence[float]) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    pv = p[ok]
    m = len(pv)
    if m == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)
    tmp = np.empty(m, dtype=float)
    tmp[order] = q
    out[ok] = tmp
    return out


def by_fdr(pvals: Sequence[float]) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    pv = p[ok]
    m = len(pv)
    if m == 0:
        return out
    c_m = np.sum(1.0 / np.arange(1, m + 1))
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * m * c_m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)
    tmp = np.empty(m, dtype=float)
    tmp[order] = q
    out[ok] = tmp
    return out


def split_regulon_name(name: str):
    m_real = re.match(r"^(.*)__REAL$", str(name))
    if m_real:
        return m_real.group(1), "REAL", 0
    m_null = re.match(r"^(.*)__NULL_(\d+)$", str(name))
    if m_null:
        return m_null.group(1), "NULL", int(m_null.group(2))
    return str(name), "", None


# ------------------------- ID mapping -------------------------

def load_entrez_symbol_from_gene_loc(path: str):
    """Load MAGMA NCBI*.gene.loc: use columns 0=Entrez, 5=symbol."""
    df = pd.read_csv(path, sep=r"\s+", header=None, engine="python", usecols=[0, 5])
    df.columns = ["entrez", "symbol"]
    df = df.dropna()
    df = df[(df["symbol"] != "-") & (df["symbol"] != "")]
    df["entrez"] = df["entrez"].astype(int)
    e2s = dict(zip(df["entrez"], df["symbol"].astype(str)))
    s2e = {}
    for e, s in e2s.items():
        s2e.setdefault(str(s).upper(), int(e))
    return e2s, s2e


def _norm_symbol(x):
    return str(x).strip().upper()


def _norm_entrez(x):
    if pd.isna(x):
        return None
    try:
        return str(int(float(x)))
    except Exception:
        s = str(x).strip()
        if re.fullmatch(r"\d+", s):
            return s
        return None


def normalize_ids(values: Iterable, id_type: str) -> List[str]:
    if id_type == "symbol":
        return [_norm_symbol(v) for v in values if str(v).strip()]
    if id_type == "entrez":
        out = []
        for v in values:
            z = _norm_entrez(v)
            if z is not None:
                out.append(z)
        return out
    raise ValueError("id_type must be symbol or entrez")


def map_ids(values: Iterable, from_type: str, to_type: str,
            e2s: Optional[Dict[int, str]] = None,
            s2e: Optional[Dict[str, int]] = None) -> List[str]:
    if from_type == to_type:
        return normalize_ids(values, to_type)
    if from_type == "symbol" and to_type == "entrez":
        if s2e is None:
            raise SystemExit("--gene-loc is required to map symbol to Entrez")
        out = []
        for v in values:
            e = s2e.get(_norm_symbol(v))
            if e is not None:
                out.append(str(e))
        return out
    if from_type == "entrez" and to_type == "symbol":
        if e2s is None:
            raise SystemExit("--gene-loc is required to map Entrez to symbol")
        out = []
        for v in values:
            e = _norm_entrez(v)
            if e is None:
                continue
            s = e2s.get(int(e))
            if s:
                out.append(_norm_symbol(s))
        return out
    raise SystemExit("Unsupported ID mapping. Please provide symbol or Entrez IDs; convert ENSG upstream.")


def load_regulons_mapped(path: str, id_type: str, work_id_type: str,
                         e2s=None, s2e=None) -> Dict[str, List[str]]:
    raw = read_gmt(path)
    out = {}
    for name, genes in raw.items():
        mapped = map_ids(genes, id_type, work_id_type, e2s=e2s, s2e=s2e)
        mapped = list(dict.fromkeys(mapped))
        if mapped:
            out[str(name)] = mapped
    return out


def write_dual_gmt_if_possible(out_prefix: str, sets_work: Dict[str, List[str]],
                               work_id_type: str, e2s=None, s2e=None,
                               save_symbol: Optional[str] = None,
                               save_entrez: Optional[str] = None):
    if work_id_type == "symbol":
        sym_path = save_symbol if save_symbol else f"{out_prefix}_REAL_PLUS_NULLS.symbol.gmt"
        write_gmt(sym_path, sets_work)
        print(f"[SAVE] REAL+NULL Symbol GMT -> {sym_path}")
        if save_entrez:
            entrez_sets = {k: map_ids(v, "symbol", "entrez", e2s=e2s, s2e=s2e) for k, v in sets_work.items()}
            entrez_sets = {k: v for k, v in entrez_sets.items() if v}
            write_gmt(save_entrez, entrez_sets)
            print(f"[SAVE] REAL+NULL Entrez GMT -> {save_entrez}")
    else:
        ent_path = save_entrez if save_entrez else f"{out_prefix}_REAL_PLUS_NULLS.entrez.gmt"
        write_gmt(ent_path, sets_work)
        print(f"[SAVE] REAL+NULL Entrez GMT -> {ent_path}")
        if save_symbol:
            symbol_sets = {k: map_ids(v, "entrez", "symbol", e2s=e2s, s2e=s2e) for k, v in sets_work.items()}
            symbol_sets = {k: v for k, v in symbol_sets.items() if v}
            write_gmt(save_symbol, symbol_sets)
            print(f"[SAVE] REAL+NULL Symbol GMT -> {save_symbol}")


# ------------------------- expression and CDS covariates -------------------------

def read_expr_stats(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")
    df = df.rename(columns={c: c.lower() for c in df.columns})
    need = {"symbol", "mean_expr", "pct_detected"}
    if not need.issubset(df.columns):
        raise ValueError(f"expr-stats file must contain columns: symbol, mean_expr, pct_detected; got {list(df.columns)}")
    df["symbol"] = df["symbol"].astype(str).str.upper()
    df["mean_expr"] = pd.to_numeric(df["mean_expr"], errors="coerce")
    df["pct_detected"] = pd.to_numeric(df["pct_detected"], errors="coerce")
    return df[["symbol", "mean_expr", "pct_detected"]].drop_duplicates("symbol")


def compute_or_load_expr_stats(args) -> Optional[pd.DataFrame]:
    if args.expr_stats and os.path.isfile(args.expr_stats):
        df = read_expr_stats(args.expr_stats)
        print(f"[ExprStats] Loaded: {args.expr_stats} | genes={df.shape[0]}")
        return df
    if args.emit_expr_stats:
        if not args.matrix_stats:
            raise SystemExit("--emit-expr-stats requires --matrix-stats")
        X = load_expression_matrix(args.matrix_stats, dtype=args.dtype)
        df = compute_expr_stats_from_df(X)
        df = df.rename(columns={c: c.lower() for c in df.columns})
        if "gene" in df.columns and "symbol" not in df.columns:
            df = df.rename(columns={"gene": "symbol"})
        df["symbol"] = df["symbol"].astype(str).str.upper()
        outp = args.expr_stats_out if args.expr_stats_out else f"{args.out}_expr_stats.tsv"
        df[["symbol", "mean_expr", "pct_detected"]].to_csv(outp, sep="\t", index=False)
        print(f"[ExprStats] Computed & saved: {outp} | genes={df.shape[0]}")
        return df[["symbol", "mean_expr", "pct_detected"]]
    return None


def attach_expr_stats(gene_df: pd.DataFrame, expr_df: Optional[pd.DataFrame],
                      work_id_type: str, e2s=None, s2e=None,
                      allow_missing: bool = False) -> pd.DataFrame:
    df = gene_df.copy()
    if expr_df is None:
        if allow_missing:
            print("[ExprStats][WARN] Missing expression covariates; setting mean_expr=pct_detected=0.")
            df["mean_expr"] = 0.0
            df["pct_detected"] = 0.0
            return df
        raise SystemExit("--expr-stats is required for ct-mode rare matched-null construction; or use --emit-expr-stats / --allow-missing-expr-stats.")

    expr = expr_df.copy()
    if work_id_type == "symbol":
        expr["gene"] = expr["symbol"].astype(str).str.upper()
    else:
        if s2e is None:
            raise SystemExit("--gene-loc is required to map expr-stats symbols to Entrez")
        expr["gene"] = expr["symbol"].map(lambda x: str(s2e.get(str(x).upper())) if s2e.get(str(x).upper()) is not None else None)
    expr = expr.dropna(subset=["gene"])[["gene", "mean_expr", "pct_detected"]].drop_duplicates("gene")
    df = df.merge(expr, on="gene", how="left")
    for c in ("mean_expr", "pct_detected"):
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)
    return df


def attach_cds_length(gene_df: pd.DataFrame, rare_df_original: pd.DataFrame, args,
                      work_id_type: str, e2s=None, s2e=None) -> pd.DataFrame:
    """Attach the coding-opportunity covariate.

    --cds-scale raw:
        args.cds_col is raw union CDS length.  The script computes
        log1p(raw length) and z-standardizes it within the analyzed gene universe.

    --cds-scale z_log:
        args.cds_col is already a z-standardized log-CDS covariate.  The script
        uses it directly as z_log_union_CDS_length, with no additional log or
        z-score transformation.  This avoids double-standardizing precomputed
        z_log_union_cds columns from external rare-variant summaries.
    """
    df = gene_df.copy()
    if args.cds_length:
        cds = read_table_auto(args.cds_length)
        if args.cds_gene_col not in cds.columns:
            raise SystemExit(f"--cds-gene-col {args.cds_gene_col!r} not found in --cds-length")
        if args.cds_col not in cds.columns:
            raise SystemExit(f"--cds-col {args.cds_col!r} not found in --cds-length")
        mapped = []
        for v in cds[args.cds_gene_col].tolist():
            m = map_ids([v], args.cds_id_type, work_id_type, e2s=e2s, s2e=s2e)
            mapped.append(m[0] if m else None)
        cds["gene"] = mapped
        cds = cds.dropna(subset=["gene"])[["gene", args.cds_col]].copy()
        cds = cds.rename(columns={args.cds_col: "cds_input_value"})
        cds["cds_input_value"] = pd.to_numeric(cds["cds_input_value"], errors="coerce")
        cds = cds.sort_values("cds_input_value", ascending=False).drop_duplicates("gene")
        df = df.merge(cds, on="gene", how="left")
    else:
        # Use CDS covariate already present in rare summary after gene harmonization.
        if args.cds_col not in rare_df_original.columns:
            raise SystemExit(f"No --cds-length provided and --cds-col {args.cds_col!r} not found in rare summary")
        tmp = rare_df_original[[args.rare_gene_col, args.cds_col]].copy()
        mapped = []
        for v in tmp[args.rare_gene_col].tolist():
            m = map_ids([v], args.rare_id_type, work_id_type, e2s=e2s, s2e=s2e)
            mapped.append(m[0] if m else None)
        tmp["gene"] = mapped
        tmp = tmp.dropna(subset=["gene"])[["gene", args.cds_col]]
        tmp = tmp.rename(columns={args.cds_col: "cds_input_value"})
        tmp["cds_input_value"] = pd.to_numeric(tmp["cds_input_value"], errors="coerce")
        tmp = tmp.sort_values("cds_input_value", ascending=False).drop_duplicates("gene")
        df = df.merge(tmp, on="gene", how="left")

    df["cds_input_value"] = pd.to_numeric(df["cds_input_value"], errors="coerce")
    if df["cds_input_value"].isna().all():
        raise SystemExit("All CDS covariate values are missing after merge. Check gene IDs and --cds-* arguments.")

    if args.cds_scale == "raw":
        med = df["cds_input_value"].median(skipna=True)
        df["union_CDS_length"] = df["cds_input_value"].fillna(med).clip(lower=1.0)
        df["log_union_CDS_length"] = np.log1p(df["union_CDS_length"])
        df["z_log_union_CDS_length"] = zscore(df["log_union_CDS_length"])
    elif args.cds_scale == "z_log":
        med = df["cds_input_value"].median(skipna=True)
        df["union_CDS_length"] = np.nan
        df["log_union_CDS_length"] = np.nan
        df["z_log_union_CDS_length"] = df["cds_input_value"].fillna(med)
    else:
        raise SystemExit("--cds-scale must be raw or z_log")

    return df


# ------------------------- rare gene scores -------------------------

def make_rare_gene_table(args, work_id_type: str, e2s=None, s2e=None) -> Tuple[pd.DataFrame, float, pd.DataFrame]:
    rare0 = read_table_auto(args.rare_summary)
    if args.rare_gene_col not in rare0.columns:
        raise SystemExit(f"--rare-gene-col {args.rare_gene_col!r} not found in rare summary")

    mapped = []
    for v in rare0[args.rare_gene_col].tolist():
        m = map_ids([v], args.rare_id_type, work_id_type, e2s=e2s, s2e=s2e)
        mapped.append(m[0] if m else None)
    rare = rare0.copy()
    rare["gene"] = mapped
    rare = rare.dropna(subset=["gene"]).copy()

    if args.score_mode == "pvalue":
        if args.p_col not in rare.columns:
            raise SystemExit(f"--p-col {args.p_col!r} not found in rare summary")
        p = pd.to_numeric(rare[args.p_col], errors="coerce").clip(lower=1e-300, upper=1.0)
        rare["rare_score_raw"] = -np.log10(p)
    elif args.score_mode == "logbf":
        if args.logbf_col not in rare.columns:
            raise SystemExit(f"--logbf-col {args.logbf_col!r} not found in rare summary")
        rare["rare_score_raw"] = pd.to_numeric(rare[args.logbf_col], errors="coerce")
    elif args.score_mode == "direct":
        if args.score_col not in rare.columns:
            raise SystemExit(f"--score-col {args.score_col!r} not found in rare summary")
        rare["rare_score_raw"] = pd.to_numeric(rare[args.score_col], errors="coerce")
    else:
        raise SystemExit("--score-mode must be pvalue, logbf, or direct")

    rare = rare.dropna(subset=["rare_score_raw"]).copy()
    rare["rare_score_win"], cutoff = winsorize_top(rare["rare_score_raw"], top=args.top_winsor)
    rare["rare_score_z"] = zscore(rare["rare_score_win"])

    # Deduplicate genes by strongest winsorized evidence.
    rare = (rare.sort_values("rare_score_win", ascending=False)
                .drop_duplicates("gene", keep="first")
                .reset_index(drop=True))

    gene_df = rare[["gene", "rare_score_raw", "rare_score_win", "rare_score_z"]].copy()
    return gene_df, cutoff, rare0


# ------------------------- competitive regression -------------------------

def ols_hc3_regression(y: np.ndarray, membership: np.ndarray, covars: np.ndarray):
    """OLS y ~ 1 + membership + covars with HC3 robust SE for membership."""
    y = np.asarray(y, dtype=float)
    membership = np.asarray(membership, dtype=float).reshape(-1, 1)
    covars = np.asarray(covars, dtype=float)
    if covars.ndim == 1:
        covars = covars.reshape(-1, 1)

    X = np.column_stack([np.ones(len(y)), membership, covars])
    ok = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y = y[ok]
    X = X[ok]
    n, p = X.shape
    if n <= p + 5:
        raise ValueError("not enough observations")
    if np.nanstd(X[:, 1]) < 1e-12:
        raise ValueError("membership has no variation")

    XtX_inv = np.linalg.pinv(X.T @ X)
    beta_hat = XtX_inv @ X.T @ y
    resid = y - X @ beta_hat
    h = np.sum((X @ XtX_inv) * X, axis=1)
    h = np.clip(h, 0, 1 - 1e-8)
    w = (resid / (1.0 - h)) ** 2
    meat = X.T @ (X * w[:, None])
    vcov = XtX_inv @ meat @ XtX_inv
    se = float(np.sqrt(max(vcov[1, 1], 0.0)))
    beta = float(beta_hat[1])
    if not np.isfinite(se) or se <= 0:
        raise ValueError("invalid SE")
    z = beta / se
    p_one = float(norm.sf(z))  # one-sided enrichment, beta > 0

    # Standardized beta, analogous in spirit to MAGMA BETA_STD.
    x = X[:, 1]
    sd_x = np.nanstd(x, ddof=1)
    sd_y = np.nanstd(y, ddof=1)
    beta_std = float(beta * sd_x / sd_y) if sd_y > 0 else np.nan
    return beta, beta_std, se, p_one, z


def run_competitive_rgs_for_sets(gene_df: pd.DataFrame,
                                 sets: Dict[str, List[str]],
                                 covar_cols: Sequence[str],
                                 min_genes: int = 0,
                                 max_regulon_frac: float = 0.5) -> pd.DataFrame:
    need = ["gene", "rare_score_z"] + list(covar_cols)
    missing = [c for c in need if c not in gene_df.columns]
    if missing:
        raise SystemExit(f"Missing columns for regression: {missing}")

    df = gene_df[need].dropna().copy()
    universe = set(df["gene"].astype(str))
    genes = df["gene"].astype(str).tolist()
    y = df["rare_score_z"].astype(float).values
    C = df[list(covar_cols)].astype(float).values
    n_universe = len(genes)

    rows = []
    for set_name, geneset in sets.items():
        targets = list(dict.fromkeys([str(g) for g in geneset if str(g) in universe]))
        n_targets = len(targets)
        if min_genes and n_targets < int(min_genes):
            continue
        if n_targets == 0:
            continue
        if n_targets / max(n_universe, 1) > max_regulon_frac:
            continue

        target_set = set(targets)
        membership = np.array([g in target_set for g in genes], dtype=int)
        try:
            beta, beta_std, se, p_one, z = ols_hc3_regression(y, membership, C)
        except Exception as e:
            # Keep failed sets out of the final table, matching MAGMA behavior better.
            continue

        # Clip the one-sided P on BOTH tails before deriving RGS_z = norm.isf(P):
        #   - lower floor 1e-300 keeps strongly enriched sets from RGS_z = +inf;
        #   - upper cap (1 - 1e-15) keeps strongly DEPLETED null sets (beta << 0,
        #     P -> 1.0) from RGS_z = norm.isf(1.0) = -inf. Depletion is not of
        #     interest for the one-sided enrichment test, so capping the negative
        #     tail at ~ -8 is harmless and keeps null_mean/null_sd/z_emp finite.
        p_safe = min(max(float(p_one), 1e-300), 1.0 - 1e-15)
        rows.append({
            "Regulon": set_name,
            "GeneSet": "RBP_regulon",
            "NGENES": int(n_targets),
            "BETA": beta,
            "BETA_STD": beta_std,
            "SE": se,
            "P": p_safe,
            "RGS_z": float(norm.isf(p_safe)),
            "RGS_mlog10P": float(-np.log10(p_safe)),
        })

    out = pd.DataFrame(rows)
    if len(out) == 0:
        raise SystemExit("No regulons were tested. Check gene ID overlap and --min_genes.")
    return out


# ------------------------- matched null regulons -------------------------

def add_rare_null_bins(bg_df: pd.DataFrame, q: int = 5) -> pd.DataFrame:
    bg = bg_df.copy()
    # Use the same coding-opportunity covariate used in rare RGS regression.
    # With --cds-scale raw this is internally computed from log1p(raw CDS);
    # with --cds-scale z_log this is the precomputed z_log covariate.
    bg["bin_cds"] = qbin(bg["z_log_union_CDS_length"], q=q)
    bg["bin_mean"] = qbin(bg["mean_expr"], q=q)
    bg["bin_pct"] = qbin(bg["pct_detected"], q=q)
    return bg


def _bucket_dict(pool: pd.DataFrame, cols: Sequence[str]) -> Dict[Tuple, np.ndarray]:
    out = {}
    for key, sub in pool.groupby(list(cols), dropna=False):
        if not isinstance(key, tuple):
            key = (key,)
        out[key] = sub["gene"].astype(str).values
    return out


def _sample_from_bucket(rng, candidates: np.ndarray, cnt: int, replace_if_needed: bool = True):
    if len(candidates) == 0:
        return []
    replace = len(candidates) < cnt if replace_if_needed else False
    return rng.choice(candidates, size=cnt, replace=replace).tolist()


def make_null_for_regulon_3d(real_genes: Sequence[str],
                             bg_binned: pd.DataFrame,
                             n_null: int = 1000,
                             rng=None,
                             exclude_self: bool = True,
                             min_bucket_size: int = 5) -> List[List[str]]:
    """Construct null regulons matched on union-CDS length, mean expression and pct detected."""
    if rng is None:
        rng = np.random.default_rng(1)

    real_genes = list(dict.fromkeys([str(g) for g in real_genes]))
    reg = bg_binned[bg_binned["gene"].isin(real_genes)].copy()
    if reg.empty:
        return []

    pool = bg_binned[~bg_binned["gene"].isin(real_genes)].copy() if exclude_self else bg_binned.copy()
    if pool.empty:
        return []

    need3 = reg.groupby(["bin_cds", "bin_mean", "bin_pct"]).size().to_dict()
    need2 = reg.groupby(["bin_cds", "bin_mean"]).size().to_dict()
    need1 = reg.groupby(["bin_cds"]).size().to_dict()

    buckets3 = _bucket_dict(pool, ["bin_cds", "bin_mean", "bin_pct"])
    buckets2 = _bucket_dict(pool, ["bin_cds", "bin_mean"])
    buckets1 = _bucket_dict(pool, ["bin_cds"])
    all_pool = pool["gene"].astype(str).values

    null_sets = []
    for _ in range(n_null):
        picked = []
        for k3, cnt in need3.items():
            if not isinstance(k3, tuple):
                k3 = (k3,)
            cand = buckets3.get(k3, np.array([], dtype=object))
            if len(cand) >= max(cnt, min_bucket_size):
                picked.extend(_sample_from_bucket(rng, cand, cnt, replace_if_needed=True))
                continue
            # 2D fallback: CDS x mean expression
            k2 = (k3[0], k3[1])
            cand = buckets2.get(k2, np.array([], dtype=object))
            if len(cand) >= max(cnt, min_bucket_size):
                picked.extend(_sample_from_bucket(rng, cand, cnt, replace_if_needed=True))
                continue
            # 1D fallback: CDS length only
            k1 = (k3[0],)
            cand = buckets1.get(k1, np.array([], dtype=object))
            if len(cand) >= max(cnt, min_bucket_size):
                picked.extend(_sample_from_bucket(rng, cand, cnt, replace_if_needed=True))
                continue
            # global fallback
            picked.extend(_sample_from_bucket(rng, all_pool, cnt, replace_if_needed=True))

        # Preserve the null regulon size as closely as possible.  If replacement
        # sampling introduced duplicates in sparse buckets, fill with additional
        # unique genes from the global matched-universe pool.  This mirrors the
        # common RGS expectation that REAL and NULL sets have the same nominal
        # regulon size before TRS consumes the REAL+NULL outputs.
        target_size = len(real_genes)
        picked = list(dict.fromkeys(map(str, picked)))
        if len(picked) < target_size:
            remaining = np.array([g for g in all_pool if str(g) not in set(picked)], dtype=object)
            need_more = target_size - len(picked)
            if len(remaining) >= need_more:
                picked.extend(rng.choice(remaining, size=need_more, replace=False).tolist())
            elif len(remaining) > 0:
                picked.extend(remaining.tolist())
        null_sets.append(picked[:target_size])
    return null_sets


def build_real_plus_null_sets(real_sets: Dict[str, List[str]],
                              bg_df: pd.DataFrame,
                              n_null: int,
                              q_bins: int,
                              seed: int,
                              exclude_self: bool,
                              min_bucket_size: int,
                              min_genes: int = 0) -> Tuple[Dict[str, List[str]], pd.DataFrame]:
    rng = np.random.default_rng(seed)
    bg_binned = add_rare_null_bins(bg_df, q=q_bins)
    universe = set(bg_binned["gene"].astype(str))

    sets_out = {}
    index_rows = []
    for rbp, genes in real_sets.items():
        real = [g for g in list(dict.fromkeys(map(str, genes))) if g in universe]
        if min_genes and len(real) < int(min_genes):
            continue
        if len(real) == 0:
            continue
        real_name = f"{rbp}__REAL"
        sets_out[real_name] = real
        index_rows.append((rbp, 0, real_name, len(real), "REAL"))
        nulls = make_null_for_regulon_3d(
            real, bg_binned, n_null=n_null, rng=rng,
            exclude_self=exclude_self, min_bucket_size=min_bucket_size,
        )
        for i, nl in enumerate(nulls, 1):
            nm = f"{rbp}__NULL_{i:04d}"
            sets_out[nm] = [str(x) for x in nl]
            index_rows.append((rbp, i, nm, len(nl), "NULL"))

    idx = pd.DataFrame(index_rows, columns=["RBP", "NULL_ID", "SET_NAME", "SIZE", "SET_KIND"])
    return sets_out, idx


# ------------------------- empirical audit -------------------------

def add_ct_split_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    splits = out["Regulon"].apply(split_regulon_name)
    out["RBP"] = splits.apply(lambda x: x[0])
    out["SET_KIND"] = splits.apply(lambda x: x[1])
    out["NULL_ID"] = splits.apply(lambda x: x[2])
    return out


def empirical_summary_from_rgs(ct_df: pd.DataFrame, score_col: str = "RGS_z",
                               fdr_method: str = "BH") -> pd.DataFrame:
    rows = []
    df = ct_df.copy()
    if "RBP" not in df.columns:
        df = add_ct_split_columns(df)
    for rbp, sub in df.groupby("RBP", sort=False):
        real = sub[sub["SET_KIND"] == "REAL"]
        null = sub[sub["SET_KIND"] == "NULL"].sort_values("NULL_ID")
        if real.empty or null.empty:
            continue
        x = float(real.iloc[0][score_col])
        znull = null[score_col].astype(float).to_numpy()
        mu = float(np.nanmean(znull))
        sd = float(np.nanstd(znull, ddof=1))
        if not np.isfinite(sd) or sd <= 1e-12:
            z_emp = np.nan
            p_param = np.nan
        else:
            z_emp = (x - mu) / sd
            p_param = float(norm.sf(z_emp))
        p_emp = float((1.0 + np.sum(znull >= x)) / (len(znull) + 1.0))
        rows.append({
            "RBP": rbp,
            "Regulon": str(real.iloc[0]["Regulon"]),
            "NGENES": int(real.iloc[0]["NGENES"]),
            "RGS_real": x,
            "null_mean": mu,
            "null_sd": sd,
            "P_empirical": p_emp,
            "z_empirical": z_emp,
            "P_param": p_param,
            "K_null": int(len(znull)),
        })
    out = pd.DataFrame(rows)
    if len(out) > 0:
        out["FDR_empirical"] = bh_fdr(out["P_empirical"].values) if fdr_method == "BH" else by_fdr(out["P_empirical"].values)
        out["FDR_param"] = bh_fdr(out["P_param"].values) if fdr_method == "BH" else by_fdr(out["P_param"].values)
    return out


# ------------------------- setup data -------------------------

def prepare_gene_covariate_table(args, work_id_type: str, e2s=None, s2e=None) -> Tuple[pd.DataFrame, float]:
    """Prepare gene-level rare score and covariates.

    Regression covariates are intentionally kept purely genetic:
        z_log_union_CDS_length

    Expression covariates are attached only for ct-mode matched-null construction
    or when explicitly provided/requested.  This mirrors common scRBP RGS: genetic
    regression first, expression-aware null matching second.
    """
    gene_df, cutoff, rare_original = make_rare_gene_table(args, work_id_type, e2s=e2s, s2e=s2e)
    gene_df = attach_cds_length(gene_df, rare_original, args, work_id_type, e2s=e2s, s2e=s2e)

    # z_log_union_CDS_length is the only default regression covariate.
    gene_df["z_log_union_CDS_length"] = pd.to_numeric(
        gene_df["z_log_union_CDS_length"], errors="coerce"
    ).fillna(0.0)

    # Expression covariates are required for ct-mode matched-null construction,
    # but not for sc-mode competitive RGS regression.
    need_expr_for_null = (args.mode == "ct")
    expr = compute_or_load_expr_stats(args)
    if need_expr_for_null or expr is not None or args.allow_missing_expr_stats:
        gene_df = attach_expr_stats(
            gene_df, expr, work_id_type, e2s=e2s, s2e=s2e,
            allow_missing=args.allow_missing_expr_stats
        )
    else:
        gene_df["mean_expr"] = 0.0
        gene_df["pct_detected"] = 0.0

    # These z-scored expression fields are retained for audit/downstream use and
    # optional sensitivity analyses, but are not part of the default regression.
    gene_df["z_mean_expr"] = zscore(np.log1p(gene_df["mean_expr"].clip(lower=0))).fillna(0.0)
    gene_df["z_pct_detected"] = zscore(gene_df["pct_detected"].clip(lower=0, upper=1)).fillna(0.0)

    for c in ["z_log_union_CDS_length", "mean_expr", "pct_detected", "z_mean_expr", "z_pct_detected"]:
        gene_df[c] = pd.to_numeric(gene_df[c], errors="coerce").fillna(0.0)
    return gene_df, cutoff


# ------------------------- CLI -------------------------

def register_subcommand(subparsers):
    ap = subparsers.add_parser(
        "rgs_rare",
        help="Rare-variant regulon genetic association score",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )

    ap.add_argument("--mode", choices=["sc", "ct"], required=True,
                    help="sc: real regulons only; ct: add matched null regulons per real set.")
    ap.add_argument("--rare-summary", required=True, help="Gene-level rare-variant summary CSV/TSV.")
    ap.add_argument("--rare-gene-col", required=True, help="Gene column in rare summary.")
    ap.add_argument("--rare-id-type", choices=["symbol", "entrez"], default="symbol")
    ap.add_argument("--score-mode", choices=["pvalue", "logbf", "direct"], required=True)
    ap.add_argument("--p-col", default=None, help="P-value column for --score-mode pvalue.")
    ap.add_argument("--logbf-col", default=None, help="logBF column for --score-mode logbf.")
    ap.add_argument("--score-col", default=None, help="Score column for --score-mode direct.")
    ap.add_argument("--top-winsor", type=float, default=0.01,
                    help="Upper-tail winsorization fraction for gene-level rare scores.")

    ap.add_argument("--sets", required=True, help="Input regulon GMT.")
    ap.add_argument("--id-type", choices=["symbol", "entrez"], default="symbol",
                    help="Gene ID type in --sets.")
    ap.add_argument("--out", required=True, help="Output prefix.")
    ap.add_argument("--min_genes", "--min-genes", dest="min_genes", type=int, default=0)
    ap.add_argument("--max-regulon-frac", type=float, default=0.5)

    # ID mapping
    ap.add_argument("--gene-loc", help="MAGMA NCBI*.gene.loc for Entrez<->Symbol mapping.")

    # CDS length
    ap.add_argument("--cds-length", default=None, help="Optional gene CDS length table CSV/TSV.")
    ap.add_argument("--cds-gene-col", default="symbol", help="Gene column in --cds-length.")
    ap.add_argument("--cds-id-type", choices=["symbol", "entrez"], default="symbol")
    ap.add_argument("--cds-col", default="union_cds_length",
                    help="CDS covariate column; searched in rare summary if --cds-length not provided.")
    ap.add_argument("--cds-scale", choices=["raw", "z_log"], default="raw",
                    help="Interpretation of --cds-col. raw: raw union CDS length, internally log1p-transformed and z-scored; z_log: precomputed z-standardized log-CDS covariate, used directly without additional transformation.")

    # Expr-stats, matching ras.py/rgs.py style
    ap.add_argument("--expr-stats", help="Precomputed expr stats TSV/CSV with symbol, mean_expr, pct_detected. Required for --mode ct matched-null construction.")
    ap.add_argument("--emit-expr-stats", type=str2bool, nargs="?", const=True, default=False,
                    help="If True and --expr-stats not provided, compute from --matrix-stats and save for ct-mode null matching.")
    ap.add_argument("--matrix-stats", help="Expression matrix to compute expr-stats when --emit-expr-stats is True.")
    ap.add_argument("--expr-stats-out", help="Where to save computed expr-stats (default: <out>_expr_stats.tsv).")
    ap.add_argument("--dtype", default="float32", help="dtype for --matrix-stats loading.")
    ap.add_argument("--allow-missing-expr-stats", action="store_true",
                    help="Allow missing mean_expr/pct_detected by setting them to zero for ct-mode null matching. Not recommended.")

    # ct-specific null controls
    ap.add_argument("--n-null", "--null", dest="n_null", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=2025)
    ap.add_argument("--q-bins", type=int, default=5,
                    help="Quantile bins for rare matched nulls. Use 10 for stricter matching if universe is large.")
    ap.add_argument("--exclude-self", type=str2bool, nargs="?", const=True, default=True)
    ap.add_argument("--min-bucket-size", type=int, default=5)
    ap.add_argument("--fdr-method", choices=["BH", "BY"], default="BH")

    # Save REAL+NULL GMT controls
    ap.add_argument("--save-null-gmt", help="Save REAL+NULL GMT in working ID type (default: <out>_REAL_PLUS_NULLS.<idtype>.gmt).")
    ap.add_argument("--save-null-gmt-symbol", help="Also save REAL+NULL Symbol GMT (requires --gene-loc if working ID is entrez).")
    ap.add_argument("--save-null-gmt-entrez", help="Also save REAL+NULL Entrez GMT (requires --gene-loc if working ID is symbol).")

    ap.set_defaults(func=main)
    return ap


def main(args):
    # Determine internal working ID type from input regulons, matching scRBP rgs behavior.
    work_id_type = args.id_type

    e2s = s2e = None
    if args.gene_loc:
        e2s, s2e = load_entrez_symbol_from_gene_loc(args.gene_loc)
    elif args.rare_id_type != work_id_type or args.cds_id_type != work_id_type or work_id_type == "entrez":
        # Entrez-only inputs do not strictly need gene-loc unless expression stats must be mapped.
        pass

    gene_df, cutoff = prepare_gene_covariate_table(args, work_id_type, e2s=e2s, s2e=s2e)
    score_path = f"{args.out}_gene_scores.tsv"
    os.makedirs(os.path.dirname(score_path) or ".", exist_ok=True)
    gene_df.to_csv(score_path, sep="\t", index=False)
    print(f"[GeneScores] Wrote: {score_path} | genes={gene_df.shape[0]} | winsor_cutoff={cutoff:.6g}")

    real_sets = load_regulons_mapped(args.sets, args.id_type, work_id_type, e2s=e2s, s2e=s2e)
    universe = set(gene_df["gene"].astype(str))
    real_sets = {k: [g for g in v if g in universe] for k, v in real_sets.items()}
    real_sets = {k: v for k, v in real_sets.items() if len(v) >= int(args.min_genes)}
    if not real_sets:
        raise SystemExit("No regulons overlap the rare gene universe. Check IDs.")

    covar_cols = ["z_log_union_CDS_length"]

    if args.mode == "sc":
        # sc mode: real regulons only, no __REAL suffix, same as common sc RGS.
        res = run_competitive_rgs_for_sets(
            gene_df, real_sets, covar_cols=covar_cols,
            min_genes=args.min_genes, max_regulon_frac=args.max_regulon_frac,
        )
        out_csv = f"{args.out}.gsa_RGS.csv"
        res = res[["Regulon", "GeneSet", "NGENES", "BETA", "BETA_STD", "SE", "P", "RGS_z", "RGS_mlog10P"]]
        res.to_csv(out_csv, index=False)
        print("Wrote:", out_csv)
        return

    # ct mode: build REAL+NULL sets and run the exact same rare competitive regression on all sets.
    sets_real_null, null_index = build_real_plus_null_sets(
        real_sets, gene_df, n_null=args.n_null, q_bins=args.q_bins, seed=args.seed,
        exclude_self=args.exclude_self, min_bucket_size=args.min_bucket_size,
        min_genes=args.min_genes,
    )
    if not sets_real_null:
        raise SystemExit("Failed to build REAL+NULL regulons.")

    null_index_path = f"{args.out}.null_index.tsv"
    null_index.to_csv(null_index_path, sep="\t", index=False)
    print("Wrote:", null_index_path)

    # Save REAL+NULL GMT for ras.py.
    if args.save_null_gmt:
        write_gmt(args.save_null_gmt, sets_real_null)
        print(f"[SAVE] REAL+NULL GMT -> {args.save_null_gmt}")
    else:
        write_dual_gmt_if_possible(
            args.out, sets_real_null, work_id_type,
            e2s=e2s, s2e=s2e,
            save_symbol=args.save_null_gmt_symbol,
            save_entrez=args.save_null_gmt_entrez,
        )

    res = run_competitive_rgs_for_sets(
        gene_df, sets_real_null, covar_cols=covar_cols,
        min_genes=args.min_genes, max_regulon_frac=args.max_regulon_frac,
    )
    res = add_ct_split_columns(res)
    res = res[[
        "Regulon", "GeneSet", "NGENES", "BETA", "BETA_STD", "SE", "P",
        "RGS_z", "RGS_mlog10P", "RBP", "SET_KIND", "NULL_ID",
    ]]
    out_csv = f"{args.out}.gsa_RGS.csv"
    res.to_csv(out_csv, index=False)
    print("Wrote:", out_csv)

    emp = empirical_summary_from_rgs(res, score_col="RGS_z", fdr_method=args.fdr_method)
    emp_out = f"{args.out}_empirical.csv"
    emp.to_csv(emp_out, index=False)
    print("Wrote:", emp_out)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    subp = p.add_subparsers()
    register_subcommand(subp)
    args = p.parse_args()
    if not hasattr(args, "func"):
        p.print_help()
        raise SystemExit(1)
    args.func(args)
