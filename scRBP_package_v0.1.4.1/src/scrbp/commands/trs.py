#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scRBP trs — Trait-Relevant Score (TRS): integrative RAS × RGS prioritisation

Integrates the Regulon Activity Score (RAS) and the Regulon-level Genetic
Association Score (RGS) into a unified Trait-Relevant Score (TRS) that
prioritises RBP regulons that are both transcriptionally active in a given
cellular context and genetically enriched for a complex trait of interest.

TRS formula
-----------
    TRS = scaled_RAS + scaled_RGS - lambda × |scaled_RAS - scaled_RGS| / sqrt(2)

where:
  - scaled_RAS and scaled_RGS are min-max normalised scores (capped at the
    99th percentile, then rescaled to [0, 1]).
  - lambda is a tunable regularisation parameter (default = 1) that controls
    the discordance penalty.
  - The penalty term (|RAS - RGS| / sqrt(2)) is the sample standard deviation
    of the two scores, and increases monotonically with their discrepancy.

Regulons with concordantly high RAS and RGS achieve the highest TRS values,
while regulons with strong activity but weak genetic support (or vice versa)
are down-weighted by the discordance penalty.

Statistical significance
------------------------
Empirical p-values (p_RAS, p_RGS, p_TRS) and z-scores are computed by
Monte Carlo sampling against matched null regulons (preserving regulon size
and genomic complexity), with per-RBP K_null null sets.  A regulon is
considered significant at FDR < 0.05 on all three metrics simultaneously.

Supports both single-cell mode (per-cell TRS) and cell-type mode (per-cell-type TRS).
"""

from __future__ import annotations
import argparse
import sys
import time
import re

import numpy as np
import pandas as pd

from scrbp.utils.stats import regulon_specificity_scores


# -------------------- small utils --------------------
def robust_minmax_vec(x: np.ndarray, q_hi: float = 0.99) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return x
    hi = np.nanquantile(x, q_hi)
    lo = np.nanmin(x)
    if not np.isfinite(hi) or not np.isfinite(lo):
        return np.zeros_like(x)
    if hi <= lo:
        hi = lo + 1e-9
    x = np.clip(x, None, hi)
    return (x - lo) / (hi - lo)


def empirical_p(real, nulls, include_equals: bool = True):
    """
    Empirical right-tailed p-value.
    real  : shape (...,) or scalar
    nulls : shape (..., K)
    include_equals: if True, use '>='; if False, use '>' to avoid tie inflation.
    """
    real  = np.asarray(real, dtype=float)
    nulls = np.asarray(nulls, dtype=float)
    if include_equals:
        ge = (nulls >= real[..., None]).sum(axis=-1)
    else:
        ge = (nulls >  real[..., None]).sum(axis=-1)
    return (1.0 + ge) / (nulls.shape[-1] + 1.0)


def robust_minmax(z, z_null, q_cap=0.99, range_eps=1e-9):
    """
    Per-position robust min-max using its own NULL distribution.

    z      : shape (...,) (e.g., CT dimension)
    z_null : shape (..., K)
    """
    z      = np.asarray(z, dtype=float)
    z_null = np.asarray(z_null, dtype=float)

    lo = np.nanmin(z_null, axis=-1)
    hi = np.nanquantile(z_null, q_cap, axis=-1)

    rng = hi - lo
    rng = np.where(rng > range_eps, rng, range_eps)

    z_cap      = np.clip(z,      lo,            hi)
    z_null_cap = np.clip(z_null, lo[..., None], hi[..., None])

    z_mm      = (z_cap      - lo) / rng
    z_null_mm = (z_null_cap - lo[..., None]) / rng[..., None]
    return z_mm, z_null_mm


def sd2(a, b):
    """Sample std of two numbers along-element (ddof=1): sd(a,b) = |a-b| / sqrt(2)."""
    return np.abs(a - b) / np.sqrt(2.0)


# -------------------- RAS I/O (CSV & LOOM) --------------------
def _load_ras_csv_regulons_by_cells(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    if df.shape[1] < 2:
        raise SystemExit(f"[RAS CSV] Expect first column=Regulon + many cell columns. Got shape={df.shape}")
    set_col = df.columns[0]
    mat = (df.set_index(set_col)
             .apply(pd.to_numeric, errors="coerce")
             .fillna(0.0))                 # (regulons x cells)
    ras = mat.T                            # -> (cells x regulons)
    ras = ras.loc[~ras.index.duplicated(keep="first")]
    ras = ras.loc[:, ~ras.columns.duplicated(keep="first")]
    if ras.shape[1] == 0:
        raise SystemExit("[RAS CSV] No regulon columns after cleaning.")
    return ras


def _load_ras_loom(path: str, dtype="float32") -> pd.DataFrame:
    try:
        import loompy
    except ImportError as e:
        raise SystemExit("Please install loompy to read .loom files: pip install loompy") from e
    with loompy.connect(path, mode="r") as ds:
        regs  = np.array(ds.ra["Regulon"], dtype=object) if "Regulon" in ds.ra else np.arange(ds.shape[0], dtype=object)
        cells = np.array(ds.ca["CellID"], dtype=object)  if "CellID"  in ds.ca else np.arange(ds.shape[1], dtype=object)
        mat = ds[:, :].astype(dtype, copy=False)         # (regulons x cells)
    ras = pd.DataFrame(mat.T, index=cells, columns=regs) # -> (cells x regulons)
    ras = ras.loc[~ras.index.duplicated(keep="first")]
    ras = ras.loc[:, ~ras.columns.duplicated(keep="first")]
    if ras.shape[1] == 0:
        raise SystemExit("[RAS LOOM] No regulon columns after cleaning.")
    return ras


def load_ras_any(path: str) -> pd.DataFrame:
    lower = path.lower()
    if lower.endswith(".loom"):
        return _load_ras_loom(path)
    if lower.endswith(".csv") or lower.endswith(".csv.gz"):
        return _load_ras_csv_regulons_by_cells(path)
    raise SystemExit(f"[RAS] Unsupported extension for {path}. Use .csv/.csv.gz or .loom")


# -------------------- RGS I/O --------------------
_MLOG10P_CANDS = ["RGS_mlog10P", "mlog10P", "neglog10P", "-log10P", "minus_log10P"]
_Z_CANDS       = ["z_RGS", "RGS_z", "z", "Z", "z_score"]


def _choose_rgs_score_series(df: pd.DataFrame, score_mode: str) -> np.ndarray:
    score_mode = score_mode.lower()
    if score_mode == "z":
        z_col = next((c for c in _Z_CANDS if c in df.columns), None)
        if z_col is None:
            raise SystemExit("[RGS] z-score requested but no z column found.")
        return df[z_col].astype(float).to_numpy()
    elif score_mode == "mlog10p":
        m_col = next((c for c in _MLOG10P_CANDS if c in df.columns), None)
        if m_col is not None:
            return df[m_col].astype(float).to_numpy()
        if "P" not in df.columns:
            raise SystemExit("[RGS] Need 'P' to compute -log10(P).")
        P = df["P"].astype(float).to_numpy()
        return -np.log10(np.clip(P, 1e-300, None))
    else:
        raise SystemExit(f"[RGS] Unknown --rgs-score: {score_mode}")


def read_rgs_sc(path: str, score_mode: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    set_col = df.columns[0]
    if "P" not in df.columns:
        raise SystemExit("[RGS] Need column 'P'.")
    score = _choose_rgs_score_series(df, score_mode)
    out = (df[[set_col, "P"]]
           .assign(SCORE=score)
           .rename(columns={set_col: "SET"})
           .dropna(subset=["SET"])
           .set_index("SET"))
    if not out.index.is_unique:
        out = (out.sort_values("P")
                 .groupby(level=0, sort=False)
                 .agg({"P": "min", "SCORE": "mean"}))
    return out  # index=RBP


def parse_rgs_ct(path: str, score_mode: str):
    df = pd.read_csv(path)
    reg_col = df.columns[0]
    if "P" not in df.columns:
        raise SystemExit("[RGS-CT] Need column 'P'.")
    score = _choose_rgs_score_series(df, score_mode)
    df2 = df[[reg_col, "P"]].copy()
    df2["SCORE"] = score

    regs = df2[reg_col].astype(str)
    is_real = regs.str.endswith("__REAL")
    is_null = regs.str.contains("__NULL_")

    s_real, p_real, s_nulls = {}, {}, {}

    # REAL
    for _, r in df2[is_real].iterrows():
        rbp = str(r[reg_col]).split("__")[0]
        s_real[rbp] = float(r["SCORE"])
        p_real[rbp] = float(r["P"])

    # NULL (stable sort)
    base = regs.str.split("__").str[0]
    for rbp, sub in df2[is_null].groupby(base[is_null]):
        tmp = sub.copy()
        tmp["_id"] = (tmp[reg_col].str.extract(r"__NULL_(\d+)$", expand=False).astype("Int64"))
        tmp = tmp.sort_values(by=["_id", reg_col], kind="stable")
        s_nulls[rbp] = tmp["SCORE"].to_numpy(dtype=float)

    RBPs = [rbp for rbp in s_real if rbp in s_nulls and len(s_nulls[rbp]) > 0]
    return RBPs, s_real, s_nulls, p_real


# rss_matrix moved to scrbp.utils.stats as regulon_specificity_scores
# Use: regulon_specificity_scores(ras_cells_x_rbp, labels)


# ---- z-scores ----
def zscore_from_null(real, nulls, eps: float = 1e-9):
    real  = np.asarray(real, dtype=float)
    nulls = np.asarray(nulls, dtype=float)
    mu  = nulls.mean(axis=-1)
    sd  = nulls.std(axis=-1, ddof=1)
    sd  = np.where(sd < eps, np.nan, sd)
    return (real - mu) / sd


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float).ravel()
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = np.empty_like(ranked)
    prev = 1.0
    for i in range(n-1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        if val > prev: val = prev
        prev = val
        q[i] = min(val, 1.0)
    q_full = np.empty_like(p)
    q_full[order] = q
    return q_full.reshape(pvals.shape)


def bh_fdr_by_row(P: np.ndarray) -> np.ndarray:
    P = np.asarray(P, float)
    Q = np.empty_like(P)
    for i in range(P.shape[0]):
        Q[i, :] = bh_fdr(P[i, :])
    return Q


# ---- single-cell z-score helpers ----
def zscore_cols(df: pd.DataFrame, eps: float = 1e-9) -> pd.DataFrame:
    X = df.to_numpy(dtype=float)
    mu = np.nanmean(X, axis=0, keepdims=True)
    sd = np.nanstd(X, ddof=1, axis=0, keepdims=True)
    sd = np.where((~np.isfinite(sd)) | (sd < eps), np.nan, sd)
    Z = (X - mu) / sd
    return pd.DataFrame(Z, index=df.index, columns=df.columns)


def zscore_vec(x: np.ndarray, eps: float = 1e-9) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    mu = np.nanmean(x)
    sd = np.nanstd(x, ddof=1)
    if (not np.isfinite(sd)) or (sd < eps):
        return np.full_like(x, np.nan, dtype=float)
    return (x - mu) / sd


# ---- progress ----
def _fmt_td(sec: float) -> str:
    s = int(sec); m, s = divmod(s, 60); h, m = divmod(m, 60)
    return f"{h:02d}:{m:02d}:{s:02d}"


def _progress(prefix: str, done: int, total: int, t0: float, bar_len: int = 28):
    done = max(0, min(done, total))
    pct = 0.0 if total == 0 else done / total
    filled = int(bar_len * pct)
    bar = "#" * filled + "-" * (bar_len - filled)
    elapsed = time.time() - t0
    sys.stdout.write(f"\r{prefix} [{bar}] {pct:6.2%}  elapsed {_fmt_td(elapsed)}")
    sys.stdout.flush()


# ---------- helpers: parse names & RAS(RSS) with REAL+NULL columns ----------
def _split_kind(name: str):
    """
    'AGO1__REAL' -> ('AGO1','REAL',None)
    'AGO1__NULL_0007' -> ('AGO1','NULL',7)
    fallback: ('name','',None)
    """
    s = str(name)
    m = re.match(r"^(.*)__REAL$", s, flags=re.I)
    if m: return m.group(1), "REAL", None
    m = re.match(r"^(.*)__NULL_(\d+)$", s, flags=re.I)
    if m: return m.group(1), "NULL", int(m.group(2))
    return s, "", None


def _rss_from_ras_with_nulls_per_rbp(ras_cells_by_reg: pd.DataFrame,
                                     labels: pd.Series):
    """
    Return per-RBP variable-length NULLs.
    Output:
      RBPs: list[str]
      ras_real_mat: (CT x R) RSS matrix for REAL columns
      ras_null_list: list length R; ras_null_list[j] has shape (CT x K_ras_j) [RAW]
      ct_levels: index of CTs
      K_ras_by: dict RBP -> K_ras_j
    """
    # align labels to cells
    if not ras_cells_by_reg.index.equals(labels.index):
        labels = labels.reindex(ras_cells_by_reg.index)
    mask = ~labels.isna()
    ras_cells_by_reg = ras_cells_by_reg.loc[mask]
    labels = labels.loc[mask]

    # parse column names
    parsed = [_split_kind(c) for c in ras_cells_by_reg.columns]
    rbp_of_col    = [p[0] for p in parsed]
    kind_of_col   = [p[1] for p in parsed]
    nullid_of_col = [p[2] for p in parsed]
    cols = ras_cells_by_reg.columns

    real_col = {}
    null_cols = {}
    for c, rbp, kind, nid in zip(cols, rbp_of_col, kind_of_col, nullid_of_col):
        if kind == "REAL":
            real_col[rbp] = c
        elif kind == "NULL" and nid is not None:
            null_cols.setdefault(rbp, []).append((nid, c))

    RBPs = sorted([r for r in real_col if r in null_cols and len(null_cols[r]) > 0])
    if not RBPs:
        raise SystemExit("[CT/GENESET-NULL] No RBP has both REAL and NULL columns in RAS.")

    # compute RSS for all columns once (RAW)
    rss_all = regulon_specificity_scores(ras_cells_by_reg, labels)  # (CT x all-reg cols)
    ct_levels = rss_all.index

    R = len(RBPs)
    ras_real_mat = np.empty((len(ct_levels), R), dtype=np.float32)
    ras_null_list = []
    K_ras_by = {}

    for j, rbp in enumerate(RBPs):
        # REAL
        ras_real_mat[:, j] = rss_all[real_col[rbp]].to_numpy(dtype=np.float32)
        # NULL: keep all ids (sorted)
        ids_cols = sorted(null_cols[rbp], key=lambda x: x[0])
        K_j = len(ids_cols)
        K_ras_by[rbp] = K_j
        if K_j == 0:
            ras_null_list.append(np.empty((len(ct_levels), 0), dtype=np.float32))
        else:
            stack = np.stack([rss_all[c].to_numpy(dtype=np.float32) for _, c in ids_cols], axis=-1)
            ras_null_list.append(stack)  # (CT x K_j) RAW

    return RBPs, ras_real_mat, ras_null_list, ct_levels, K_ras_by


# -------------------- Modes --------------------
def run_mode_sc(args):
    _sc_t0 = time.time()
    ras = load_ras_any(args.ras)     # CSV or LOOM via --ras
    rgs = read_rgs_sc(args.rgs_csv, score_mode=args.rgs_score)

    # Align regulons
    common = [r for r in ras.columns if r in rgs.index]
    if len(common) == 0:
        raise SystemExit("[SC] No overlapping regulons between RAS and RGS.")
    dropped = len(ras.columns) - len(common)
    if dropped > 0:
        print(f"[SC][WARN] {dropped} RAS regulons not found in RGS; keeping {len(common)} overlap.")
    ras = ras[common]

    # RAS_mm across cells
    ras_mm = ras.apply(lambda col: robust_minmax_vec(col.values, q_hi=float(args.q_hi_ras)),
                       axis=0, result_type="broadcast")
    ras_mm = pd.DataFrame(ras_mm.to_numpy(), index=ras.index, columns=ras.columns)

    # RGS_mm across regulons
    score_vec = rgs.loc[ras.columns, "SCORE"].astype(float).values
    rgs_mm_vec = robust_minmax_vec(score_vec, q_hi=float(args.q_hi_rgs))
    rgs_mm = pd.Series(rgs_mm_vec, index=ras.columns, name="RGS_mm")

    # TRS = RAS + RGS - lambda * sd(RAS,RGS)
    lam = float(args.lambda_penalty)
    trs = ras_mm.to_numpy() + rgs_mm.values[None, :] - lam * sd2(ras_mm.to_numpy(), rgs_mm.values[None, :])
    trs_df = pd.DataFrame(trs, index=ras.index, columns=ras.columns)

    # z-scores
    z_ras_df  = zscore_cols(ras_mm)
    z_trs_df  = zscore_cols(trs_df)
    z_rgs_vec = zscore_vec(rgs_mm.values)

    # Outputs
    trs_df.to_csv(f"{args.out_prefix}.sc.TRS_matrix.csv")
    pd.DataFrame({
        "cell":   np.repeat(trs_df.index.values, trs_df.shape[1]),
        "RBP":    np.tile(trs_df.columns.values, trs_df.shape[0]),
        "RAS_mm": ras_mm.to_numpy().ravel(),
        "RGS_mm": np.tile(rgs_mm.values, trs_df.shape[0]),
        "TRS_mm": trs_df.to_numpy().ravel(),
        "z_RAS":  z_ras_df.to_numpy().ravel(),
        "z_RGS":  np.tile(z_rgs_vec, trs_df.shape[0]),
        "z_TRS":  z_trs_df.to_numpy().ravel(),
    }).to_csv(f"{args.out_prefix}.sc.TRS_long.csv", index=False)
    print("[SC] Wrote:",
          f"{args.out_prefix}.sc.TRS_matrix.csv",
          f"{args.out_prefix}.sc.TRS_long.csv")
    print(f"[SC] Elapsed {_fmt_td(time.time() - _sc_t0)}")


def run_mode_ct(args):
    """
    CT mode (gene-set NULL, per-RBP K):
      - RAS input must carry columns RBP__REAL and RBP__NULL_#### (AUCell/RSS-ready).
      - RGS CSV must carry rows RBP__REAL and RBP__NULL_#### (GSA output).
      - For each RBP j, K_j = min(K_ras_j, K_rgs_j).
      - p_RAS/p_RGS on RAW vs RAW NULL; p_TRS on min-max TRS vs min-max TRS NULL (strict '>').
    """
    _ct_t0 = time.time()

    # 1) RAS (cells x regulons, with REAL+NULL columns)
    ras_raw = load_ras_any(args.ras)

    # 2) Cell-type labels & min cells per CT
    ct = pd.read_csv(args.celltypes_csv)
    if ct.shape[1] < 2:
        raise SystemExit("[CT] --celltypes-csv must have two columns: cell_id,cell_type")
    ct.columns = ["cell_id", "cell_type"] + list(ct.columns[2:])
    labels = pd.Series(ct["cell_type"].values, index=ct["cell_id"].values, name="cell_type")
    if not ras_raw.index.equals(labels.index):
        labels = labels.reindex(ras_raw.index)
    mask = ~labels.isna()
    ras_raw = ras_raw.loc[mask]
    labels  = labels.loc[mask]
    vc = labels.value_counts()
    min_ct = int(getattr(args, "min_cells_pert_ct", 25))
    keep_cts = vc[vc >= min_ct].index
    keep_mask = labels.isin(keep_cts)
    ras_raw = ras_raw.loc[keep_mask]
    labels  = labels.loc[keep_mask]

    # 3) From RAS with REAL+NULL columns -> per-RBP ras_real & ras_null_list (RAW)
    RBPs_ras, ras_real_mat, ras_null_list, ct_levels, K_ras_by = _rss_from_ras_with_nulls_per_rbp(ras_raw, labels)

    # 4) Parse RGS (REAL+NULL rows)
    RBPs_rgs, s_real_map, s_null_map, _p_real_map = parse_rgs_ct(args.rgs_csv, score_mode=args.rgs_score)

    # 5) Overlap RBPs with available NULLs on both sides
    RBPs = [r for r in RBPs_ras if (r in RBPs_rgs and K_ras_by.get(r, 0) > 0 and len(s_null_map.get(r, [])) > 0)]
    if not RBPs:
        raise SystemExit("[CT/GENESET-NULL] No overlap RBPs with available NULLs on both RAS and RGS sides.")

    pos_ras = {rbp: i for i, rbp in enumerate(RBPs_ras)}

    # 6) Prepare outputs (CT x R)
    R = len(RBPs)
    C = len(ct_levels)
    ras_mm_mat = np.full((C, R), np.nan, dtype=float)
    rgs_mm_mat = np.full((C, R), np.nan, dtype=float)
    trs_mm_mat = np.full((C, R), np.nan, dtype=float)

    p_ras_mat  = np.full((C, R), np.nan, dtype=float)
    p_trs_mat  = np.full((C, R), np.nan, dtype=float)
    z_ras_mat  = np.full((C, R), np.nan, dtype=float)
    z_trs_mat  = np.full((C, R), np.nan, dtype=float)

    p_rgs_vec  = np.full((R,), np.nan, dtype=float)
    z_rgs_vec  = np.full((R,), np.nan, dtype=float)
    K_vec      = np.full((R,), 0,       dtype=int)

    lam = float(args.lambda_penalty)
    q_cap_ras = float(getattr(args, "q_cap_ras", 0.99))
    q_cap_rgs = float(getattr(args, "q_cap_rgs", 0.99))

    # 7) Per-RBP computation
    for j, rbp in enumerate(RBPs):
        j_ras = pos_ras[rbp]
        ras_real_j = ras_real_mat[:, j_ras]          # (CT,) RAW RSS (REAL)
        ras_null_j = ras_null_list[j_ras]            # (CT x K_ras_j) RAW RSS (NULL)
        K_ras_j = ras_null_j.shape[1]

        s_real_j = float(s_real_map[rbp])            # scalar RAW RGS (REAL)
        s_null_j = np.asarray(s_null_map[rbp], dtype=float)  # (K_rgs_j,) RAW RGS (NULL)
        K_rgs_j = s_null_j.shape[0]

        K_j = int(min(K_ras_j, K_rgs_j))
        if K_j <= 0:
            continue
        K_vec[j] = K_j

        # truncate to common K_j (RAW)
        ras_null_j = ras_null_j[:, :K_j]             # (CT x K_j)
        s_null_j   = s_null_j[:K_j]                  # (K_j,)

        # ---------- TRUE empirical p on RAW ----------
        p_ras_mat[:, j] = empirical_p(ras_real_j, ras_null_j, include_equals=True)   # RAW RSS
        p_rgs_vec[j]    = empirical_p(s_real_j,  s_null_j,   include_equals=True)    # RAW RGS

        # ---------- min-max for TRS construction & TRS p ----------
        ras_mm_j,  ras_null_mm_j  = robust_minmax(ras_real_j, ras_null_j, q_cap=q_cap_ras)  # (CT,), (CT x K_j)
        rgs_mm_j,  rgs_null_mm_j  = robust_minmax(s_real_j,  s_null_j,   q_cap=q_cap_rgs)   # ( ),  (K_j,)

        rgs_mm_ct_j   = np.broadcast_to(rgs_mm_j, (C,))                    # (CT,)
        rgs_null_ct_j = np.broadcast_to(rgs_null_mm_j[None, :], (C, K_j))  # (CT x K_j)

        trs_real_j = ras_mm_j + rgs_mm_ct_j - lam * sd2(ras_mm_j, rgs_mm_ct_j)                    # (CT,)
        trs_null_j = ras_null_mm_j + rgs_null_ct_j - lam * sd2(ras_null_mm_j, rgs_null_ct_j)      # (CT x K_j)

        # p_TRS: strict right tail to avoid cap ties
        p_trs_mat[:, j] = empirical_p(trs_real_j, trs_null_j, include_equals=False)

        # z-scores
        z_ras_mat[:, j] = zscore_from_null(ras_real_j, ras_null_j)   # RAW
        z_rgs_vec[j]    = zscore_from_null(s_real_j,  s_null_j)      # RAW
        z_trs_mat[:, j] = zscore_from_null(trs_real_j, trs_null_j)   # min-max TRS

        # store mm values for output
        ras_mm_mat[:, j] = ras_mm_j
        rgs_mm_mat[:, j] = rgs_mm_ct_j
        trs_mm_mat[:, j] = trs_real_j

    # 8) FDR (row-wise BH)
    do_fdr = bool(getattr(args, "do_fdr", True))
    if do_fdr:
        q_TRS = bh_fdr_by_row(p_trs_mat)
        q_RAS = bh_fdr_by_row(p_ras_mat)
        q_RGS_v = bh_fdr(p_rgs_vec)
        q_RGS = np.broadcast_to(q_RGS_v[None, :], (C, R))
    else:
        q_TRS = q_RAS = q_RGS = None

    # 9) Output
    ct_index = list(ct_levels)
    rows = []
    for i, ct_name in enumerate(ct_index):
        for j, rbp in enumerate(RBPs):
            row = {
                "cell_type": ct_name,
                "RBP": rbp,
                "RAS_mm": float(ras_mm_mat[i, j]) if np.isfinite(ras_mm_mat[i, j]) else np.nan,
                "RGS_mm": float(rgs_mm_mat[i, j]) if np.isfinite(rgs_mm_mat[i, j]) else np.nan,
                "TRS_mm": float(trs_mm_mat[i, j]) if np.isfinite(trs_mm_mat[i, j]) else np.nan,
                "p_RAS":  float(p_ras_mat[i, j])  if np.isfinite(p_ras_mat[i, j])  else np.nan,   # RAW-based p
                "p_RGS":  float(p_rgs_vec[j])     if np.isfinite(p_rgs_vec[j])     else np.nan,   # RAW-based p
                "p_TRS":  float(p_trs_mat[i, j])  if np.isfinite(p_trs_mat[i, j])  else np.nan,   # strict '>'
                "z_RAS":  float(z_ras_mat[i, j])  if np.isfinite(z_ras_mat[i, j])  else np.nan,   # RAW-based z
                "z_RGS":  float(z_rgs_vec[j])     if np.isfinite(z_rgs_vec[j])     else np.nan,   # RAW-based z
                "z_TRS":  float(z_trs_mat[i, j])  if np.isfinite(z_trs_mat[i, j])  else np.nan,
                "RGS_score_type": args.rgs_score,
                "K_null": int(K_vec[j]),
            }
            if do_fdr:
                row.update({
                    "q_RAS": float(q_RAS[i, j]) if np.isfinite(q_RAS[i, j]) else np.nan,
                    "q_RGS": float(q_RGS[i, j]) if np.isfinite(q_RGS[i, j]) else np.nan,
                    "q_TRS": float(q_TRS[i, j]) if np.isfinite(q_TRS[i, j]) else np.nan,
                })
            rows.append(row)

    out_long = f"{args.out_prefix}.ct.TRS_long.csv"
    out_mat  = f"{args.out_prefix}.ct.TRS_matrix.csv"
    pd.DataFrame(rows).to_csv(out_long, index=False)
    pd.DataFrame(trs_mm_mat, index=ct_index, columns=RBPs).to_csv(out_mat)

    # per-RBP K_null for audit
    pd.DataFrame({"RBP": RBPs, "K_null": K_vec}).to_csv(f"{args.out_prefix}.ct.Knull_perRBP.csv", index=False)

    print("[CT/GENESET-NULL, per-RBP K | TRUE empirical p] Wrote:",
          out_long, out_mat, f"{args.out_prefix}.ct.Knull_perRBP.csv")
    print(f"[CT/GENESET-NULL, per-RBP K | TRUE empirical p] Elapsed {_fmt_td(time.time() - _ct_t0)}")


def register_subcommand(subparsers):
    """Register this subcommand with the main parser."""
    p = subparsers.add_parser(
        'trs',
        help='Compute Trait Relevance Scores (RAS + RGS)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )
    p.add_argument("--mode", choices=["sc", "ct"], required=True)

    # RAS can be CSV or LOOM
    p.add_argument("--ras", required=True,
                   help="RAS file (.csv/.csv.gz with first column=Regulon; or .loom with rows=regulons, cols=cells).")

    # RGS inputs & options
    p.add_argument("--rgs-csv", required=True,
                   help="[SC] RGS with RBP rows; [CT] RGS with RBP__REAL and RBP__NULL_#### rows (e.g., GSA).")
    p.add_argument("--rgs-score", choices=["mlog10p", "z"], default="mlog10p",
                   help="RGS score: 'mlog10p' uses -log10(P) (compute from P if missing); 'z' uses z_RGS / RGS_z.")

    p.add_argument("--out-prefix", required=True)
    p.add_argument("--lambda-penalty", type=float, default=1.0)

    # SC-only robust min-max caps
    p.add_argument("--q-hi-ras", dest="q_hi_ras", type=float, default=0.99,
                   help="[SC] Upper quantile for robust min-max of RAS across cells.")
    p.add_argument("--q-hi-rgs", dest="q_hi_rgs", type=float, default=0.99,
                   help="[SC] Upper quantile for robust min-max of RGS across regulons.")

    # CT controls
    p.add_argument("--do-fdr", type=int, default=1, help="[CT] Apply BH-FDR (0/1).")
    p.add_argument("--celltypes-csv", help="[CT] CSV with two columns: cell_id,cell_type")
    p.add_argument("--min_cells_pert_ct", type=int, default=25,
                   help="[CT] Minimum cells per cell-type; CTs with fewer cells are dropped.")
    p.add_argument("--q-cap-ras", dest="q_cap_ras", type=float, default=0.99,
                   help="[CT] Upper quantile cap for RAS NULL-based min-max.")
    p.add_argument("--q-cap-rgs", dest="q_cap_rgs", type=float, default=0.99,
                   help="[CT] Upper quantile cap for RGS NULL-based min-max.")
    p.add_argument("--save-nulls-prefix", default=None,
                   help="[CT] Ragged K makes 3D saving tricky; omitted in v14.")

    p.set_defaults(func=main)
    return p


def main(args):
    """Main execution function."""
    if args.mode == "ct" and not args.celltypes_csv:
        sys.exit("[CT] Please provide --celltypes-csv (two columns: cell_id,cell_type).")
    if args.mode == "sc":
        run_mode_sc(args)
    else:
        run_mode_ct(args)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    subp = p.add_subparsers()
    register_subcommand(subp)
    args = p.parse_args()
    main(args)
