#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scRBP rgs — Regulon-level Genetic Association Score (RGS) via MAGMA gene-set analysis

Integrates genetic association signals from GWAS into each RBP regulon by
running MAGMA competitive gene-set regression.  For each regulon, MAGMA
evaluates whether its member genes carry stronger polygenic signals than
background genes, after adjusting for gene size, SNP density, and effective
gene size (and their log-transformed terms).

The resulting Regulon-level Genetic Association Score (RGS) is a standardised
estimate of the relative enrichment of GWAS signal within a focal regulon
compared with genes assigned to other regulons.  A large positive RGS
indicates that the RBP regulon is significantly enriched for disease-associated
genetic variants.

In cell-type mode (``--mode ct``), size-matched null regulons are automatically
generated and used to compute empirical p-values and z-scores, using a 4D null
distribution matched on gene count, SNP count (NSNPS), number of principal
components (NPARAM), mean expression, and detection rate (pct_detected).

Outputs
-------
- Real regulon MAGMA GSA results: ``<out>_real.csv``
- Null regulon MAGMA GSA results: ``<out>_null.csv`` (cell-type mode only)
- Empirical p-values and z-scores: ``<out>_empirical.csv``

Prerequisites: MAGMA binary (≥v1.10) and a precomputed ``.genes.raw`` file
from MAGMA gene analysis using GWAS summary statistics.
"""

import os
import shlex
import subprocess
import argparse
import io
import re
import glob
import tempfile
import gzip
import shutil

import numpy as np
import pandas as pd
from scipy.stats import norm

# Optional (for fast CSV): polars
try:
    import polars as pl
except Exception:
    pl = None


from scrbp.utils.io import (
    read_gmt,
    write_gmt,
    load_expression_matrix,
    compute_expr_stats_from_df,
)

# ------------------------- Common utils -------------------------

def str2bool(v):
    if isinstance(v, bool): return v
    v = str(v).strip().lower()
    if v in ("y", "yes", "t", "true", "1"):  return True
    if v in ("n", "no", "f", "false", "0"):  return False
    raise argparse.ArgumentTypeError("Boolean value expected (True/False).")


def read_gsa_out(gsa_out_path: str, one_sided: bool = True) -> pd.DataFrame:
    with open(gsa_out_path, "r") as f:
        raw = [ln.rstrip("\n") for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    try:
        start = next(i for i, s in enumerate(raw) if s.strip().upper().startswith("VARIABLE"))
    except StopIteration:
        raise ValueError(f"Cannot find header line 'VARIABLE' in {gsa_out_path}")
    txt = "\n".join(raw[start:])
    df = pd.read_fwf(io.StringIO(txt), header=0)
    df = df[df["VARIABLE"].astype(str).str.upper() != "VARIABLE"].copy()

    for col in ["NGENES", "BETA", "BETA_STD", "SE", "P"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    if "P" not in df.columns:
        raise ValueError(f"'P' not found after parsing {gsa_out_path}. Got columns: {list(df.columns)}")

    p = df["P"].clip(lower=np.finfo(float).tiny)
    if one_sided:
        df["RGS_z"] = norm.isf(p)
    else:
        sign = np.sign(df.get("BETA", 1.0).fillna(1.0))
        df["RGS_z"] = norm.isf(p/2.0) * sign
    df["RGS_mlog10P"] = -np.log10(p)

    if "VARIABLE" in df.columns: df = df.rename(columns={"VARIABLE": "Regulon"})
    if "TYPE" in df.columns:     df = df.rename(columns={"TYPE": "GeneSet"})
    if "GeneSet" in df.columns:
        df["GeneSet"] = df["GeneSet"].astype(str).apply(
            lambda x: "RBP_regulon" if x.strip().upper() == "SET" else x
        )
    return df


def read_genes_raw_head(path):
    covar = []
    with open(path) as f:
        for ln in f:
            if ln.startswith("# COVAR"):
                covar = ln.split("=", 1)[1].split()
                break
    need = 7 + len(covar)
    rows = []
    with open(path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"): continue
            toks = ln.strip().split()
            rows.append(toks[:need])
    cols = ["GENE", "CHR", "START", "STOP", "NSNPS", "NPARAM", "N"] + covar
    df = pd.DataFrame(rows, columns=cols)
    for c in ["GENE", "NSNPS", "NPARAM"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["GENE"] = df["GENE"].astype(int)
    return df[["GENE", "NSNPS", "NPARAM"]]


def load_entrez_symbol_from_gene_loc(path):
    df = pd.read_csv(path, sep=r"\s+", header=None, engine="python", usecols=[0, 5])
    df.columns = ["entrez", "symbol"]
    df = df.dropna()
    df = df[(df["symbol"] != "-") & (df["symbol"] != "")]
    df["entrez"] = df["entrez"].astype(int)
    e2s = dict(zip(df["entrez"], df["symbol"].astype(str)))
    s2e = {}
    for e, s in e2s.items():
        s2e.setdefault(s.upper(), e)
    return e2s, s2e


# read_gmt, write_gmt moved to scrbp.utils.io

def qbin(x, q=10):
    qq = pd.qcut(x, q=q, labels=False, duplicates="drop")
    return qq.fillna(0).astype(int)


# load_expression_matrix, compute_expr_stats_from_df moved to scrbp.utils.io

def read_expr_stats(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")  # auto-sep
    need = {"symbol", "mean_expr", "pct_detected"}
    # normalize column names
    lower_map = {c: c.lower() for c in df.columns}
    df = df.rename(columns=lower_map)
    if not need.issubset(df.columns):
        raise ValueError(f"expr-stats file must contain columns: symbol, mean_expr, pct_detected; got {list(df.columns)}")
    df["symbol"] = df["symbol"].astype(str).str.upper()
    return df[["symbol", "mean_expr", "pct_detected"]].copy()


# ------------------------- NULL construction (4D bins) -------------------------

def make_null_for_regulon_4d(entrez_list,
                             bg_df,                # columns: GENE, NSNPS, NPARAM, mean_expr, pct_detected
                             n_null=1000, q=10, rng=None, exclude_self=True):
    if rng is None:
        rng = np.random.default_rng(1)

    reg = bg_df[bg_df["GENE"].isin(entrez_list)].copy()

    # Build bins: NSNPS, NPARAM, mean_expr, pct_detected
    for k, col in [("bin_nsnps", "NSNPS"), ("bin_nparam", "NPARAM"),
                   ("bin_mean", "mean_expr"), ("bin_pct", "pct_detected")]:
        bg_df[k] = qbin(bg_df[col], q=q)
        reg[k]   = qbin(reg[col],   q=q)

    # target composition
    need4 = reg.groupby(["bin_nsnps", "bin_nparam", "bin_mean", "bin_pct"]).size().to_dict()
    need2 = reg.groupby(["bin_nsnps", "bin_nparam"]).size().to_dict()  # 2D fallback

    pool = bg_df[~bg_df["GENE"].isin(entrez_list)] if exclude_self else bg_df

    # Precompute buckets
    buckets4 = {k: pool[(pool["bin_nsnps"]==k[0]) & (pool["bin_nparam"]==k[1]) &
                        (pool["bin_mean"]==k[2]) & (pool["bin_pct"]==k[3])]["GENE"].values
                for k in need4.keys()}
    buckets2 = {k: pool[(pool["bin_nsnps"]==k[0]) & (pool["bin_nparam"]==k[1])]["GENE"].values
                for k in need2.keys()}

    null_sets = []
    for _ in range(n_null):
        picked = []
        # try to match 4D composition first
        for k4, cnt in need4.items():
            cand = buckets4.get(k4, np.array([], dtype=int))
            if len(cand) >= cnt and cnt > 0:
                pick = rng.choice(cand, size=cnt, replace=False)
                picked.extend(pick.tolist())
            else:
                # fallback: use 2D bucket (nsnps x nparam)
                k2 = (k4[0], k4[1])
                cand2 = buckets2.get(k2, np.array([], dtype=int))
                if len(cand2) >= cnt and cnt > 0:
                    pick = rng.choice(cand2, size=cnt, replace=False)
                    picked.extend(pick.tolist())
                else:
                    # last fallback: global pool with replacement
                    cand_all = pool["GENE"].values
                    pick = rng.choice(cand_all, size=cnt, replace=True)
                    picked.extend(pick.tolist())
        null_sets.append(picked)
    return null_sets


def split_regulon_name(name: str):
    m_real = re.match(r"^(.*)__REAL$", name)
    if m_real:
        return m_real.group(1), "REAL", 0
    m_null = re.match(r"^(.*)__NULL_(\d+)$", name)
    if m_null:
        return m_null.group(1), "NULL", int(m_null.group(2))
    return name, "", None


def register_subcommand(subparsers):
    """Register this subcommand with the main parser."""
    ap = subparsers.add_parser(
        'rgs',
        help='Regulon Gene Set analysis with MAGMA',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )

    # Core GSA
    ap.add_argument("--mode", choices=["sc", "ct"], required=True,
                    help="sc: real regulons only; ct: add matched null regulons per real set.")
    ap.add_argument("--magma", required=True, help="Path to MAGMA binary")
    ap.add_argument("--genes-raw", required=True, help="MAGMA <prefix>.genes.raw from gene analysis")
    ap.add_argument("--sets", required=True, help="Input regulon GMT (symbol or entrez; see --id-type)")
    ap.add_argument("--id-type", choices=["entrez", "symbol"], default="entrez",
                    help="Gene ID type in --sets")
    ap.add_argument("--out", required=True, help="Output prefix")
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--cleanup-out", type=str2bool, nargs="?", const=True, default=True)
    ap.add_argument("--extra", default="", help='Extra MAGMA args, e.g. "--model condition-hide"')
    ap.add_argument("--min_genes", "--min-genes", dest="min_genes", type=int, default=0)

    # ct-specific
    ap.add_argument("--gene-loc", help="MAGMA NCBI*.gene.loc for Entrez<->Symbol (needed if we must map expr_stats symbols)")
    ap.add_argument("--n-null", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=2025)
    ap.add_argument("--q-bins", type=int, default=10)
    ap.add_argument("--exclude-self", type=str2bool, nargs="?", const=True, default=True)

    # Expr-stats (for ct NULL)
    ap.add_argument("--expr-stats", help="Precomputed expr stats TSV (symbol, mean_expr, pct_detected).")
    ap.add_argument("--emit-expr-stats", type=str2bool, nargs="?", const=True, default=False,
                    help="If True and --expr-stats not provided, compute from --matrix-stats and save.")
    ap.add_argument("--matrix-stats",
                    help="Expression matrix to compute expr-stats when --emit-expr-stats is True.")
    ap.add_argument("--expr-stats-out",
                    help="Where to save computed expr-stats (default: <out>_expr_stats.tsv).")
    ap.add_argument("--dtype", default="float32", help="dtype for --matrix-stats loading")

    # Persist REAL+NULL GMT for scRBP_trs
    ap.add_argument("--save-null-gmt",
                    help="Save REAL+NULL Entrez GMT for reuse (default: <out>_REAL_PLUS_NULLS.entrez.gmt).")
    ap.add_argument("--save-null-gmt-symbol",
                    help="Also save REAL+NULL Symbol GMT (requires --gene-loc).")

    ap.set_defaults(func=main)
    return ap


def main(args):
    """Main execution function."""
    # Sanity checks
    if not (os.path.isfile(args.magma) and os.access(args.magma, os.X_OK)):
        raise SystemExit(f"MAGMA not executable: {args.magma}")
    if not os.path.isfile(args.genes_raw):
        raise SystemExit(f".genes.raw not found: {args.genes_raw}")
    if not os.path.isfile(args.sets):
        raise SystemExit(f"Gene-set file not found: {args.sets}")

    # ---------- Prepare input GMT ----------
    tmp_dir = tempfile.mkdtemp(prefix="scrbp_magma_")
    combined_gmt = None
    null_index_path = None
    temp_paths = []
    short_to_full = {}   # ct mode: short ID -> full regulon name (internal, never saved)

    if args.mode == "sc":
        # sc: no need for NULL, no forced expr-stats
        if args.id_type == "entrez":
            combined_gmt = args.sets
        else:
            if not args.gene_loc:
                raise SystemExit("--gene-loc is required when --id-type=symbol")
            _, s2e = load_entrez_symbol_from_gene_loc(args.gene_loc)
            bg = read_genes_raw_head(args.genes_raw)
            bg_entrez = set(bg["GENE"].tolist())
            real_sets = read_gmt(args.sets)
            real_entrez = {}
            for name, genes in real_sets.items():
                ents = []
                for g in genes:
                    e = s2e.get(str(g).upper())
                    if e is not None and e in bg_entrez:
                        ents.append(str(e))
                if ents: real_entrez[name] = ents
            combined_gmt = os.path.join(tmp_dir, "real_entrez.gmt")
            write_gmt(combined_gmt, real_entrez); temp_paths.append(combined_gmt)

    else:  # -------- ct mode: build REAL + matched NULL --------
        bg = read_genes_raw_head(args.genes_raw)          # GENE, NSNPS, NPARAM
        bg_entrez = set(bg["GENE"].tolist())

        # Load REAL sets as Entrez
        if args.id_type == "entrez":
            real_entrez = read_gmt(args.sets)
            cleaned = {}
            for name, genes in real_entrez.items():
                ents = []
                for g in genes:
                    try: e = int(g)
                    except Exception: continue
                    if e in bg_entrez: ents.append(e)
                if ents: cleaned[name] = ents
            real_entrez = cleaned
            e2s, s2e = load_entrez_symbol_from_gene_loc(args.gene_loc) if args.gene_loc else ({}, {})
        else:
            if not args.gene_loc:
                raise SystemExit("--gene-loc is required in ct mode when --id-type=symbol")
            e2s, s2e = load_entrez_symbol_from_gene_loc(args.gene_loc)
            real_symbol = read_gmt(args.sets)
            real_entrez = {}
            for name, genes in real_symbol.items():
                ents = []
                for g in genes:
                    e = s2e.get(str(g).upper())
                    if e is not None and e in bg_entrez:
                        ents.append(e)
                if ents: real_entrez[name] = ents

        # ---- expr-stats: read or compute (only needed in ct mode) ----
        expr_stats_df = None
        if args.expr_stats and os.path.isfile(args.expr_stats):
            expr_stats_df = read_expr_stats(args.expr_stats)
            print(f"[ExprStats] Loaded from: {args.expr_stats} | genes={expr_stats_df.shape[0]}")
        elif args.emit_expr_stats:
            if not args.matrix_stats:
                print("[ExprStats][WARN] --emit-expr-stats is True but no --matrix-stats provided; "
                      "NULL matching will fall back to NSNPS/NPARAM only.")
            else:
                X = load_expression_matrix(args.matrix_stats, dtype=args.dtype)
                expr_stats_df = compute_expr_stats_from_df(X)
                outp = args.expr_stats_out if args.expr_stats_out else f"{args.out}_expr_stats.tsv"
                expr_stats_df.to_csv(outp, sep="\t", index=False)
                print(f"[ExprStats] Computed & saved: {outp} | genes={expr_stats_df.shape[0]}")

        # Merge expr-stats into bg (by symbol)
        if expr_stats_df is not None:
            if not e2s:
                if not args.gene_loc:
                    raise SystemExit("--gene-loc is required to map expr-stats symbols to Entrez in ct mode.")
                e2s, s2e = load_entrez_symbol_from_gene_loc(args.gene_loc)
            map_df = pd.DataFrame({"GENE": list(e2s.keys()),
                                   "symbol": [e2s[e] for e in e2s.keys()]})
            map_df["symbol"] = map_df["symbol"].astype(str).str.upper()
            bg = bg.merge(map_df, on="GENE", how="left")
            bg = bg.merge(expr_stats_df, on="symbol", how="left")
            for c in ("mean_expr", "pct_detected"):
                if c in bg.columns:
                    bg[c] = pd.to_numeric(bg[c], errors="coerce").fillna(0.0)
            have_expr_cov = True
        else:
            have_expr_cov = False
            bg["mean_expr"] = 0.0
            bg["pct_detected"] = 0.0

        rng = np.random.default_rng(args.seed)
        sets_entrez = {}        # full name -> gene list (for saving / indexing)
        _sets_short  = {}       # short ID  -> gene list (for MAGMA only)
        index_rows = []
        _sid = 0                # short-ID counter

        for rbp, ent_list in real_entrez.items():
            ent_list = list(dict.fromkeys(int(x) for x in ent_list))
            full_name = f"{rbp}__REAL"
            genes_str = [str(e) for e in ent_list]
            sets_entrez[full_name] = genes_str

            # Assign a short ID that MAGMA cannot truncate (max 9 chars)
            _sid += 1
            sid = f"GS{_sid:07d}"
            _sets_short[sid] = genes_str
            short_to_full[sid] = full_name

            nulls = make_null_for_regulon_4d(
                ent_list, bg, n_null=args.n_null, q=args.q_bins,
                rng=rng, exclude_self=args.exclude_self
            )
            for i, nl in enumerate(nulls, 1):
                full_name = f"{rbp}__NULL_{i:04d}"
                genes_str = [str(e) for e in nl]
                sets_entrez[full_name] = genes_str
                index_rows.append((rbp, i, full_name, len(nl)))

                _sid += 1
                sid = f"GS{_sid:07d}"
                _sets_short[sid] = genes_str
                short_to_full[sid] = full_name

        # Write SHORT-ID GMT for MAGMA (avoids fixed-width column truncation)
        combined_gmt = os.path.join(tmp_dir, "real_plus_nulls.entrez.gmt")
        write_gmt(combined_gmt, _sets_short)
        temp_paths.append(combined_gmt)
        print(f"[SHORT-ID] Mapped {len(short_to_full)} gene sets to short IDs "
              f"(GS0000001..GS{_sid:07d}) for MAGMA")

        # 保存 NULL 索引
        null_index_path = f"{args.out}.null_index.tsv"
        pd.DataFrame(index_rows, columns=["RBP", "NULL_ID", "SET_NAME", "SIZE"]).to_csv(
            null_index_path, sep="\t", index=False
        )
        print("Wrote:", null_index_path)

        # ------- Save REAL+NULL GMT (full names, not short IDs) -------
        # 1) Entrez GMT
        save_entrez_path = args.save_null_gmt if args.save_null_gmt \
            else f"{args.out}_REAL_PLUS_NULLS.entrez.gmt"
        write_gmt(save_entrez_path, sets_entrez)
        print(f"[SAVE] REAL+NULL Entrez GMT -> {save_entrez_path}")

        # 2) Symbol GMT (optional)
        if args.save_null_gmt_symbol:
            if not args.gene_loc and not e2s:
                raise SystemExit("--save-null-gmt-symbol requires --gene-loc to map Entrez->Symbol.")
            if not e2s:
                e2s, _ = load_entrez_symbol_from_gene_loc(args.gene_loc)
            # convert Entrez sets to Symbol
            sets_symbol = {}
            for name, genes in sets_entrez.items():
                syms = []
                for g in genes:
                    try:
                        e = int(g)
                        s = e2s.get(e)
                        if s: syms.append(str(s))
                    except Exception:
                        pass
                if syms:
                    sets_symbol[name] = syms
            save_sym_path = args.save_null_gmt_symbol
            write_gmt(save_sym_path, sets_symbol)
            print(f"[SAVE] REAL+NULL Symbol GMT -> {save_sym_path}")

    # ---------- Run MAGMA ----------
    cmd = [args.magma, "--gene-results", args.genes_raw, "--set-annot", combined_gmt, "--out", args.out]
    if args.threads: cmd += ["--threads", str(args.threads)]
    if args.extra:   cmd += shlex.split(args.extra)

    print("Running:", " ".join(shlex.quote(c) for c in cmd))
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(proc.stdout)
    if proc.returncode != 0:
        raise SystemExit(f"[MAGMA] GSA failed with code {proc.returncode}")

    # ---------- Parse .gsa.out ----------
    gsa_out = f"{args.out}.gsa.out"
    df = read_gsa_out(gsa_out, one_sided=True)

    # ---- Restore full regulon names from short IDs (ct mode) ----
    if short_to_full:
        n_before = len(df)
        n_mapped = df["Regulon"].isin(short_to_full).sum()
        df["Regulon"] = df["Regulon"].map(short_to_full).fillna(df["Regulon"])
        print(f"[SHORT-ID] Restored {n_mapped}/{n_before} regulon names from short IDs")

    if args.min_genes and "NGENES" in df.columns:
        before = len(df)
        df = df[df["NGENES"].fillna(0) >= int(args.min_genes)].copy()
        after = len(df)
        print(f"[FILTER] min_genes={args.min_genes}: kept {after}/{before} regulons.")
    elif args.min_genes and "NGENES" not in df.columns:
        print("[FILTER][WARN] NGENES column not found in .gsa.out; --min_genes ignored.")

    if args.mode == "ct":
        splits = df["Regulon"].apply(split_regulon_name)
        df["RBP"]      = splits.apply(lambda x: x[0])
        df["SET_KIND"] = splits.apply(lambda x: x[1])
        df["NULL_ID"]  = splits.apply(lambda x: x[2])

    out_csv = f"{args.out}.gsa_RGS.csv"
    df.to_csv(out_csv, index=False)
    print("Wrote:", out_csv)

    # ---------- Cleanup ----------
    if args.cleanup_out:
        removed = 0
        for f in glob.glob(f"{args.out}*.out"):
            try: os.remove(f); removed += 1
            except OSError as e: print(f"[CLEANUP][WARN] failed to remove {f}: {e}")
        # do NOT delete explicitly saved GMT files
        for p in temp_paths:
            try: os.remove(p)
            except OSError: pass
        try: os.rmdir(tmp_dir)
        except OSError: pass
        if removed == 0:
            print(f"[CLEANUP] no files matched: {args.out}*.out")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    subp = p.add_subparsers()
    register_subcommand(subp)
    args = p.parse_args()
    main(args)
