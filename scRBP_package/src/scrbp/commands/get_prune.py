#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scRBP getPrune — Rank-based motif-enrichment pruning of RBP regulon candidates

Refines co-expression-derived RBP–target modules by applying a rank-based
motif-enrichment framework (via ctxcore) that retains only high-fidelity
direct RBP–target interactions supported by physical RBP-binding evidence.

For each RBP regulon candidate module (from getModule), enrichment against
precomputed transcriptome-wide motif-target rankings is evaluated.  Only
targets ranked above the leading-edge index — defined as the rank position
where the deviation between the motif-specific recovery curve and the
background reference exceeds two standard deviations — are retained as
direct targets in the final regulon.  Modules without an enriched matched
motif are discarded.

Output: a merged CSV (scRBP_prune_matched_results.csv) containing motif
enrichment statistics (NES, AUC) and leading-edge gene sets for each retained
RBP–strategy combination.

Auto-detects ranking feather index column (``features`` or ``motifs``) and
motif-regulator annotation column (``RBP`` or ``TF``).
"""

import os
import argparse
import multiprocessing as mp
from collections import defaultdict

import numpy as np
import pandas as pd
from tqdm import tqdm

from ctxcore.genesig import GeneSignature
from ctxcore.recovery import enrichment4features, leading_edge4row
from ctxcore.rnkdb import RankingDatabase


# =========================
# Internal utilities
# =========================
def _flatten_cols(cols):
    """Flatten column names (possibly MultiIndex tuples) to plain strings for Parquet/CSV output."""
    flat = []
    for c in cols:
        if isinstance(c, tuple):
            flat.append("-".join(map(str, c)))
        else:
            flat.append(str(c))
    return flat


def _find_nes_col(df):
    """Locate the NES column after flattening. Prefers 'Enrichment-NES'; falls back to any column ending in 'nes'."""
    if "Enrichment-NES" in df.columns:
        return "Enrichment-NES"
    cands = [c for c in df.columns if isinstance(c, str) and c.lower().endswith("nes")]
    return cands[0] if cands else None


def _safe_for_parquet(df: pd.DataFrame) -> pd.DataFrame:
    """Convert object-dtype columns to strings (list/ndarray → comma-joined) to prevent pyarrow serialisation errors."""
    out = df.copy()
    for col in out.columns:
        if pd.api.types.is_object_dtype(out[col]):
            out[col] = out[col].map(
                lambda x: ",".join(map(str, x)) if isinstance(x, (list, tuple, np.ndarray))
                else ("" if x is None else str(x))
            ).astype("string")
    return out


def read_ranking_feather(
    path: str,
    index_candidates=("features", "motifs", "motif", "motif_id", "Motif")
) -> pd.DataFrame:
    """
    Read motif×gene ranking feather and set motif id as index.
    Supports both RBP-style (features) and TF-style (motifs).
    """
    df = pd.read_feather(path)

    # 1) find index column (exact match)
    idx_col = None
    for c in index_candidates:
        if c in df.columns:
            idx_col = c
            break

    # 2) case-insensitive match
    if idx_col is None:
        cols_lower = {str(c).lower(): c for c in df.columns}
        for c in index_candidates:
            key = str(c).lower()
            if key in cols_lower:
                idx_col = cols_lower[key]
                break

    if idx_col is None:
        raise SystemExit(
            f"[ERROR] Cannot find motif id column in ranking feather: {path}\n"
            f"Expect one of {list(index_candidates)}, but got columns like: {list(df.columns)[:30]} ..."
        )

    print(f"[INFO] ranking index column detected: {idx_col}")

    df[idx_col] = df[idx_col].astype(str).str.strip()
    df = df.set_index(idx_col)

    # gene columns as strings
    df.columns = df.columns.astype(str)

    return df


def read_motif_regulator_links(path: str):
    """
    Read motif-regulator links (annotation) file.
    Supports columns like:
      - Motif,RBP
      - motif,rbp
      - motifs,TF
      - motif_id,tf
      - motif_id,regulator
    Returns: (normalized_df, regulator_column_name)
      normalized_df columns include: motif_id, <reg_col>
    """
    # read with auto delimiter
    try:
        df = pd.read_csv(path)  # common CSV
    except Exception:
        df = pd.read_csv(path, sep=None, engine="python")

    # normalize headers
    df.columns = (
        df.columns.astype(str).str.strip()
        .str.replace(r"\s+", "_", regex=True).str.lower()
    )

    # find motif column
    if "motif_id" not in df.columns:
        for cand in ("motif", "motifs", "feature", "features"):
            if cand in df.columns:
                df = df.rename(columns={cand: "motif_id"})
                break

    # find regulator column
    reg_col = None
    for cand in ("rbp", "tf", "regulator", "factor", "gene", "regulator_name"):
        if cand in df.columns:
            reg_col = cand
            break

    if "motif_id" not in df.columns or reg_col is None:
        raise SystemExit(
            f"[ERROR] Annotation file must contain motif column (motif/motifs/features/motif_id) "
            f"and regulator column (rbp/tf/regulator/gene). Got columns: {list(df.columns)}"
        )

    # clean values
    df["motif_id"] = df["motif_id"].astype(str).str.strip()
    df[reg_col] = df[reg_col].astype(str).str.strip().str.upper()

    print(f"[INFO] annotation columns detected: motif_id + {reg_col}")

    return df, reg_col


# =========================
# In-memory RankingDatabase
# =========================
class SimpleMemoryRankingDatabase(RankingDatabase):
    def __init__(self, name, df):
        super().__init__(name)
        self._df = df

    @property
    def total_genes(self):
        return self._df.shape[1]

    @property
    def genes(self):
        return tuple(self._df.columns)

    def load_full(self):
        return self._df

    def load(self, gs):
        # Subset to only the columns (genes) required by this GeneSignature to reduce memory transfer
        return self._df.loc[:, self._df.columns.intersection(gs.genes)]


# =========================
# Build GeneSignature list
# =========================
def build_gene_signatures(
    df_modules: pd.DataFrame,
    db,
    weight_mode: str = "importance",   # "importance" or "equal"
    min_genes: int = 10,
) -> list:
    """
    Build a list of GeneSignature objects from an RBP-target module table.

    Parameters
    ----------
    df_modules : pd.DataFrame
        Module table with columns ``RBP``, ``Gene``, ``Strategy``,
        and optionally ``Importance``.
    db : RankingDatabase
        Motif-target ranking database used to filter to ranked genes.
    weight_mode : str
        ``"importance"`` — use the ``Importance`` column as gene weights;
        falls back to equal weighting if the column is absent.
        ``"equal"``      — assign weight 1.0 to all genes.
    min_genes : int
        Minimum number of genes (after filtering to ranked genes) required
        to retain a (RBP, Strategy) module.
    """
    weight_mode = (weight_mode or "equal").lower()
    assert weight_mode in {"importance", "equal"}, \
        f"weight_mode must be 'importance' or 'equal', got: {weight_mode}"

    gene_signatures = []
    genes_in_db = set(db.genes)

    use_importance = (weight_mode == "importance") and ("Importance" in df_modules.columns)

    if (weight_mode == "importance") and not use_importance:
        print("[WARN] weight_mode='importance' selected but 'Importance' column not found; falling back to equal.")

    # Build one GeneSignature per (RBP, Strategy) pair
    for (rbp, strategy), sub_df in df_modules.groupby(["RBP", "Strategy"], sort=False):
        if use_importance:
            w = pd.to_numeric(sub_df["Importance"], errors="coerce").fillna(0.0).astype(float).values
            g = sub_df["Gene"].astype(str).values
            gene2weight_all = dict(zip(g, w))
        else:
            gene2weight_all = {g: 1.0 for g in sub_df["Gene"].astype(str).values}

        # Retain only genes present in the ranking database
        filtered = {g: w for g, w in gene2weight_all.items() if g in genes_in_db}

        if len(filtered) < min_genes:
            continue

        gs_obj = GeneSignature(name=f"{rbp}__{strategy}", gene2weight=filtered)
        gene_signatures.append(gs_obj)

    return gene_signatures


# =========================
# Parallel worker (subprocess)
# =========================
_G_DB = None           # Ranking database shared across worker processes
_G_MOTIF2REGS = None   # motif -> set(REGULATOR) annotation map; REGULATOR names are uppercased


def _init_worker(db, motif2regs):
    global _G_DB, _G_MOTIF2REGS
    _G_DB = db
    _G_MOTIF2REGS = motif2regs


def _process_one_signature(gs_obj: GeneSignature,
                           rank_threshold: int,
                           auc_threshold: float,
                           nes_threshold: float,
                           part_dir: str,
                           part_id: int):
    """
    Worker function: run enrichment + leading-edge extraction + NES filtering +
    motif-regulator annotation matching for one GeneSignature; write results to
    a Parquet shard file.

    Returns
    -------
    tuple : (part_path or None, error_msg or None)
    """
    try:
        df_enrich = enrichment4features(
            _G_DB, gs_obj, rank_threshold=rank_threshold, auc_threshold=auc_threshold
        )
        if df_enrich.empty:
            rbp, strategy = gs_obj.name.split("__", 1)
            return None, f"{rbp}__{strategy}: empty enrichment"

        # Background reference curve: mean recovery + 2 × std across all motifs
        recovery_matrix = np.stack(df_enrich["Recovery"].values)
        avg2stdrcc = np.mean(recovery_matrix, axis=0) + 2 * np.std(recovery_matrix, axis=0)

        # Compute leading-edge genes for each enriched motif
        ranking_cols = [c for c in df_enrich.columns if isinstance(c, tuple) and c[0] == "Ranking"]
        genes = [c[1] for c in ranking_cols]
        df_le = df_enrich.apply(
            leading_edge4row,
            axis=1,
            args=(),
            avg2stdrcc=avg2stdrcc,
            genes=np.array(genes),
        )

        # Keep only NES/AUC enrichment scores and leading-edge metrics
        keep_cols = []
        for col in [("Enrichment", "NES"), ("Enrichment", "AUC")]:
            if col in df_enrich.columns:
                keep_cols.append(col)
        df_small = pd.concat([df_enrich[keep_cols], df_le], axis=1)

        # motif_id 列
        df_small = df_small.reset_index().rename(columns={"index": "motif_id"})
        df_small["motif_id"] = df_small["motif_id"].astype(str).str.strip()

        # 扁平列名
        df_small.columns = _flatten_cols(df_small.columns)

        # NES 过滤
        nes_col = _find_nes_col(df_small)
        if nes_col is None:
            rbp, strategy = gs_obj.name.split("__", 1)
            return None, f"{rbp}__{strategy}: NES column not found"
        df_small = df_small[df_small[nes_col] > nes_threshold]
        if df_small.empty:
            rbp, strategy = gs_obj.name.split("__", 1)
            return None, f"{rbp}__{strategy}: empty after NES>{nes_threshold}"

        # Annotate with RBP and strategy labels; regulator names uppercased for matching
        rbp, strategy = gs_obj.name.split("__", 1)
        reg_key = rbp.strip().upper()  # regulator name for annotation lookup (TF names also stored in RBP column)
        df_small["RBP"] = rbp          # output column named RBP for downstream compatibility
        df_small["Strategy"] = strategy

        # 注释匹配：当前 regulator 是否在 motif 的注释集合中
        mask = df_small["motif_id"].map(lambda m: reg_key in _G_MOTIF2REGS.get(m, ()))
        df_small = df_small[mask]
        if df_small.empty:
            return None, f"{rbp}__{strategy}: no matched motifs after annotation filter"

        # Sanitise types before writing to Parquet (prevents pyarrow serialisation errors)
        df_small = _safe_for_parquet(df_small)
        part_path = os.path.join(part_dir, f"part_{part_id:06d}.parquet")
        df_small.to_parquet(part_path, index=False)
        return part_path, None

    except Exception as e:
        rbp, strategy = gs_obj.name.split("__", 1)
        return None, f"{rbp}__{strategy}: {repr(e)}"


def _worker_task(args):
    """Top-level picklable wrapper around _process_one_signature (replaces lambda/partial)."""
    (gs_obj, part_id, rank_threshold, auc_threshold, nes_threshold, part_dir) = args
    return _process_one_signature(
        gs_obj,
        rank_threshold=rank_threshold,
        auc_threshold=auc_threshold,
        nes_threshold=nes_threshold,
        part_dir=part_dir,
        part_id=part_id,
    )


def register_subcommand(subparsers):
    """Register this subcommand with the main parser."""
    parser = subparsers.add_parser(
        'getPrune',
        help='Motif-based network pruning with ctxcore',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )
    parser.add_argument("--rbp_targets", required=True,
                        help="TSV/CSV: columns=[RBP, Gene, Strategy, (optional) Importance]. "
                             "NOTE: For TF workflow, you can still use column name 'RBP' for TF names.")
    parser.add_argument("--motif_rbp_links", required=True,
                        help="TSV/CSV: motif-regulator links. Supports columns like Motif,RBP or Motif,TF.")
    parser.add_argument("--motif_target_ranks", required=True,
                        help="Feather: motif x gene ranking matrix (index column can be 'features' or 'motifs')")
    parser.add_argument("--save_dir", required=True, help="Output directory")
    parser.add_argument("--rank_threshold", type=int, default=1500)
    parser.add_argument("--auc_threshold", type=float, default=0.05)
    parser.add_argument("--min_genes", type=int, default=20)
    parser.add_argument("--nes_threshold", type=float, default=3.0)
    parser.add_argument("--n_jobs", type=int, default=os.cpu_count())
    parser.add_argument("--chunksize", type=int, default=4, help="imap_unordered chunksize")
    # 调试辅助
    parser.add_argument("--only_rbp", type=str, default=None)
    parser.add_argument("--only_strategy", type=str, default=None)
    parser.add_argument("--max_sigs", type=int, default=None)
    parser.add_argument("--keep_parts", action="store_true", help="Keep parquet parts (for debugging)")
    parser.add_argument(
        "--weight_mode",
        type=str,
        choices=["importance", "equal"],
        default="importance",
        help="Gene weight mode: importance=use Importance column, equal=all weights equal"
    )
    parser.set_defaults(func=main)
    return parser


def main(args):
    """Main execution function."""
    os.makedirs(args.save_dir, exist_ok=True)
    part_dir = os.path.join(args.save_dir, "parts_matched")
    os.makedirs(part_dir, exist_ok=True)

    # Step 1: load motif-target ranking matrix and downcast dtypes to reduce memory usage
    df_rank = read_ranking_feather(args.motif_target_ranks)

    for c in df_rank.columns:
        s = df_rank[c]
        if pd.api.types.is_integer_dtype(s):
            df_rank[c] = pd.to_numeric(s, downcast="integer")
        elif pd.api.types.is_float_dtype(s):
            df_rank[c] = pd.to_numeric(s, downcast="float")

    db = SimpleMemoryRankingDatabase("motifdb", df_rank)

    # Step 2: load RBP-target modules (TF workflows can reuse this format with TF names in the RBP column)
    try:
        df_modules = pd.read_csv(args.rbp_targets, sep="\t", engine="python")
    except Exception:
        df_modules = pd.read_csv(args.rbp_targets, sep=None, engine="python")

    # Step 3: load motif-regulator annotation and build motif_id -> set(REGULATOR) map; check coverage
    df_annot, reg_col = read_motif_regulator_links(args.motif_rbp_links)

    # Pre-check: motif overlap between ranking matrix and annotation file
    rank_motifs = pd.Index(df_rank.index.astype(str).str.strip())
    annot_motifs = pd.Index(df_annot["motif_id"].unique())
    overlap = rank_motifs.intersection(annot_motifs)
    print(f"[CHECK] motifs in ranking={len(rank_motifs)}, in annotation={len(annot_motifs)}, overlap={len(overlap)}")
    if len(overlap) == 0:
        raise SystemExit(
            "No overlap between ranking index and annotation 'Motif'. "
            "Please check that both sides use the same Motif ID set."
        )

    # Build motif -> REGULATOR set mapping (regulator names uppercased for consistent matching)
    motif2regs = defaultdict(set)
    for m, r in zip(df_annot["motif_id"], df_annot[reg_col]):
        if m and r and m.lower() != "nan" and r.lower() != "nan":
            motif2regs[m].add(r)

    # Step 4: build GeneSignature objects for all (RBP, Strategy) modules
    gene_signatures = build_gene_signatures(
        df_modules,
        db,
        weight_mode=args.weight_mode,
        min_genes=args.min_genes
    )
    if not gene_signatures:
        raise SystemExit("No GeneSignature built -- check columns/IDs/min_genes.")

    # Debug filters: restrict to a specific RBP/strategy or limit total signatures
    if args.only_rbp or args.only_strategy:
        keep = []
        for gs_obj in gene_signatures:
            rbp, strat = gs_obj.name.split("__", 1)
            if args.only_rbp and rbp != args.only_rbp:
                continue
            if args.only_strategy and strat != args.only_strategy:
                continue
            keep.append(gs_obj)
        gene_signatures = keep
        print(f"[DEBUG] filtered to {len(gene_signatures)} signatures for RBP={args.only_rbp}, strategy={args.only_strategy}")
    if args.max_sigs:
        gene_signatures = gene_signatures[:args.max_sigs]
        print(f"[DEBUG] limiting to first {len(gene_signatures)} signatures")

    # Step 5: parallel enrichment — NES filtering + annotation matching + shard writing
    errors = []
    parts = []
    try:
        ctx = mp.get_context("fork")
    except ValueError:
        ctx = mp.get_context()  # Windows/macOS spawn fallback

    with ctx.Pool(processes=args.n_jobs, initializer=_init_worker, initargs=(db, motif2regs)) as pool:
        task_iter = (
            (gs_obj, i, args.rank_threshold, args.auc_threshold, args.nes_threshold, part_dir)
            for i, gs_obj in enumerate(gene_signatures)
        )
        it = pool.imap_unordered(_worker_task, task_iter, chunksize=args.chunksize)

        for part_path, err in tqdm(it, total=len(gene_signatures),
                                   desc=f"Enrichment->LE->NES->match (n_jobs={args.n_jobs})"):
            if err:
                errors.append(err)
                continue
            if part_path:
                parts.append(part_path)

    if errors:
        with open(os.path.join(args.save_dir, "errors.log"), "w") as f:
            f.write("\n".join(errors))
        print(f"[WARN] {len(errors)} GeneSignatures had warnings/errors. See errors.log")

    if not parts:
        raise SystemExit("No matched parts written -- all signatures filtered out by NES/annotation.")

    # Step 6: aggregate all matched shards into the final CSV (streamed to avoid peak memory)
    out_csv = os.path.join(args.save_dir, "scRBP_prune_matched_results.csv")
    wrote_header = False
    for pth in sorted(parts):
        d = pd.read_parquet(pth)

        # Rename numeric columns from ctxcore leading-edge output ("0", "1") to descriptive names
        d.columns = d.columns.map(str)
        d = d.rename(columns={"0": "LeadingEdge-Genes", "1": "LeadingEdge-RankAtMax"})

        # desired order
        desired = [
            "motif_id",
            "RBP",  # NOTE: even TF names will be stored here to keep compatibility with downstream scripts
            "Enrichment-AUC",
            "Enrichment-NES",
            "Strategy",
            "LeadingEdge-RankAtMax",
            "LeadingEdge-Genes",
        ]

        front = [c for c in desired if c in d.columns]
        tail = [c for c in d.columns if c not in front]
        d = d[front + tail]

        d.to_csv(out_csv, index=False, mode="a", header=not wrote_header)
        wrote_header = True

    # 清理分片
    if not args.keep_parts:
        for p in parts:
            try:
                os.remove(p)
            except Exception:
                pass
        try:
            os.rmdir(part_dir)
        except Exception:
            pass

    print(f"Saved final matched results to: {out_csv}")


if __name__ == "__main__":
    try:
        mp.set_start_method("fork")
    except RuntimeError:
        pass
    p = argparse.ArgumentParser()
    subp = p.add_subparsers()
    register_subcommand(subp)
    args = p.parse_args()
    main(args)
