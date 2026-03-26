#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scRBP getMerge_GRN — Merge multi-seed GRN runs into a consensus RBP regulon network

scRBP recommends running getGRN with multiple random seeds (e.g., 30 seeds) to
reduce the stochasticity inherent in tree-based network inference.  This command
aggregates all seed runs into a single consensus network by averaging edge scores
across runs and computing cross-seed stability metrics.

For each RBP–Gene edge:
  - Importance and Correlation are averaged across all runs (edges absent in a
    given run contribute 0 to the average, i.e., missing = 0 convention).
  - Mode (activating / repressing / unknown) is re-assigned from the averaged
    Spearman correlation using the shared assign_mode utility.
  - n_present counts the number of seed runs containing the edge.
  - presence_rate = n_present / N_runs quantifies cross-seed stability.

Optional filters allow retaining only high-confidence edges:
  --corr-threshold   Remove edges where |averaged Correlation| <= threshold.
  --n_present        Keep edges present in at least this many seed runs.
  --present_rate     Keep edges with presence_rate >= this value (e.g., 0.5 = 50% of runs).
"""

import os
import glob
import argparse

import numpy as np
import pandas as pd

from scrbp.utils.stats import assign_mode  # shared with get_grn


def register_subcommand(subparsers):
    """Register this subcommand with the main parser."""
    p = subparsers.add_parser(
        'getMerge_GRN',
        help='Merge multiple GRN runs into consensus network',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )
    p.add_argument(
        "--pattern",
        required=True,
        help=(
            "Glob pattern to match input GRN files, e.g. "
            "'z_GRNBoost2_8autoimmunediseases_170Kcells_seed*_scRBP_GRNs.tsv'"
        ),
    )
    p.add_argument(
        "--output",
        required=True,
        help="Output TSV file for merged GRN, e.g. 'merged_scRBP_GRNs.tsv'.",
    )
    p.add_argument(
        "--corr-threshold",
        type=float,
        default=0.0,
        help=(
            "Filter edges with abs(averaged Correlation) <= threshold. "
            "Default: 0 (no additional filtering)."
        ),
    )
    p.add_argument(
        "--n_present",
        type=int,
        default=0,
        help=(
            "Minimum number of runs in which an edge must appear "
            "(n_present >= this value). Default: 0 (no filtering)."
        ),
    )
    p.add_argument(
        "--present_rate",
        type=float,
        default=0.0,
        help=(
            "Minimum presence rate (n_present / N_runs) required to "
            "keep an edge. Default: 0.0 (no filtering)."
        ),
    )
    p.set_defaults(func=main)
    return p


def main(args):
    """Main execution function."""
    files = sorted(glob.glob(args.pattern))
    if not files:
        raise SystemExit(f"[ERROR] No files matched pattern: {args.pattern}")

    print(f"[INFO] Found {len(files)} GRN files:")
    for f in files:
        print("  -", os.path.basename(f))

    dfs = []
    for f in files:
        df = pd.read_csv(f, sep="\t")
        missing_cols = {"RBP", "Gene", "Importance", "Correlation"} - set(df.columns)
        if missing_cols:
            raise ValueError(
                f"File {f} is missing required columns: {sorted(missing_cols)}"
            )
        df = df[["RBP", "Gene", "Importance", "Correlation"]].copy()
        df["run"] = os.path.basename(f)
        dfs.append(df)

    all_df = pd.concat(dfs, ignore_index=True)
    n_runs = len(files)
    print(f"[INFO] Total edges across all runs: {all_df.shape[0]}")
    print(f"[INFO] Number of runs (seeds): {n_runs}")

    # Sum + count for each edge
    grouped = (
        all_df.groupby(["RBP", "Gene"], as_index=False)
        .agg(
            Importance_sum=("Importance", "sum"),
            Correlation_sum=("Correlation", "sum"),
            n_present=("Importance", "size"),
        )
    )

    # Presence rate
    grouped["presence_rate"] = grouped["n_present"] / float(n_runs)

    # Average across all runs (missing edges treated as 0)
    grouped["Importance"] = grouped["Importance_sum"] / float(n_runs)
    grouped["Correlation"] = grouped["Correlation_sum"] / float(n_runs)

    # Mode from averaged Correlation
    grouped["Mode"] = grouped["Correlation"].apply(assign_mode)

    # Filter 1: remove edges with weak averaged correlation
    if args.corr_threshold > 0:
        before = grouped.shape[0]
        grouped = grouped[grouped["Correlation"].abs() > args.corr_threshold]
        after = grouped.shape[0]
        print(
            f"[INFO] Filtered edges with |Correlation| <= {args.corr_threshold}: "
            f"{before} -> {after}"
        )

    # Filter 2: remove edges with insufficient seed-run count
    if args.n_present > 0:
        before = grouped.shape[0]
        grouped = grouped[grouped["n_present"] >= args.n_present]
        after = grouped.shape[0]
        print(
            f"[INFO] Filtered edges with n_present < {args.n_present}: "
            f"{before} -> {after}"
        )

    # Filter 3: remove edges below presence rate threshold
    if args.present_rate > 0:
        before = grouped.shape[0]
        grouped = grouped[grouped["presence_rate"] >= args.present_rate]
        after = grouped.shape[0]
        print(
            f"[INFO] Filtered edges with presence_rate < {args.present_rate}: "
            f"{before} -> {after}"
        )

    # Arrange final output columns and sort by averaged importance
    out_df = grouped[
        ["RBP", "Gene", "Importance", "Correlation",
         "Mode", "n_present", "presence_rate"]
    ].copy()
    out_df = out_df.sort_values("Importance", ascending=False)

    out_df.to_csv(args.output, sep="\t", index=False)
    print(f"[INFO] Merged GRN saved to: {args.output}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    subp = p.add_subparsers()
    register_subcommand(subp)
    args = p.parse_args()
    main(args)
