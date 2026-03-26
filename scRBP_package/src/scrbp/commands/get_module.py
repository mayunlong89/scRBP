#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scRBP getModule — Extract RBP regulon candidate modules from a consensus GRN

Applies three complementary selection strategies to identify high-confidence
RBP–target edges from the consensus gene regulatory network (GRN):

  TopN-per-gene      For each target gene, retain the top-N highest-importance RBPs.
  TopN-per-RBP       For each RBP, retain its top-N highest-importance target genes.
  Percentile-per-RBP For each RBP, retain targets above the Xth importance percentile.

All strategies are applied independently and merged into a single output file
with a ``Strategy`` column labelling each edge's origin.  The resulting modules
serve as foreground gene sets for downstream motif-enrichment pruning (getPrune).
"""

import sys
import argparse
import warnings
import re

import pandas as pd


def _parse_int_list(s, name, min_val=1):
    if isinstance(s, (list, tuple)):
        vals = s
    else:
        parts = re.split(r"[,\s]+", str(s).strip())
        vals = [p for p in parts if p != ""]
    try:
        out = sorted({int(x) for x in vals})
    except Exception:
        raise ValueError(f"--{name} needs integer list, e.g.: 5,10,50")
    bad = [v for v in out if v < min_val]
    if bad:
        raise ValueError(f"--{name} values must be >= {min_val}, found: {bad}")
    return out


def _parse_float_list(s, name, lo=0.0, hi=1.0):
    parts = re.split(r"[,\s]+", str(s).strip())
    vals = [p for p in parts if p != ""]
    try:
        out = sorted({float(x) for x in vals})
    except Exception:
        raise ValueError(f"--{name} needs float list, e.g.: 0.75,0.9,0.95")
    bad = [v for v in out if not (lo < v <= hi)]
    if bad:
        raise ValueError(f"--{name} values must be in ({lo}, {hi}], found: {bad}")
    return out


def load_data(input_path, threshold, verbose=False):
    try:
        df = pd.read_csv(input_path, sep="\t")
        for col in ["RBP", "Gene", "Importance"]:
            if col not in df.columns:
                raise ValueError(f"Input file missing required column: {col}")
        if verbose:
            print(f"Loaded {len(df)} records from {input_path}")
        return df[df['Importance'] > threshold]
    except Exception as e:
        print(f"Failed to load input file: {e}")
        sys.exit(1)


def save_top_n_per_rbp(df, n, strategy_name, verbose):
    result = (df.sort_values(['RBP', 'Importance'], ascending=[True, False])
                .groupby('RBP', group_keys=False).head(n).copy())
    result['Strategy'] = strategy_name
    if verbose:
        print(f"Processed Top-{n} targets per RBP")
    return result


def save_top_n_per_gene(df, n, strategy_name, verbose):
    result = (df.sort_values(['Gene', 'Importance'], ascending=[True, False])
                .groupby('Gene', group_keys=False).head(n).copy())
    result['Strategy'] = strategy_name
    if verbose:
        print(f"Processed Top-{n} RBPs per gene")
    return result


def save_percentile_per_rbp(df, percentile, strategy_name, verbose):
    def filter_group(group):
        thr = group['Importance'].quantile(percentile)
        return group[group['Importance'] >= thr]
    result = df.groupby('RBP', group_keys=False).apply(filter_group).reset_index(drop=True)
    result['Strategy'] = strategy_name
    if verbose:
        print(f"Processed Top-{int(percentile*100)} percentile targets per RBP")
    return result


def register_subcommand(subparsers):
    """Register this subcommand with the main parser."""
    parser = subparsers.add_parser(
        'getModule',
        help='Extract GRN submodules using multiple strategies',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )
    parser.add_argument('--input', required=True, help='Input GRN file (TSV)')
    parser.add_argument('--importance_threshold', type=float, default=0.005)
    # All three strategy parameters accept comma- or space-separated lists
    parser.add_argument('--top_n_list', type=str, default="5,10,50",
                        help='Top-N RBPs per gene, e.g. "5,10,50,100"')
    parser.add_argument('--target_top_n', type=str, default="50",
                        help='Top-N targets per RBP, e.g. "50,100,200"')
    parser.add_argument('--percentile', type=str, default="0.75,0.9",
                        help='Percentiles per RBP, e.g. "0.75,0.9,0.95"')
    parser.add_argument('--output_merged', required=True, help='Output merged file name (TSV)')
    parser.add_argument('--verbose', action='store_true')
    parser.set_defaults(func=main)
    return parser


def main(args):
    """Main execution function."""
    warnings.simplefilter(action='ignore', category=DeprecationWarning)

    # Parse multi-value strategy parameters
    top_n_list = _parse_int_list(args.top_n_list, "top_n_list", min_val=1)
    target_top_n_list = _parse_int_list(args.target_top_n, "target_top_n", min_val=1)
    percentile_list = _parse_float_list(args.percentile, "percentile", lo=0.0, hi=1.0)

    if args.verbose:
        print(f"Top-N per gene: {top_n_list}")
        print(f"Top-N per RBP:  {target_top_n_list}")
        print(f"Percentiles:    {percentile_list}")

    # Load and filter data
    df_filtered = load_data(args.input, args.importance_threshold, args.verbose)

    all_modules = []

    # Strategy 1: top-N highest-importance RBPs per target gene
    for n in top_n_list:
        strategy_name = f"Top{n}perTarget"
        all_modules.append(save_top_n_per_gene(df_filtered, n, strategy_name, args.verbose))

    # Strategy 2: top-N highest-importance target genes per RBP
    for n in target_top_n_list:
        strategy_name = f"Top{n}perRBP"
        all_modules.append(save_top_n_per_rbp(df_filtered, n, strategy_name, args.verbose))

    # Strategy 3: importance percentile threshold per RBP
    for p in percentile_list:
        strategy_name = f"TopPercentile{int(p*100)}"
        all_modules.append(save_percentile_per_rbp(df_filtered, p, strategy_name, args.verbose))

    merged_df = pd.concat(all_modules, ignore_index=True)
    # Retain all available columns in canonical order
    cols = [c for c in ['RBP', 'Gene', 'Importance', 'Correlation', 'Mode', 'Strategy'] if c in merged_df.columns]
    final_df = merged_df[cols]
    final_df.to_csv(args.output_merged, sep="\t", index=False)
    print(f"All GRN submodule strategies merged and saved to {args.output_merged}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    subp = p.add_subparsers()
    register_subcommand(subp)
    args = p.parse_args()
    main(args)
