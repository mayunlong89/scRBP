#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scRBP — Single-cell RNA Binding Protein regulon inference

Command-line interface entry point.  Dispatches to one of ten subcommands:

  getSketch      Stratified GeoSketch downsampling of large single-cell datasets
  getGRN         RBP → Gene/Isoform network inference (GRNBoost2 or GENIE3)
  getMerge_GRN   Merge multi-seed GRN runs into a consensus network
  getModule      Extract regulon candidate modules (Top-N / percentile)
  getPrune       Rank-based motif-enrichment pruning via ctxcore
  getRegulon     Export pruned regulons to GMT format (symbol + Entrez)
  mergeRegulons  Merge region-specific GMT files across transcript regions
  ras            Regulon Activity Score (RAS) computation via AUCell
  rgs            Regulon-level Genetic Association Score (RGS) via MAGMA
  trs            Trait Relevance Score (TRS) integrating RAS and RGS
"""

import argparse
import sys

from scrbp.commands import (
    get_sketch,
    get_grn,
    get_merge_grn,
    get_module,
    get_prune,
    get_regulon,
    merge_regulons,
    ras,
    rgs,
    trs,
)


def main():
    parser = argparse.ArgumentParser(
        prog="scRBP",
        description="scRBP — Single-cell RNA Binding Protein regulon inference",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s 0.1.3.11",
    )

    subparsers = parser.add_subparsers(
        title="subcommands",
        dest="subcommand",
        metavar="<command>",
    )

    # Register all subcommands
    get_sketch.register_subcommand(subparsers)
    get_grn.register_subcommand(subparsers)
    get_merge_grn.register_subcommand(subparsers)
    get_module.register_subcommand(subparsers)
    get_prune.register_subcommand(subparsers)
    get_regulon.register_subcommand(subparsers)
    merge_regulons.register_subcommand(subparsers)
    ras.register_subcommand(subparsers)
    rgs.register_subcommand(subparsers)
    trs.register_subcommand(subparsers)

    args = parser.parse_args()

    if args.subcommand is None:
        parser.print_help()
        sys.exit(0)

    args.func(args)


if __name__ == "__main__":
    main()
