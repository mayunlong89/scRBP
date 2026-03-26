#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scRBP getRegulon — Export pruned RBP regulons to GMT format

Converts the motif-enrichment pruning output (from getPrune) into two GMT
files used by downstream steps:

  1. Symbol GMT  — leading-edge target genes identified by gene symbol,
                   used by ``ras`` for Regulon Activity Score (RAS) computation.
  2. Entrez GMT  — target genes mapped to Entrez IDs,
                   used by ``rgs`` for MAGMA gene-set analysis (RGS computation).

Symbol-to-Entrez mapping is resolved from one or more user-provided
reference files (MAGMA gene.loc, HGNC complete set, or NCBI gene_info).
The MAGMA gene.loc file (``--map-custom``) is recommended as it ensures
consistency with the reference panel used for MAGMA gene analysis.

Sets with fewer than ``--min_genes`` targets after mapping are discarded.
"""

import argparse
import os
import gzip
import re
import ast
from collections import defaultdict
from typing import Dict

import pandas as pd


########################
# Parse the leading-edge gene column (handles dict/set/list string formats)
########################
def parse_genes(cell: str):
    if cell is None:
        return []
    s = str(cell).strip()
    if not s:
        return []
    if s.startswith("{") and s.endswith("}"):
        try:
            vals = list(ast.literal_eval(s))
            return [str(x).strip() for x in vals if str(x).strip()]
        except Exception:
            pass
    if "(" in s and ")" in s and "'" in s:
        try:
            items = ast.literal_eval("[" + s + "]")
            out = []
            for it in items:
                if isinstance(it, (list, tuple)) and len(it) >= 1:
                    out.append(str(it[0]).strip())
            if out:
                return out
        except Exception:
            pass
        return [m for m in re.findall(r"'([^']+)'", s) if m.strip()]
    return [g.strip() for g in s.split(",") if g.strip()]


########################
# Symbol -> Entrez mapping loaders
########################
def load_map_from_magma_gene_loc(path: str) -> Dict[str, str]:
    """
    Load a symbol-to-Entrez mapping from a MAGMA NCBI*.gene.loc file
    (whitespace-delimited; column 1 = Entrez ID, column 6 = gene symbol).
    """
    df = pd.read_csv(path, sep=r"\s+", header=None, engine="python", usecols=[0, 5])
    df.columns = ["entrez", "symbol"]
    df = df.dropna()
    df["entrez"] = df["entrez"].astype(str)
    df["symbol"] = df["symbol"].astype(str)
    df = df[(df["symbol"] != "-") & (df["symbol"] != "")]
    m: Dict[str, str] = {}
    for sym, eid in df[["symbol", "entrez"]].itertuples(index=False):
        m.setdefault(sym.upper(), eid)
    return m


def load_map_from_hgnc(path: str) -> Dict[str, str]:
    # HGNC complete set (TSV); lowercase column names for cross-version compatibility
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    df.columns = [c.lower() for c in df.columns]
    if "symbol" not in df.columns:
        raise ValueError(f"{path} does not contain 'symbol' column")
    if "entrez_id" not in df.columns:
        raise ValueError(f"{path} does not contain 'entrez_id' column (please confirm download is hgnc_complete_set.txt)")

    m: Dict[str, str] = {}
    # Primary symbol
    main = df[["symbol", "entrez_id"]].dropna()
    for sym, entrez in main.itertuples(index=False):
        if sym and entrez:
            m[sym.upper()] = str(entrez)

    # Alias and previous symbols (columns may not exist in all HGNC versions; handled gracefully)
    alias_col = "alias_symbol" if "alias_symbol" in df.columns else None
    prev_col  = "prev_symbol"  if "prev_symbol"  in df.columns else None
    if alias_col or prev_col:
        for _, row in df.iterrows():
            entrez = row.get("entrez_id")
            if not isinstance(entrez, str) or not entrez:
                continue
            for col in (alias_col, prev_col):
                if not col:
                    continue
                bucket = row.get(col)
                if not isinstance(bucket, str) or not bucket:
                    continue
                for a in re.split(r"[,\|; ]+", bucket):
                    a = a.strip()
                    if a:
                        m.setdefault(a.upper(), str(entrez))
    return m


def load_map_from_ncbi_gene_info(path: str, tax_id: int = 9606) -> Dict[str, str]:
    # Read NCBI gene_info (do NOT pass comment="#" — it would drop the '#tax_id' header row)
    opener = gzip.open if str(path).endswith(".gz") else open
    df = pd.read_csv(opener(path, "rt"), sep="\t", dtype=str, low_memory=False)
    # Normalise column names: lowercase and strip leading '#'
    df.columns = [c.lower().lstrip("#") for c in df.columns]
    need = {"tax_id", "geneid", "symbol", "synonyms"}
    if not need.issubset(df.columns):
        raise ValueError(f"{path} missing columns: {need - set(df.columns)}")

    df = df[df["tax_id"].astype(str) == str(tax_id)]
    m: Dict[str, str] = {}
    for _, r in df.iterrows():
        eid = str(r["geneid"])
        sym = (r.get("symbol") or "").strip()
        if sym:
            m.setdefault(sym.upper(), eid)
        syn = (r.get("synonyms") or "")
        if syn and syn != "-":
            for a in syn.split("|"):
                a = a.strip()
                if a:
                    m.setdefault(a.upper(), eid)
    return m


def compose_symbol2entrez_map(args) -> Dict[str, str]:
    # Merge mapping sources in order of increasing priority: NCBI < HGNC < custom (MAGMA gene.loc)
    m: Dict[str, str] = {}
    if args.map_ncbi:
        m.update(load_map_from_ncbi_gene_info(args.map_ncbi, tax_id=args.taxid))
    if args.map_hgnc:
        m.update(load_map_from_hgnc(args.map_hgnc))
    if args.map_custom:
        m.update(load_map_from_magma_gene_loc(args.map_custom))
    if not m:
        raise SystemExit("Please provide --map-custom (NCBI*.gene.loc) or --map-hgnc / --map-ncbi (at least one).")
    return m


from scrbp.utils.io import write_gmt  # shared with rgs


def register_subcommand(subparsers):
    """Register this subcommand with the main parser."""
    ap = subparsers.add_parser(
        'getRegulon',
        help='Convert context scores to GMT regulon format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )
    ap.add_argument("--input", required=True, help="ctx/prune CSV file")
    ap.add_argument("--out-symbol", required=True, help="output GMT (gene symbols)")
    ap.add_argument("--out-entrez", required=True, help="output GMT (Entrez IDs)")
    ap.add_argument("--rbp_col", default="RBP", help="column name for RBP (default: RBP)")
    ap.add_argument("--genes_col", default=None,
                    help="column name for LeadingEdge genes (auto-detect if not set)")
    ap.add_argument("--min_genes", type=int, default=1, help="min targets to keep an RBP line")

    # 映射来源（仅保留这三种；推荐 --map-custom=NCBI38.gene.loc 或 NCBI37.3.gene.loc）
    ap.add_argument("--map-custom", default=None,
                    help="MAGMA gene.loc file (e.g., NCBI38.gene.loc). Column 1 = Entrez ID, Column 6 = gene symbol.")
    ap.add_argument("--map-hgnc", default=None,
                    help="HGNC complete set TSV (must contain 'symbol' and 'entrez_id' columns).")
    ap.add_argument("--map-ncbi", default=None,
                    help="NCBI Homo_sapiens.gene_info or Homo_sapiens.gene_info.gz.")
    ap.add_argument("--taxid", type=int, default=9606, help="NCBI taxonomy ID (default: 9606 = human).")

    # Additional options
    ap.add_argument("--drop-unmapped-genes", action="store_true",
                    help="drop genes without Entrez mapping in the Entrez GMT (default: keep them out)")
    ap.add_argument("--drop-empty-sets", action="store_true",
                    help="drop sets that become empty (< min_genes) after Entrez mapping")
    ap.add_argument("--sanitize-rbp", action="store_true",
                    help="replace whitespace in RBP name with '_'")
    ap.set_defaults(func=main)
    return ap


def main(args):
    """Main execution function."""
    df = pd.read_csv(args.input)

    # Auto-detect the leading-edge genes column
    genes_col = args.genes_col
    if genes_col is None:
        for cand in ["LeadingEdge-Genes", "Leading_Edge_Genes", "LeadingEdge_Genes"]:
            if cand in df.columns:
                genes_col = cand
                break
    if genes_col is None:
        raise ValueError("Cannot find LeadingEdge gene column (tried: LeadingEdge-Genes / Leading_Edge_Genes / LeadingEdge_Genes)")
    if args.rbp_col not in df.columns:
        raise ValueError(f"Cannot find RBP column: {args.rbp_col}")

    # Aggregate: RBP -> set of leading-edge target genes (symbol)
    rbp2genes_symbol = defaultdict(set)
    for _, row in df.iterrows():
        rbp = str(row[args.rbp_col]).strip()
        if args.sanitize_rbp:
            rbp = re.sub(r"\s+", "_", rbp)
        genes = parse_genes(row[genes_col])
        for g in genes:
            if g:
                rbp2genes_symbol[rbp].add(g)

    # Write symbol GMT (filtered by min_genes)
    rbp2genes_symbol_filtered = {k: {g for g in v if g} for k, v in rbp2genes_symbol.items()}
    rbp2genes_symbol_filtered = {k: v for k, v in rbp2genes_symbol_filtered.items() if len(v) >= args.min_genes}
    kept_sym = write_gmt(args.out_symbol, rbp2genes_symbol_filtered)
    print(f"[SYMBOL] Wrote {kept_sym} regulons to {args.out_symbol}")

    # Build symbol-to-Entrez mapping (recommended: --map-custom=NCBI*.gene.loc)
    sym2ent = compose_symbol2entrez_map(args)

    # Convert symbol regulons to Entrez regulons
    total_symbols = 0
    mapped_symbols = 0
    rbp2genes_entrez = {}
    for rbp, genes in rbp2genes_symbol.items():
        entrez_set = set()
        for g in genes:
            total_symbols += 1
            e = sym2ent.get(g.upper())
            if e:
                mapped_symbols += 1
                entrez_set.add(str(e))
            elif not args.drop_unmapped_genes:
                # Unmapped symbols are excluded from the Entrez set (original symbol not retained)
                pass
        if len(entrez_set) >= args.min_genes or (not args.drop_empty_sets and len(entrez_set) > 0):
            rbp2genes_entrez[rbp] = entrez_set

    cov = 0.0 if total_symbols == 0 else mapped_symbols / total_symbols
    kept_ent = write_gmt(args.out_entrez, rbp2genes_entrez)
    print(f"[ENTREZ] Wrote {kept_ent} regulons to {args.out_entrez}")
    print(f"[MAPPING] Symbol->Entrez coverage: {mapped_symbols}/{total_symbols} = {cov:.1%}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    subp = p.add_subparsers()
    register_subcommand(subp)
    args = p.parse_args()
    main(args)
