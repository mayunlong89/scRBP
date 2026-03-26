#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scRBP mergeRegulons - Merge region-specific GMT regulon files into a combined GMT file"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import List, Optional, Sequence, Tuple, Dict, Any

import pandas as pd


def setup_logger(verbose: bool = True) -> None:
    """
    Configure logging.

    Parameters
    ----------
    verbose : bool, default=True
        If True, use INFO logging level. Otherwise use WARNING.
    """
    level = logging.INFO if verbose else logging.WARNING
    logging.basicConfig(
        level=level,
        format="[%(asctime)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )


def extract_region_name(dirname: str, region_order: Sequence[str]) -> str:
    """
    Extract transcript region name from directory basename.

    Priority
    --------
    1. If basename ends with one of the known region names, return that region.
    2. Otherwise return the last underscore-delimited token.

    Parameters
    ----------
    dirname : str
        Directory path or basename.
    region_order : sequence of str
        Known region names.

    Returns
    -------
    str
        Extracted region name.
    """
    base = Path(dirname).name

    for region in region_order:
        if base.endswith(f"_{region}"):
            return region

    if "_" in base:
        return base.rsplit("_", 1)[-1]

    return base


def sort_region_dirs(region_dirs: Sequence[Path], region_order: Sequence[str]) -> List[Path]:
    """
    Sort region directories according to user-defined region order.

    Unknown regions are placed at the end.

    Parameters
    ----------
    region_dirs : sequence of pathlib.Path
        Region-specific directories.
    region_order : sequence of str
        Preferred region order.

    Returns
    -------
    list of pathlib.Path
        Sorted directories.
    """
    order_map = {region: i for i, region in enumerate(region_order)}

    def _key(path: Path) -> Tuple[int, str]:
        region = extract_region_name(str(path), region_order)
        return (order_map.get(region, 999), str(path))

    return sorted(region_dirs, key=_key)


def deduplicate_preserve_order(lines: Sequence[str]) -> List[str]:
    """
    Deduplicate lines while preserving original order.

    Parameters
    ----------
    lines : sequence of str
        Input lines.

    Returns
    -------
    list of str
        Deduplicated lines.
    """
    seen = set()
    output = []
    for line in lines:
        if line not in seen:
            seen.add(line)
            output.append(line)
    return output


def read_gmt_lines(
    infile: Path,
    region: str,
    append_region_to_setname: bool = True,
    min_fields: int = 3,
) -> Tuple[List[str], int, int]:
    """
    Read GMT file and optionally append region suffix to set name.

    Expected GMT format
    -------------------
    set_name \\t description \\t gene1 \\t gene2 ...

    Parameters
    ----------
    infile : pathlib.Path
        Input GMT file.
    region : str
        Region suffix to append.
    append_region_to_setname : bool, default=True
        Whether to append '_{region}' to the set name.
    min_fields : int, default=3
        Minimum number of fields required for a valid GMT line.

    Returns
    -------
    tuple
        (processed_lines, n_total_lines, n_valid_lines)
    """
    output_lines = []
    n_total = 0
    n_valid = 0

    with infile.open("r", encoding="utf-8") as f:
        for line in f:
            n_total += 1
            line = line.rstrip("\n")
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) < min_fields:
                continue

            if append_region_to_setname:
                fields[0] = f"{fields[0]}_{region}"

            output_lines.append("\t".join(fields))
            n_valid += 1

    return output_lines, n_total, n_valid


def merge_one_gmt_type(
    parent_dir: Path,
    region_dirs: Sequence[Path],
    input_name: str,
    output_name: str,
    region_order: Sequence[str],
    append_region_to_setname: bool = True,
    dedup_lines: bool = True,
    min_fields: int = 3,
    overwrite: bool = True,
) -> Dict[str, Any]:
    """
    Merge one GMT file type across all region-specific directories.

    Parameters
    ----------
    parent_dir : pathlib.Path
        Parent analysis directory.
    region_dirs : sequence of pathlib.Path
        Region-specific directories.
    input_name : str
        Input GMT filename inside each region directory.
    output_name : str
        Output merged GMT filename written to parent_dir.
    region_order : sequence of str
        Preferred region ordering.
    append_region_to_setname : bool, default=True
        Whether to append region suffix to set names.
    dedup_lines : bool, default=True
        Whether to deduplicate identical merged lines.
    min_fields : int, default=3
        Minimum number of fields for valid GMT lines.
    overwrite : bool, default=True
        Whether to overwrite existing output file.

    Returns
    -------
    dict
        Summary dictionary for this merge task.
    """
    all_lines = []
    region_stats = []

    for rdir in region_dirs:
        infile = rdir / input_name
        region = extract_region_name(str(rdir), region_order)

        if not infile.exists():
            logging.warning(f"Missing file: {infile}")
            region_stats.append({
                "parent_dir": str(parent_dir),
                "region_dir": str(rdir),
                "region": region,
                "input_name": input_name,
                "output_name": output_name,
                "input_file_exists": False,
                "n_total_lines": 0,
                "n_valid_lines": 0,
                "n_output_lines_added": 0,
            })
            continue

        logging.info(f"Merging {infile.name} from region: {region}")
        lines, n_total, n_valid = read_gmt_lines(
            infile=infile,
            region=region,
            append_region_to_setname=append_region_to_setname,
            min_fields=min_fields,
        )
        all_lines.extend(lines)

        region_stats.append({
            "parent_dir": str(parent_dir),
            "region_dir": str(rdir),
            "region": region,
            "input_name": input_name,
            "output_name": output_name,
            "input_file_exists": True,
            "n_total_lines": n_total,
            "n_valid_lines": n_valid,
            "n_output_lines_added": len(lines),
        })

    n_before = len(all_lines)

    if dedup_lines:
        all_lines = deduplicate_preserve_order(all_lines)

    n_after = len(all_lines)
    n_removed = n_before - n_after

    out_file = parent_dir / output_name

    if out_file.exists() and not overwrite:
        raise FileExistsError(f"Output file exists and overwrite=False: {out_file}")

    if len(all_lines) == 0:
        logging.warning(f"No valid GMT lines found for {input_name} under {parent_dir}")
        return {
            "parent_dir": str(parent_dir),
            "input_name": input_name,
            "output_name": output_name,
            "output_file": str(out_file),
            "output_written": False,
            "n_regions_scanned": len(region_dirs),
            "n_regions_with_existing_input": sum(x["input_file_exists"] for x in region_stats),
            "n_merged_before_dedup": n_before,
            "n_merged_after_dedup": n_after,
            "n_removed_by_dedup": n_removed,
            "region_stats": region_stats,
        }

    with out_file.open("w", encoding="utf-8") as f:
        for line in all_lines:
            f.write(line + "\n")

    logging.info(f"OK: wrote {out_file}")

    return {
        "parent_dir": str(parent_dir),
        "input_name": input_name,
        "output_name": output_name,
        "output_file": str(out_file),
        "output_written": True,
        "n_regions_scanned": len(region_dirs),
        "n_regions_with_existing_input": sum(x["input_file_exists"] for x in region_stats),
        "n_merged_before_dedup": n_before,
        "n_merged_after_dedup": n_after,
        "n_removed_by_dedup": n_removed,
        "region_stats": region_stats,
    }


def scRBP_mergeRegulons(
    base_dir: str,
    input_name: str,
    output_name: str,
    recursive: bool = False,
    tissue_glob: str = "z_GRNBoost2_*_30times",
    region_glob: str = "Results_final_*_RBP_top1500_*",
    append_region_to_setname: bool = True,
    dedup_lines: bool = True,
    region_order: Optional[Sequence[str]] = None,
    min_fields: int = 3,
    overwrite: bool = True,
    return_region_level_summary: bool = False,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Merge region-specific regulon GMT files into one combined GMT file.

    Parameters
    ----------
    base_dir : str
        Base directory to process.
    input_name : str
        Input GMT filename expected inside each region-specific subdirectory.
    output_name : str
        Output merged GMT filename written to each parent directory.
    recursive : bool, default=False
        If True, process multiple parent directories under base_dir matching tissue_glob.
        If False, process base_dir itself.
    tissue_glob : str, default="z_GRNBoost2_*_30times"
        Glob pattern used to identify parent directories when recursive=True.
    region_glob : str, default="Results_final_*_RBP_top1500_*"
        Glob pattern used to identify region-specific subdirectories.
    append_region_to_setname : bool, default=True
        Whether to append region suffix to set names.
    dedup_lines : bool, default=True
        Whether to deduplicate identical merged lines.
    region_order : sequence of str, optional
        Preferred order of regions.
    min_fields : int, default=3
        Minimum number of tab-delimited fields required for valid GMT lines.
    overwrite : bool, default=True
        Whether to overwrite existing output files.
    return_region_level_summary : bool, default=False
        If True, return region-level summary. Otherwise return task-level summary.
    verbose : bool, default=True
        Whether to print INFO logs.

    Returns
    -------
    pandas.DataFrame
        Summary dataframe.
    """
    setup_logger(verbose=verbose)

    if region_order is None:
        region_order = ["3UTR", "5UTR", "CDS", "Introns"]

    base_path = Path(base_dir)
    if not base_path.exists():
        raise FileNotFoundError(f"Base directory does not exist: {base_path}")

    if recursive:
        parent_dirs = [p for p in base_path.glob(tissue_glob) if p.is_dir()]
    else:
        parent_dirs = [base_path]

    if len(parent_dirs) == 0:
        raise FileNotFoundError("No target parent directories found.")

    task_records = []
    region_records = []

    logging.info(f"BASE_DIR = {base_path}")
    logging.info(f"INPUT = {input_name}")
    logging.info(f"OUTPUT = {output_name}")
    logging.info(f"RECURSIVE = {int(recursive)}")
    logging.info(f"TISSUE_GLOB = {tissue_glob}")
    logging.info(f"REGION_GLOB = {region_glob}")
    logging.info(f"REGION_ORDER = {list(region_order)}")

    for parent_dir in parent_dirs:
        logging.info("=" * 60)
        logging.info(f"Processing: {parent_dir}")
        logging.info("=" * 60)

        region_dirs = [p for p in parent_dir.glob(region_glob) if p.is_dir()]
        if len(region_dirs) == 0:
            logging.warning(f"No region dirs matched: {parent_dir / region_glob}")
            continue

        region_dirs = sort_region_dirs(region_dirs, region_order)

        logging.info("Found region dirs:")
        for d in region_dirs:
            logging.info(f"  - {d}")

        summary = merge_one_gmt_type(
            parent_dir=parent_dir,
            region_dirs=region_dirs,
            input_name=input_name,
            output_name=output_name,
            region_order=region_order,
            append_region_to_setname=append_region_to_setname,
            dedup_lines=dedup_lines,
            min_fields=min_fields,
            overwrite=overwrite,
        )

        task_records.append({
            "parent_dir": summary["parent_dir"],
            "input_name": summary["input_name"],
            "output_name": summary["output_name"],
            "output_file": summary["output_file"],
            "output_written": summary["output_written"],
            "n_regions_scanned": summary["n_regions_scanned"],
            "n_regions_with_existing_input": summary["n_regions_with_existing_input"],
            "n_merged_before_dedup": summary["n_merged_before_dedup"],
            "n_merged_after_dedup": summary["n_merged_after_dedup"],
            "n_removed_by_dedup": summary["n_removed_by_dedup"],
        })

        region_records.extend(summary["region_stats"])

    if return_region_level_summary:
        return pd.DataFrame(region_records)

    return pd.DataFrame(task_records)


def register_subcommand(subparsers):
    """Register this subcommand with the main parser."""
    parser = subparsers.add_parser(
        'mergeRegulons',
        help='Merge region-specific GMT regulon files into a combined GMT file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )

    parser.add_argument(
        "--base_dir",
        required=True,
        help="Base directory to process."
    )
    parser.add_argument(
        "--input",
        required=True,
        dest="input_name",
        help="Input GMT filename expected inside each region directory."
    )
    parser.add_argument(
        "--output",
        required=True,
        dest="output_name",
        help="Output merged GMT filename."
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Recursively process multiple parent directories under base_dir."
    )
    parser.add_argument(
        "--tissue_glob",
        default="z_GRNBoost2_*_30times",
        help="Glob pattern for parent dirs when --recursive is set."
    )
    parser.add_argument(
        "--region_glob",
        default="Results_final_*_RBP_top1500_*",
        help="Glob pattern for region-specific subdirectories."
    )
    parser.add_argument(
        "--append_region_to_setname",
        action="store_true",
        help="Append region suffix to set names."
    )
    parser.add_argument(
        "--dedup_lines",
        action="store_true",
        help="Deduplicate identical merged GMT lines."
    )
    parser.add_argument(
        "--region_order",
        nargs="+",
        default=["3UTR", "5UTR", "CDS", "Introns"],
        help="Preferred region order."
    )
    parser.add_argument(
        "--min_fields",
        type=int,
        default=3,
        help="Minimum number of fields required for a valid GMT line."
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files."
    )
    parser.add_argument(
        "--return_region_level_summary",
        action="store_true",
        help="Return region-level summary."
    )
    parser.add_argument(
        "--summary_out",
        default=None,
        help="Optional output TSV for summary table."
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress INFO logs."
    )

    parser.set_defaults(func=main)
    return parser


def main(args) -> None:
    """Main execution function."""
    summary_df = scRBP_mergeRegulons(
        base_dir=args.base_dir,
        input_name=args.input_name,
        output_name=args.output_name,
        recursive=args.recursive,
        tissue_glob=args.tissue_glob,
        region_glob=args.region_glob,
        append_region_to_setname=args.append_region_to_setname,
        dedup_lines=args.dedup_lines,
        region_order=args.region_order,
        min_fields=args.min_fields,
        overwrite=args.overwrite,
        return_region_level_summary=args.return_region_level_summary,
        verbose=not args.quiet,
    )

    if args.summary_out is not None:
        summary_path = Path(args.summary_out)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(summary_path, sep="\t", index=False)
        logging.info(f"Summary saved to: {summary_path}")
    else:
        print(summary_df.to_csv(sep="\t", index=False))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    subp = p.add_subparsers()
    register_subcommand(subp)
    args = p.parse_args()
    main(args)
