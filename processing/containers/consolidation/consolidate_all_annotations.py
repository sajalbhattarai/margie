#!/usr/bin/env python3
"""Merge per-tool TSV outputs into a single consolidated table.

This script is intentionally tolerant of heterogeneous tool outputs. It unions
columns across discovered TSV files and appends two provenance columns:
`source_tool` and `source_file`.
"""

import argparse
from pathlib import Path
import sys
from typing import List

import pandas as pd


def discover_tsvs(input_root: Path) -> List[Path]:
    files: List[Path] = []
    for path in sorted(input_root.rglob("*.tsv")):
        if path.name == "consolidated_all_genomes.tsv":
            continue
        files.append(path)
    return files


def source_tool_from_path(path: Path, input_root: Path) -> str:
    rel = path.relative_to(input_root)
    if len(rel.parts) > 1:
        return rel.parts[0]
    return "unknown"


def merge_tables(files: List[Path], input_root: Path) -> pd.DataFrame:
    frames: List[pd.DataFrame] = []
    for tsv in files:
        try:
            df = pd.read_csv(tsv, sep="\t", dtype=str)
        except Exception:
            # Skip malformed/empty files to keep consolidation non-blocking.
            continue
        if df.empty:
            continue

        df["source_tool"] = source_tool_from_path(tsv, input_root)
        df["source_file"] = str(tsv)
        frames.append(df)

    if not frames:
        return pd.DataFrame(columns=["source_tool", "source_file"])

    return pd.concat(frames, ignore_index=True, sort=False)


def main() -> int:
    parser = argparse.ArgumentParser(description="Consolidate annotation TSV files")
    parser.add_argument(
        "--input-root",
        default="/container/output",
        help="Root directory containing per-tool output folders (default: /container/output)",
    )
    parser.add_argument(
        "--output",
        default="/container/output/consolidated_all_genomes.tsv",
        help="Output TSV path (default: /container/output/consolidated_all_genomes.tsv)",
    )
    args = parser.parse_args()

    input_root = Path(args.input_root)
    output_path = Path(args.output)

    if not input_root.exists():
        print(f"[ERROR] Input root not found: {input_root}", file=sys.stderr)
        return 1

    files = discover_tsvs(input_root)
    print(f"[INFO] Found {len(files)} TSV files under {input_root}")

    merged = merge_tables(files, input_root)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_path, sep="\t", index=False)
    print(f"[INFO] Wrote {len(merged)} rows to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
