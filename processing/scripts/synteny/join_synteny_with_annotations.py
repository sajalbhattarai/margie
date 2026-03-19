#!/usr/bin/env python3
"""Join MCScanX syntenic gene pairs with consolidated annotation tables."""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


CONSOLIDATED_PATTERNS = [
    "*/consolidated_annotations.tsv",
    "*/consolidated_*.tsv",
]


def load_consolidated_tables(consolidated_dir: str) -> pd.DataFrame:
    root = Path(consolidated_dir)

    files = []
    for pattern in CONSOLIDATED_PATTERNS:
        files.extend(root.glob(pattern))
    files = sorted(set(files))

    if not files:
        raise FileNotFoundError(f"No consolidated TSV files found under {consolidated_dir}")

    frames = []
    for fpath in files:
        try:
            df = pd.read_csv(fpath, sep="\t")
            if "feature_id" not in df.columns:
                continue
            df["genome_name"] = fpath.parent.name
            frames.append(df)
        except Exception:
            continue

    if not frames:
        raise RuntimeError("No readable consolidated tables with feature_id were found")

    merged = pd.concat(frames, ignore_index=True)
    merged = merged.drop_duplicates(subset=["feature_id"], keep="first")
    return merged


def main() -> None:
    parser = argparse.ArgumentParser(description="Join syntenic pairs with consolidated annotations")
    parser.add_argument("--pairs", required=True, help="Synteny pair TSV from synteny_pairs_from_collinearity.py")
    parser.add_argument("--consolidated-dir", required=True, help="Directory with per-genome consolidated outputs")
    parser.add_argument("--output", required=True, help="Output TSV")
    args = parser.parse_args()

    pairs = pd.read_csv(args.pairs, sep="\t")
    if pairs.empty:
        pairs.to_csv(args.output, sep="\t", index=False)
        print(f"No pairs to join; wrote empty output to {args.output}")
        return

    annotations = load_consolidated_tables(args.consolidated_dir)

    left = annotations.add_suffix("_a").rename(columns={"feature_id_a": "gene_a"})
    right = annotations.add_suffix("_b").rename(columns={"feature_id_b": "gene_b"})

    out = pairs.merge(left, on="gene_a", how="left").merge(right, on="gene_b", how="left")

    col_a = "consensus_annotation_a"
    col_b = "consensus_annotation_b"
    if col_a in out.columns and col_b in out.columns:
        out["consensus_match"] = out[col_a].fillna("") == out[col_b].fillna("")

    out.to_csv(args.output, sep="\t", index=False)
    print(f"Joined {len(out)} syntenic pairs to annotations: {args.output}")


if __name__ == "__main__":
    main()
