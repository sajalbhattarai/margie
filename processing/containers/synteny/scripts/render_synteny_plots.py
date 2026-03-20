#!/usr/bin/env python3
"""Generate simple diagrammatic summaries from synteny pairs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Render synteny summary plots")
    parser.add_argument("--pairs", required=True, help="Input synteny pairs TSV")
    parser.add_argument("--output-dir", required=True, help="Directory for plots")
    args = parser.parse_args()

    pairs = pd.read_csv(args.pairs, sep="\t")
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if pairs.empty or "block_id" not in pairs.columns:
        # Create a placeholder image to indicate no syntenic blocks were detected.
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No syntenic pairs detected", ha="center", va="center", fontsize=14)
        ax.axis("off")
        fig.tight_layout()
        fig.savefig(out_dir / "synteny_block_sizes.png", dpi=150)
        plt.close(fig)
        print(f"No synteny pairs available; wrote placeholder plot to {out_dir / 'synteny_block_sizes.png'}")
        return

    counts = pairs.groupby("block_id", as_index=False).size().sort_values("size", ascending=False)
    top = counts.head(30)

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(top["block_id"].astype(str), top["size"], color="#2a6f97")
    ax.set_xlabel("Synteny block ID")
    ax.set_ylabel("Anchored gene pairs")
    ax.set_title("Top synteny blocks by anchored pair count")
    ax.tick_params(axis="x", rotation=90)
    fig.tight_layout()
    fig.savefig(out_dir / "synteny_block_sizes.png", dpi=150)
    plt.close(fig)

    counts.to_csv(out_dir / "synteny_block_sizes.tsv", sep="\t", index=False)
    print(f"Wrote plot: {out_dir / 'synteny_block_sizes.png'}")
    print(f"Wrote table: {out_dir / 'synteny_block_sizes.tsv'}")


if __name__ == "__main__":
    main()
