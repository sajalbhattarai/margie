#!/usr/bin/env python3
"""Add genome/contig coordinate context to MCScanX gene-pair output."""

from __future__ import annotations

import argparse
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Annotate synteny pairs with loci metadata")
    parser.add_argument("--pairs", required=True, help="Input synteny pairs TSV")
    parser.add_argument("--loci", required=True, help="Gene loci TSV generated during MCScanX input prep")
    parser.add_argument("--output", required=True, help="Output TSV")
    args = parser.parse_args()

    pairs = pd.read_csv(args.pairs, sep="\t")
    loci = pd.read_csv(args.loci, sep="\t")

    if pairs.empty:
        pairs.to_csv(args.output, sep="\t", index=False)
        print(f"No pairs to annotate; wrote empty output to {args.output}")
        return

    left = loci.add_suffix("_a").rename(columns={"gene_id_a": "gene_a"})
    right = loci.add_suffix("_b").rename(columns={"gene_id_b": "gene_b"})

    out = pairs.merge(left, on="gene_a", how="left").merge(right, on="gene_b", how="left")
    out.to_csv(args.output, sep="\t", index=False)
    print(f"Annotated {len(out)} syntenic pairs with loci metadata: {args.output}")


if __name__ == "__main__":
    main()
