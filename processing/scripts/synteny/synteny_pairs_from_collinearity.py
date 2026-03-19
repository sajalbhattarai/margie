#!/usr/bin/env python3
"""Parse MCScanX .collinearity output into a flat syntenic gene-pair table."""

from __future__ import annotations

import argparse
import re
import pandas as pd

BLOCK_RE = re.compile(r"^##\s+Alignment\s+(\d+):")
PAIR_RE = re.compile(r"^\s*\d+-\s*\d+:\s+(\S+)\s+(\S+)\s+(\S+)")


def parse_collinearity(path: str) -> pd.DataFrame:
    block_id = None
    records = []

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")

            m_block = BLOCK_RE.match(line)
            if m_block:
                block_id = int(m_block.group(1))
                continue

            m_pair = PAIR_RE.match(line)
            if m_pair and block_id is not None:
                gene_a = m_pair.group(1)
                gene_b = m_pair.group(2)
                evalue = m_pair.group(3)
                records.append(
                    {
                        "block_id": block_id,
                        "gene_a": gene_a,
                        "gene_b": gene_b,
                        "anchor_evalue": evalue,
                    }
                )

    return pd.DataFrame(records)


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract MCScanX syntenic gene pairs")
    parser.add_argument("--collinearity", required=True, help="Input .collinearity file")
    parser.add_argument("--output", required=True, help="Output TSV for gene pairs")
    args = parser.parse_args()

    df = parse_collinearity(args.collinearity)
    if df.empty:
        df = pd.DataFrame(columns=["block_id", "gene_a", "gene_b", "anchor_evalue"])

    df.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(df)} syntenic gene pairs to {args.output}")


if __name__ == "__main__":
    main()
