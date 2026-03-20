#!/usr/bin/env python3
"""Filter BLAST/DIAMOND tabular hits to make MCScanX tractable on large datasets.

Input format expected:
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

COLS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


def load_gene_to_genome(loci_tsv: Path) -> dict[str, str]:
    loci = pd.read_csv(loci_tsv, sep="\t", dtype=str, usecols=["gene_id", "genome"])
    loci = loci.dropna(subset=["gene_id", "genome"]).drop_duplicates(subset=["gene_id"], keep="first")
    return dict(zip(loci["gene_id"], loci["genome"]))


def main() -> None:
    parser = argparse.ArgumentParser(description="Filter MCScanX BLAST input for large runs")
    parser.add_argument("--input", required=True, help="Raw blast/diamond outfmt6 file")
    parser.add_argument("--loci", required=True, help="mcscanx_all.gene_loci.tsv")
    parser.add_argument("--output", required=True, help="Filtered output path")
    parser.add_argument("--max-hits-per-query", type=int, default=30, help="Keep top N hits per qseqid by bitscore")
    parser.add_argument("--intergenome-only", action="store_true", help="Drop same-genome hits")
    parser.add_argument("--min-identity", type=float, default=25.0, help="Minimum percent identity")
    parser.add_argument("--min-aln-len", type=int, default=40, help="Minimum alignment length")
    args = parser.parse_args()

    in_path = Path(args.input)
    loci_path = Path(args.loci)
    out_path = Path(args.output)

    if not in_path.exists():
        raise FileNotFoundError(f"Missing input file: {in_path}")
    if not loci_path.exists():
        raise FileNotFoundError(f"Missing loci file: {loci_path}")

    g2g = load_gene_to_genome(loci_path)

    df = pd.read_csv(in_path, sep="\t", header=None, names=COLS, dtype={"qseqid": str, "sseqid": str})
    if df.empty:
        out_path.write_text("", encoding="utf-8")
        print(f"Wrote empty filtered file: {out_path}")
        return

    df = df[df["qseqid"] != df["sseqid"]]
    df = df[pd.to_numeric(df["pident"], errors="coerce") >= float(args.min_identity)]
    df = df[pd.to_numeric(df["length"], errors="coerce") >= int(args.min_aln_len)]

    if args.intergenome_only:
        qg = df["qseqid"].map(g2g)
        sg = df["sseqid"].map(g2g)
        df = df[qg.notna() & sg.notna() & (qg != sg)]

    df["bitscore_num"] = pd.to_numeric(df["bitscore"], errors="coerce").fillna(-1)
    df["evalue_num"] = pd.to_numeric(df["evalue"], errors="coerce").fillna(1e9)

    df = df.sort_values(["qseqid", "bitscore_num", "evalue_num"], ascending=[True, False, True])
    topn = max(1, int(args.max_hits_per_query))
    df = df.groupby("qseqid", as_index=False, sort=False).head(topn)

    df = df[COLS]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False, header=False)
    print(f"Filtered hits written: {out_path}")
    print(f"Rows retained: {len(df)}")


if __name__ == "__main__":
    main()
