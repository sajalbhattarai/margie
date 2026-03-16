#!/usr/bin/env python3
"""Fallback COG classifier that tolerates unmapped CDD IDs.

This script runs rpsblast against Cog_LE and builds COG outputs while skipping
CDD IDs that do not map to a COG accession in cddid.tbl(.gz).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

from cogclassifier.cog import CogDefinitionRecord, CogFuncCategoryRecord
import cogclassifier


def load_cdd_to_cog_map(cddid_path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    opener = gzip.open if cddid_path.suffix == ".gz" else open
    with opener(cddid_path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            cdd_id, acc_id = row[0], row[1]
            if acc_id.startswith("COG"):
                mapping[cdd_id] = acc_id
    return mapping


def run_rpsblast(input_faa: Path, db_prefix: Path, out_tsv: Path, threads: int, evalue: float) -> None:
    cmd = [
        "rpsblast",
        "-query",
        str(input_faa),
        "-db",
        str(db_prefix),
        "-outfmt",
        "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-out",
        str(out_tsv),
        "-evalue",
        str(evalue),
        "-num_threads",
        str(threads),
        "-mt_mode",
        "1",
    ]
    subprocess.run(cmd, check=True)


def load_best_hits(rps_tsv: Path) -> Dict[str, Tuple[str, float, float]]:
    best: Dict[str, Tuple[str, float, float]] = {}
    with open(rps_tsv, "r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            query = parts[0]
            saccver = parts[1].replace("CDD:", "")
            pident = float(parts[2])
            evalue = float(parts[10])
            bitscore = float(parts[11])
            prev = best.get(query)
            if prev is None or evalue < prev[1] or (evalue == prev[1] and bitscore > prev[2]):
                best[query] = (saccver, evalue, pident)
    return best


def count_query_sequences(input_faa: Path) -> int:
    count = 0
    with open(input_faa, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def main() -> int:
    parser = argparse.ArgumentParser(description="Safe COG fallback classifier")
    parser.add_argument("--input", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--db", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--evalue", type=float, default=0.01)
    args = parser.parse_args()

    input_faa = Path(args.input)
    outdir = Path(args.outdir)
    db_dir = Path(args.db)

    outdir.mkdir(parents=True, exist_ok=True)

    db_prefix = db_dir / "Cog_LE" / "Cog"
    cddid_tbl = db_dir / "cddid.tbl.gz"
    if not cddid_tbl.exists():
        cddid_tbl = db_dir / "cddid.tbl"

    if not input_faa.exists():
        raise FileNotFoundError(f"Input FASTA not found: {input_faa}")
    if not (db_prefix.with_suffix(".pin").exists() or db_prefix.with_suffix(".phr").exists()):
        raise FileNotFoundError(f"RPS-BLAST DB not found: {db_prefix}")
    if not cddid_tbl.exists():
        raise FileNotFoundError(f"CDD mapping file not found: {cddid_tbl}")

    resource_dir = Path(cogclassifier.__file__).resolve().parent / "resources"
    cog_fc = CogFuncCategoryRecord(resource_dir / "cog_func_category.tsv")
    cog_def = CogDefinitionRecord(resource_dir / "cog_definition.tsv")
    cdd_to_cog = load_cdd_to_cog_map(cddid_tbl)

    rps_out = outdir / "rpsblast.tsv"
    run_rpsblast(input_faa, db_prefix, rps_out, args.threads, args.evalue)

    best_hits = load_best_hits(rps_out)

    classify_rows: List[Tuple[str, str, str, str, str, str, str, str, str]] = []
    letter_counts: Counter[str] = Counter()

    for query_id, (cdd_id, evalue, pident) in best_hits.items():
        cog_id = cdd_to_cog.get(cdd_id)
        if not cog_id:
            continue
        cog_info = cog_def.get(cog_id)
        if cog_info is None:
            continue
        letter = cog_info.one_letter
        try:
            fc = cog_fc.get(letter)
            desc = fc.desc
        except Exception:
            desc = ""
        classify_rows.append(
            (
                query_id,
                cog_id,
                cdd_id,
                f"{evalue:.2e}" if evalue > 0 else "0.0",
                f"{pident:.2f}",
                cog_info.gene_name,
                cog_info.cog_name,
                letter,
                desc,
            )
        )
        letter_counts[letter] += 1

    classify_tsv = outdir / "cog_classify.tsv"
    with open(classify_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "QUERY_ID",
                "COG_ID",
                "CDD_ID",
                "EVALUE",
                "IDENTITY",
                "GENE_NAME",
                "COG_NAME",
                "COG_LETTER",
                "COG_DESCRIPTION",
            ]
        )
        writer.writerows(classify_rows)

    count_tsv = outdir / "cog_count.tsv"
    with open(count_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["LETTER", "COUNT", "GROUP", "COLOR", "DESCRIPTION"])
        for fc in cog_fc.get_all():
            writer.writerow([fc.letter, letter_counts.get(fc.letter, 0), fc.group, fc.color, fc.desc])

    query_count = count_query_sequences(input_faa)
    classify_count = len(classify_rows)
    ratio = (100.0 * classify_count / query_count) if query_count else 0.0
    print(f"[safe] Query sequences: {query_count}")
    print(f"[safe] Classified sequences: {classify_count} ({ratio:.2f}%)")
    print(f"[safe] Wrote: {rps_out}")
    print(f"[safe] Wrote: {classify_tsv}")
    print(f"[safe] Wrote: {count_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
