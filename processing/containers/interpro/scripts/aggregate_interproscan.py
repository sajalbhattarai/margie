#!/usr/bin/env python3
"""Aggregate InterProScan TSV output to Margie interpro.tsv schema."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate InterProScan TSV per feature")
    parser.add_argument("--organism", required=True)
    parser.add_argument("--interproscan", required=True, help="Path to InterProScan TSV output")
    parser.add_argument("--output", required=True, help="Path to aggregated interpro.tsv")
    args = parser.parse_args()

    source = Path(args.interproscan)
    out = Path(args.output)

    if not source.exists() or source.stat().st_size == 0:
        raise FileNotFoundError(f"InterProScan TSV not found or empty: {source}")

    ipr_by_feature: dict[str, set[str]] = defaultdict(set)

    with source.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            feature_id = (row[0] if len(row) > 0 else "").strip()
            ipr_acc = (row[11] if len(row) > 11 else "").strip()
            if feature_id and ipr_acc and ipr_acc != "-":
                ipr_by_feature[feature_id].add(ipr_acc)

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "organism",
                "feature_id",
                "INTERPRO_ipr_pfam_ids",
                "INTERPRO_ipr_tigrfam_ids",
                "INTERPRO_ipr_overlap_count",
                "INTERPRO_ipr_support_flags",
                "INTERPRO_pfam_signatures",
                "INTERPRO_tigrfam_signatures",
                "INTERPRO_pfam_descriptions",
                "INTERPRO_tigrfam_descriptions",
            ]
        )

        for feature_id in sorted(ipr_by_feature):
            iprs = sorted(ipr_by_feature[feature_id])
            writer.writerow(
                [
                    args.organism,
                    feature_id,
                    ";".join(iprs),
                    "",
                    0,
                    "IPR_INTERPROSCAN",
                    "",
                    "",
                    "",
                    "",
                ]
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
