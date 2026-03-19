#!/usr/bin/env python3
"""Apply concordance scoring to one consolidated annotation table."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import pandas as pd

from concordance_common import (
    build_signature,
    calculate_concordance_score,
    extract_gene_set,
    get_operon_info,
)

TOOL_COLUMN_MAP = {
    "RAST": ["feature_id"],
    "COG": ["COG_id"],
    "PFAM": ["PFAM_id", "PFAM_accession"],
    "TIGRFAM": ["TIGRFAM_id"],
    "TIGRFAM_GENPROP": ["TIGRFAM_genprop_ids"],
    "INTERPRO": ["INTERPRO_ipr_support_flags"],
    "KEGG": ["KEGG_KO_IDs"],
    "EGGNOG": ["EGGNOG_hit_protein", "EGGNOG_og"],
    "UNIPROT": ["UNIPROT_accession"],
    "TCDB": ["TCDB_id"],
    "MEROPS": ["MEROPS_id"],
    "DBCAN": ["DBCAN_families"],
}


def build_operon_coregulation_scores(df: pd.DataFrame, concordance_scores: list[int]) -> list[int]:
    columns = list(df.columns)
    operon_id_col = None
    for candidate in ["OPERON_operon_id", "OPERON_id", "operon_id", "Operon_ID"]:
        if candidate in columns:
            operon_id_col = candidate
            break

    if operon_id_col is None:
        return [0] * len(df)

    coreg_scores = [0] * len(df)

    # Build index groups by operon id
    operon_to_indexes = {}
    for idx, value in enumerate(df[operon_id_col].tolist()):
        if pd.isna(value):
            continue
        operon_id = str(value).strip()
        if not operon_id:
            continue
        operon_to_indexes.setdefault(operon_id, []).append(idx)

    for indexes in operon_to_indexes.values():
        if len(indexes) < 2:
            continue

        scores = [concordance_scores[i] for i in indexes]
        avg = sum(scores) / len(scores)
        high_conf = len([s for s in scores if s >= 70])

        if avg >= 80:
            avg_component = 50
        elif avg >= 60:
            avg_component = 40
        elif avg >= 40:
            avg_component = 30
        elif avg >= 20:
            avg_component = 20
        else:
            avg_component = 10

        high_pct = (high_conf / len(scores)) * 100.0
        if high_pct >= 75:
            high_component = 30
        elif high_pct >= 50:
            high_component = 25
        elif high_pct >= 25:
            high_component = 15
        else:
            high_component = 5

        variance = sum((x - avg) ** 2 for x in scores) / len(scores)
        stddev = variance ** 0.5
        if stddev <= 10:
            consistency_component = 20
        elif stddev <= 20:
            consistency_component = 15
        elif stddev <= 30:
            consistency_component = 10
        else:
            consistency_component = 5

        total = min(100, avg_component + high_component + consistency_component)
        for i in indexes:
            coreg_scores[i] = total

    return coreg_scores


def apply_scoring(consolidated_file: Path, model_file: Path, output_file: Path) -> None:
    with open(model_file, "r", encoding="utf-8") as handle:
        model = json.load(handle)

    patterns = model.get("patterns", {})
    total_genomes = int(model.get("metadata", {}).get("total_genomes", 1))

    df = pd.read_csv(consolidated_file, sep="\t", low_memory=False)
    columns = list(df.columns)

    concordance_scores = []
    hit_counts = []
    signatures = []
    pattern_frequency = []
    genomes_with_pattern = []
    operon_consistency = []

    for _, row in df.iterrows():
        gene_set = extract_gene_set(row, TOOL_COLUMN_MAP, columns)
        operon_info = get_operon_info(row, columns)

        hit_count = len(gene_set)
        signature = build_signature(gene_set)
        score = calculate_concordance_score(signature, patterns, gene_set, operon_info, total_genomes)

        hit_counts.append(hit_count)
        signatures.append(signature if signature else "NO_ANNOTATIONS")
        concordance_scores.append(score)

        if signature and signature in patterns:
            pdata = patterns[signature]
            pattern_frequency.append(int(pdata.get("total_count", 0)))
            genomes_with_pattern.append(int(pdata.get("num_genomes", 0)))

            operon_freq = float(pdata.get("operon_frequency", 0.0))
            if operon_info["has_operon"] and operon_freq >= 0.7:
                operon_consistency.append("CONSISTENT_OPERON")
            elif (not operon_info["has_operon"]) and operon_freq <= 0.3:
                operon_consistency.append("CONSISTENT_NON_OPERON")
            elif 0.3 < operon_freq < 0.7:
                operon_consistency.append("MIXED_PATTERN")
            else:
                operon_consistency.append("INCONSISTENT")
        else:
            pattern_frequency.append(0)
            genomes_with_pattern.append(0)
            operon_consistency.append("NOVEL_PATTERN" if signature else "NO_ANNOTATIONS")

    coreg_scores = build_operon_coregulation_scores(df, concordance_scores)

    df["concordance_hit_count"] = hit_counts
    df["concordance_signature"] = signatures
    df["concordance_pattern_frequency"] = pattern_frequency
    df["concordance_genomes_with_pattern"] = genomes_with_pattern
    df["concordance_operon_consistency"] = operon_consistency
    df["concordance_score"] = concordance_scores
    df["coregulation_score"] = coreg_scores

    output_file.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_file, sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Apply concordance scoring")
    parser.add_argument("--consolidated", required=True, help="Input consolidated_annotations.tsv")
    parser.add_argument("--model", required=True, help="Concordance model JSON")
    parser.add_argument("--output", required=True, help="Output scored TSV path")
    args = parser.parse_args()

    apply_scoring(Path(args.consolidated), Path(args.model), Path(args.output))
    print(f"Scored concordance output written to {args.output}")


if __name__ == "__main__":
    main()
