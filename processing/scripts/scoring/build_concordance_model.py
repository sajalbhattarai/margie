#!/usr/bin/env python3
"""Build concordance reference model from consolidated annotation files."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import pandas as pd

from concordance_common import (
    build_signature,
    extract_gene_set,
    get_operon_info,
    now_utc_iso,
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


def discover_inputs(consolidated_dir: Path):
    return sorted(consolidated_dir.glob("*/consolidated_annotations.tsv"))


def build_model(consolidated_dir: Path, output_file: Path) -> dict:
    files = discover_inputs(consolidated_dir)
    if not files:
        raise FileNotFoundError(f"No consolidated annotations found under {consolidated_dir}")

    pattern_data = {}

    total_genes_processed = 0
    genome_list = []

    for path in files:
        genome = path.parent.name
        genome_list.append(genome)

        df = pd.read_csv(path, sep="\t", low_memory=False)
        columns = list(df.columns)

        for _, row in df.iterrows():
            gene_set = extract_gene_set(row, TOOL_COLUMN_MAP, columns)
            if not gene_set:
                continue

            signature = build_signature(gene_set)
            operon_info = get_operon_info(row, columns)
            hits = len(gene_set)

            if signature not in pattern_data:
                pattern_data[signature] = {
                    "total_count": 0,
                    "genome_counts": {},
                    "operon_occurrences": 0,
                    "non_operon_occurrences": 0,
                    "operon_sizes": [],
                    "hit_counts": [],
                    "example_genes": [],
                }

            pdata = pattern_data[signature]
            pdata["total_count"] += 1
            pdata["genome_counts"][genome] = pdata["genome_counts"].get(genome, 0) + 1
            pdata["hit_counts"].append(hits)

            if operon_info["has_operon"]:
                pdata["operon_occurrences"] += 1
                try:
                    operon_size = int(float(str(operon_info.get("operon_size", 0))))
                except ValueError:
                    operon_size = 0
                if operon_size > 0:
                    pdata["operon_sizes"].append(operon_size)
            else:
                pdata["non_operon_occurrences"] += 1

            if len(pdata["example_genes"]) < 5:
                pdata["example_genes"].append(
                    {
                        "genome": genome,
                        "gene_id": str(row.get("feature_id", "")),
                        "operon": operon_info,
                    }
                )

            total_genes_processed += 1

    patterns = {}
    for signature, pdata in pattern_data.items():
        total_count = pdata["total_count"]
        operon_occurrences = pdata["operon_occurrences"]
        hit_counts = pdata["hit_counts"]
        operon_sizes = pdata["operon_sizes"]

        patterns[signature] = {
            "total_count": total_count,
            "genome_counts": dict(pdata["genome_counts"]),
            "num_genomes": len(pdata["genome_counts"]),
            "operon_occurrences": operon_occurrences,
            "non_operon_occurrences": pdata["non_operon_occurrences"],
            "operon_frequency": (operon_occurrences / total_count) if total_count else 0.0,
            "avg_operon_size": (sum(operon_sizes) / len(operon_sizes)) if operon_sizes else 0.0,
            "avg_hit_count": (sum(hit_counts) / len(hit_counts)) if hit_counts else 0.0,
            "max_hit_count": max(hit_counts) if hit_counts else 0,
            "min_hit_count": min(hit_counts) if hit_counts else 0,
            "example_genes": pdata["example_genes"],
        }

    model = {
        "metadata": {
            "created_at_utc": now_utc_iso(),
            "consolidated_dir": str(consolidated_dir),
            "total_genes_processed": total_genes_processed,
            "total_genomes": len(files),
            "total_patterns": len(patterns),
            "genome_list": genome_list,
        },
        "patterns": patterns,
    }

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w", encoding="utf-8") as handle:
        json.dump(model, handle, indent=2)

    return model


def main() -> None:
    parser = argparse.ArgumentParser(description="Build concordance reference model")
    parser.add_argument("--consolidated-dir", required=True, help="Path to output/consolidated")
    parser.add_argument("--output", required=True, help="Output model JSON path")
    args = parser.parse_args()

    model = build_model(Path(args.consolidated_dir), Path(args.output))
    print(
        f"Built concordance model with {model['metadata']['total_patterns']} patterns "
        f"from {model['metadata']['total_genomes']} genomes"
    )


if __name__ == "__main__":
    main()
