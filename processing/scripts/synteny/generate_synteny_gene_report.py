#!/usr/bin/env python3
"""Generate scalable gene-centric synteny reports as JSON and TXT.

This reporter is designed for large runs (hundreds/thousands of genomes) and
focuses on per-gene records rather than giant flat tables.

Core synteny in this report is based only on positional collinearity evidence
from MCScanX anchors/blocks and sequence homology evidence from DIAMOND/BLAST.
Operon information is optional and used only as an additional bonus context.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd


def load_rast_functions(rasttk_root: Path, genomes: set[str]) -> dict[tuple[str, str], str]:
    fn_map: dict[tuple[str, str], str] = {}
    for genome in sorted(genomes):
        p = rasttk_root / genome / "rast.tsv"
        if not p.exists():
            continue
        try:
            df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
        except Exception:
            continue
        if "feature_id" not in df.columns:
            continue
        desc_col = "RAST_description" if "RAST_description" in df.columns else None
        if desc_col is None:
            continue
        for row in df[["feature_id", desc_col]].itertuples(index=False):
            fid = str(row[0])
            desc = "" if pd.isna(row[1]) else str(row[1])
            fn_map[(genome, fid)] = desc
    return fn_map


def load_operon_tables(operon_root: Path, genomes: set[str]) -> tuple[dict[tuple[str, str], dict[str, Any]], dict[tuple[str, str], set[str]]]:
    operon_by_gene: dict[tuple[str, str], dict[str, Any]] = {}
    operon_members: dict[tuple[str, str], set[str]] = {}

    for genome in sorted(genomes):
        p = operon_root / genome / "operons.tsv"
        if not p.exists():
            continue
        try:
            df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
        except Exception:
            continue
        if "feature_id" not in df.columns or "operon_id" not in df.columns:
            continue

        for row in df.itertuples(index=False):
            r = row._asdict()
            fid = str(r.get("feature_id", ""))
            op_id = str(r.get("operon_id", ""))
            members_raw = str(r.get("member_genes", ""))
            members = {m.strip() for m in members_raw.split(",") if m.strip()}
            if not members and fid:
                members = {fid}

            if op_id:
                operon_members[(genome, op_id)] = members

            if fid:
                operon_by_gene[(genome, fid)] = {
                    "operon_id": op_id,
                    "operon_size": r.get("operon_size", ""),
                    "operon_position": r.get("operon_position", ""),
                    "member_genes": sorted(members),
                    "contig": r.get("contig", ""),
                    "start": r.get("start", ""),
                    "end": r.get("end", ""),
                    "strand": r.get("strand", ""),
                }

    return operon_by_gene, operon_members


def safe_float(v: Any) -> float | None:
    try:
        if pd.isna(v):
            return None
        return float(v)
    except Exception:
        return None


def directed_edges(df: pd.DataFrame) -> list[dict[str, Any]]:
    edges: list[dict[str, Any]] = []
    for r in df.itertuples(index=False):
        d = r._asdict()

        edge_ab = {
            "query_gene": d.get("gene_a", ""),
            "query_genome": d.get("genome_a", ""),
            "query_contig": d.get("contig_a", ""),
            "query_start": d.get("start_a", ""),
            "query_end": d.get("end_a", ""),
            "query_strand": d.get("strand_a", ""),
            "target_gene": d.get("gene_b", ""),
            "target_genome": d.get("genome_b", ""),
            "target_contig": d.get("contig_b", ""),
            "target_start": d.get("start_b", ""),
            "target_end": d.get("end_b", ""),
            "target_strand": d.get("strand_b", ""),
            "block_id": d.get("block_id", ""),
            "anchor_evalue": d.get("anchor_evalue", ""),
            "blast_identity_pct": safe_float(d.get("blast_identity_pct")),
            "blast_bitscore": safe_float(d.get("blast_bitscore")),
            "blast_evalue": safe_float(d.get("blast_evalue")),
        }
        edge_ba = {
            "query_gene": d.get("gene_b", ""),
            "query_genome": d.get("genome_b", ""),
            "query_contig": d.get("contig_b", ""),
            "query_start": d.get("start_b", ""),
            "query_end": d.get("end_b", ""),
            "query_strand": d.get("strand_b", ""),
            "target_gene": d.get("gene_a", ""),
            "target_genome": d.get("genome_a", ""),
            "target_contig": d.get("contig_a", ""),
            "target_start": d.get("start_a", ""),
            "target_end": d.get("end_a", ""),
            "target_strand": d.get("strand_a", ""),
            "block_id": d.get("block_id", ""),
            "anchor_evalue": d.get("anchor_evalue", ""),
            "blast_identity_pct": safe_float(d.get("blast_identity_pct")),
            "blast_bitscore": safe_float(d.get("blast_bitscore")),
            "blast_evalue": safe_float(d.get("blast_evalue")),
        }
        if edge_ab["query_gene"] and edge_ab["target_gene"]:
            edges.append(edge_ab)
        if edge_ba["query_gene"] and edge_ba["target_gene"]:
            edges.append(edge_ba)

    return edges


def build_gene_report(
    edges: list[dict[str, Any]],
    fn_map: dict[tuple[str, str], str],
    operon_by_gene: dict[tuple[str, str], dict[str, Any]],
    operon_members: dict[tuple[str, str], set[str]],
    max_links_per_gene: int,
) -> dict[str, Any]:
    gene_edges: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    # Fast lookup for "where does a given gene map in a given target genome"
    map_lookup: dict[tuple[str, str, str], set[str]] = defaultdict(set)

    for e in edges:
        gk = (str(e["query_genome"]), str(e["query_gene"]))
        gene_edges[gk].append(e)
        map_lookup[(str(e["query_genome"]), str(e["query_gene"]), str(e["target_genome"]))].add(str(e["target_gene"]))

    genes_payload: list[dict[str, Any]] = []

    for (genome, gene), glinks in sorted(gene_edges.items(), key=lambda x: (x[0][0], x[0][1])):
        first = glinks[0]
        q_operon = operon_by_gene.get((genome, gene), {})
        q_op_id = str(q_operon.get("operon_id", ""))
        q_members = operon_members.get((genome, q_op_id), set()) if q_op_id else set()
        q_mates = sorted([m for m in q_members if m != gene])

        cotranscribed_with = []
        for mate in q_mates:
            mate_op = operon_by_gene.get((genome, mate), {})
            cotranscribed_with.append(
                {
                    "gene": mate,
                    "function": fn_map.get((genome, mate), ""),
                    "contig": mate_op.get("contig", ""),
                    "start": mate_op.get("start", ""),
                    "end": mate_op.get("end", ""),
                    "strand": mate_op.get("strand", ""),
                }
            )

        syntenic_with = []
        for e in glinks[:max_links_per_gene]:
            tgt_genome = str(e["target_genome"])
            tgt_gene = str(e["target_gene"])
            tgt_op = operon_by_gene.get((tgt_genome, tgt_gene), {})
            tgt_op_id = str(tgt_op.get("operon_id", ""))
            tgt_members = operon_members.get((tgt_genome, tgt_op_id), set()) if tgt_op_id else set()

            conserved_mates = []
            clustered_elsewhere = []
            for mate in q_mates:
                mapped_targets = map_lookup.get((genome, mate, tgt_genome), set())
                in_same_target_operon = sorted([x for x in mapped_targets if x in tgt_members])
                if in_same_target_operon:
                    conserved_mates.append(
                        {
                            "query_operon_mate": mate,
                            "target_operon_mate_candidates": in_same_target_operon,
                        }
                    )
                elif mapped_targets:
                    clustered_elsewhere.append(
                        {
                            "query_operon_mate": mate,
                            "mapped_elsewhere_in_target_genome": sorted(mapped_targets),
                        }
                    )

            denom = max(1, len(q_mates))
            operon_conservation = {
                "query_operon_id": q_op_id,
                "target_operon_id": tgt_op_id,
                "query_operon_size": q_operon.get("operon_size", ""),
                "target_operon_size": tgt_op.get("operon_size", ""),
                "query_operon_mates_total": len(q_mates),
                "query_operon_mates_conserved_with_target_operon": len(conserved_mates),
                "query_operon_mate_conservation_pct": round(100.0 * len(conserved_mates) / denom, 2),
                "conserved_mates": conserved_mates,
                "mates_clustered_elsewhere": clustered_elsewhere,
                "both_genes_part_of_operons": bool(q_op_id and tgt_op_id),
            }

            syntenic_with.append(
                {
                    "target_gene": tgt_gene,
                    "target_genome": tgt_genome,
                    "target_function": fn_map.get((tgt_genome, tgt_gene), ""),
                    "target_position": {
                        "contig": e.get("target_contig", ""),
                        "start": e.get("target_start", ""),
                        "end": e.get("target_end", ""),
                        "strand": e.get("target_strand", ""),
                    },
                    # Core synteny evidence (independent of operon membership/strand constraints)
                    "core_synteny_evidence": {
                        "block_id": e.get("block_id", ""),
                        "anchor_evalue": e.get("anchor_evalue", ""),
                        "blast_identity_pct": e.get("blast_identity_pct"),
                        "blast_bitscore": e.get("blast_bitscore"),
                        "blast_evalue": e.get("blast_evalue"),
                    },
                    "target_cotranscribed_with": sorted([m for m in tgt_members if m != tgt_gene]),
                    # Optional bonus layer; does not affect block detection/order consistency.
                    "optional_operon_bonus": operon_conservation,
                }
            )

        avg_conservation = None
        if syntenic_with:
            vals = [
                s["optional_operon_bonus"]["query_operon_mate_conservation_pct"]
                for s in syntenic_with
                if s["optional_operon_bonus"]["query_operon_mates_total"] > 0
            ]
            if vals:
                avg_conservation = round(sum(vals) / len(vals), 2)

        genes_payload.append(
            {
                "gene": gene,
                "housing_genome": genome,
                "location_in_genome": {
                    "contig": first.get("query_contig", ""),
                    "start": first.get("query_start", ""),
                    "end": first.get("query_end", ""),
                    "strand": first.get("query_strand", ""),
                },
                "function": fn_map.get((genome, gene), ""),
                "operon_context": {
                    "operon_id": q_op_id,
                    "operon_size": q_operon.get("operon_size", ""),
                    "operon_position": q_operon.get("operon_position", ""),
                },
                "cotranscribed_with": cotranscribed_with,
                "syntenic_with": syntenic_with,
                "operon_conservation_summary": {
                    "syntenic_links_reported": len(syntenic_with),
                    "query_operon_mates_total": len(q_mates),
                    "average_query_operon_mate_conservation_pct": avg_conservation,
                },
            }
        )

    return {
        "schema_version": "1.0",
        "description": "Gene-centric synteny and operon-conservation report with genomic coordinates and function annotations.",
        "core_synteny_principles": {
            "genomic_position": True,
            "relative_order": True,
            "neighborhood_conservation": True,
            "requires_same_strand": False,
            "requires_operon_membership": False,
            "assumes_functional_linkage": False,
        },
        "optional_bonus_layers": {
            "operon_aware_scoring": True,
            "operon_layer_is_non_blocking": True,
        },
        "genes": genes_payload,
    }


def write_txt_summary(payload: dict[str, Any], out_txt: Path, max_genes: int) -> None:
    genes = payload.get("genes", [])
    principles = payload.get("core_synteny_principles", {})
    bonus = payload.get("optional_bonus_layers", {})
    with out_txt.open("w", encoding="utf-8") as f:
        f.write("Gene-Centric Synteny and Operon Conservation Report\n")
        f.write("==================================================\n\n")
        f.write("Method Summary\n")
        f.write("--------------\n")
        f.write("Core synteny is inferred from positional collinearity and neighborhood order conservation.\n")
        f.write(f"- Genomic position used: {principles.get('genomic_position', True)}\n")
        f.write(f"- Relative order used: {principles.get('relative_order', True)}\n")
        f.write(f"- Neighborhood conservation used: {principles.get('neighborhood_conservation', True)}\n")
        f.write(f"- Requires same strand: {principles.get('requires_same_strand', False)}\n")
        f.write(f"- Requires operon membership: {principles.get('requires_operon_membership', False)}\n")
        f.write(f"- Assumes functional linkage: {principles.get('assumes_functional_linkage', False)}\n")
        f.write(f"- Operon-aware scoring is optional bonus only: {bonus.get('operon_aware_scoring', True)}\n")
        f.write(f"- Operon layer is non-blocking for synteny calls: {bonus.get('operon_layer_is_non_blocking', True)}\n\n")
        f.write(f"Total genes in JSON payload: {len(genes)}\n")
        f.write("This report describes per-gene syntenic matches, positions, and operon conservation behavior.\n\n")

        for i, g in enumerate(genes):
            if i >= max_genes:
                f.write(f"\n... truncated in TXT after {max_genes} genes. Full details are in JSON.\n")
                break

            loc = g.get("location_in_genome", {})
            f.write(f"GENE: {g.get('gene', '')}\n")
            f.write(f"  Housing_genome: {g.get('housing_genome', '')}\n")
            f.write(f"  Location_in_genome: {loc.get('contig', '')}:{loc.get('start', '')}-{loc.get('end', '')}\n")
            f.write(f"  Strand: {loc.get('strand', '')}\n")
            f.write(f"  Function: {g.get('function', '')}\n")
            op = g.get("operon_context", {})
            f.write(f"  Operon: {op.get('operon_id', '')} (size={op.get('operon_size', '')}, position={op.get('operon_position', '')})\n")

            cotx = g.get("cotranscribed_with", [])
            f.write("  cotranscribed_with:\n")
            if not cotx:
                f.write("    - none\n")
            else:
                for c in cotx:
                    f.write(
                        f"    - {c.get('gene', '')} | {c.get('contig', '')}:{c.get('start', '')}-{c.get('end', '')} ({c.get('strand', '')}) | {c.get('function', '')}\n"
                    )

            f.write("  syntenic_with:\n")
            syn = g.get("syntenic_with", [])
            if not syn:
                f.write("    - none\n")
            else:
                for s in syn:
                    p = s.get("target_position", {})
                    oc = s.get("optional_operon_bonus", {})
                    f.write(
                        f"    - {s.get('target_gene', '')} in {s.get('target_genome', '')} at "
                        f"{p.get('contig', '')}:{p.get('start', '')}-{p.get('end', '')} ({p.get('strand', '')}); "
                        f"identity={s.get('core_synteny_evidence', {}).get('blast_identity_pct', '')}; "
                        f"query_operon_mate_conservation={oc.get('query_operon_mate_conservation_pct', '')}%\n"
                    )
                    conserved = oc.get("conserved_mates", [])
                    moved = oc.get("mates_clustered_elsewhere", [])
                    if conserved:
                        f.write("      conserved_operon_mates:\n")
                        for c in conserved:
                            f.write(
                                f"        - {c.get('query_operon_mate', '')} -> {', '.join(c.get('target_operon_mate_candidates', []))}\n"
                            )
                    if moved:
                        f.write("      operon_mates_clustered_elsewhere:\n")
                        for m in moved:
                            f.write(
                                f"        - {m.get('query_operon_mate', '')} -> {', '.join(m.get('mapped_elsewhere_in_target_genome', []))}\n"
                            )
            f.write("\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate gene-centric synteny JSON/TXT reports")
    parser.add_argument("--pairs-detailed", required=True, help="Path to synteny_gene_pairs_detailed.tsv")
    parser.add_argument("--rasttk-root", required=True, help="Path to output/rasttk")
    parser.add_argument("--operon-root", required=True, help="Path to output/operon")
    parser.add_argument("--output-json", required=True, help="Output JSON path")
    parser.add_argument("--output-txt", required=True, help="Output TXT path")
    parser.add_argument("--max-links-per-gene", type=int, default=100, help="Cap syntenic links retained per gene in JSON")
    parser.add_argument("--max-genes-in-txt", type=int, default=200, help="Cap number of genes fully rendered in TXT")
    args = parser.parse_args()

    pairs_path = Path(args.pairs_detailed)
    if not pairs_path.exists():
        raise FileNotFoundError(f"Missing pairs file: {pairs_path}")

    df = pd.read_csv(pairs_path, sep="\t", low_memory=False)
    if "genome_a" not in df.columns or "genome_b" not in df.columns:
        raise ValueError("pairs-detailed TSV missing expected genome columns")

    inter = df[df["genome_a"] != df["genome_b"]].copy()
    genomes = set(inter["genome_a"].dropna().astype(str).tolist()) | set(inter["genome_b"].dropna().astype(str).tolist())

    fn_map = load_rast_functions(Path(args.rasttk_root), genomes)
    operon_by_gene, operon_members = load_operon_tables(Path(args.operon_root), genomes)

    edges = directed_edges(inter)
    payload = build_gene_report(
        edges=edges,
        fn_map=fn_map,
        operon_by_gene=operon_by_gene,
        operon_members=operon_members,
        max_links_per_gene=max(1, args.max_links_per_gene),
    )

    out_json = Path(args.output_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    out_txt = Path(args.output_txt)
    out_txt.parent.mkdir(parents=True, exist_ok=True)
    write_txt_summary(payload, out_txt, max_genes=max(1, args.max_genes_in_txt))

    print(f"Wrote JSON report: {out_json}")
    print(f"Wrote TXT report: {out_txt}")
    print(f"Genes in JSON: {len(payload.get('genes', []))}")


if __name__ == "__main__":
    main()
