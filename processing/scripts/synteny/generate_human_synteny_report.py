#!/usr/bin/env python3
"""Generate human-readable synteny reports with BLAST and operon perspectives."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


BLAST_COLS = [
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


def pair_key(g1: str, g2: str) -> str:
    a, b = sorted((str(g1), str(g2)))
    return f"{a}||{b}"


def get_genome_labels(df: pd.DataFrame) -> tuple[str, str]:
    if df.empty:
        return "Genome_A", "Genome_B"
    a_vals = sorted(df["genome_a"].dropna().astype(str).unique().tolist())
    b_vals = sorted(df["genome_b"].dropna().astype(str).unique().tolist())
    a_label = a_vals[0] if len(a_vals) == 1 else "Genome_A"
    b_label = b_vals[0] if len(b_vals) == 1 else "Genome_B"
    return a_label, b_label


def pick_existing(*paths: Path) -> Path:
    for path in paths:
        if path.exists():
            return path
    return paths[0]


def infer_synteny_root(input_path: Path) -> Path:
    if input_path.parent.name == "processed":
        return input_path.parent.parent
    return input_path.parent


def clear_header_tables(inter: pd.DataFrame, blocks: pd.DataFrame, operon_pairs: pd.DataFrame, operon_ready: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pair_rename = {
        "block_id": "Synteny_Block_ID",
        "gene_a": "Genome_A_Feature_ID",
        "genome_a": "Genome_A_Name",
        "contig_a": "Genome_A_Contig_ID",
        "start_a": "Genome_A_Start_bp",
        "end_a": "Genome_A_End_bp",
        "strand_a": "Genome_A_Strand",
        "midpoint_a": "Genome_A_Midpoint_bp",
        "gene_b": "Genome_B_Feature_ID",
        "genome_b": "Genome_B_Name",
        "contig_b": "Genome_B_Contig_ID",
        "start_b": "Genome_B_Start_bp",
        "end_b": "Genome_B_End_bp",
        "strand_b": "Genome_B_Strand",
        "midpoint_b": "Genome_B_Midpoint_bp",
        "anchor_evalue": "MCScanX_Anchor_Evalue",
        "blast_evalue": "BLASTP_Evalue",
        "blast_identity_pct": "BLASTP_Percent_Identity",
        "blast_bitscore": "BLASTP_BitScore",
        "blast_aln_len": "BLASTP_Alignment_Length_aa",
        "gene_a_in_operon": "Genome_A_In_Operon",
        "gene_b_in_operon": "Genome_B_In_Operon",
        "both_operonic": "Both_Genes_In_Operons",
        "operon_pair_id": "Operon_Pair_ID",
        "op_operon_id_a": "Genome_A_Operon_ID",
        "op_operon_position_a": "Genome_A_Operon_Position",
        "op_operon_size_a": "Genome_A_Operon_Size",
        "op_operon_id_b": "Genome_B_Operon_ID",
        "op_operon_position_b": "Genome_B_Operon_Position",
        "op_operon_size_b": "Genome_B_Operon_Size",
        "operon_conservation_degree": "Operon_Conservation_Degree_pct",
    }

    block_rename = {
        "block_id": "Synteny_Block_ID",
        "genome_a": "Genome_A_Name",
        "genome_b": "Genome_B_Name",
        "pairs": "Matched_Gene_Pairs",
        "a_start_min": "Genome_A_Block_Start_bp",
        "a_end_max": "Genome_A_Block_End_bp",
        "b_start_min": "Genome_B_Block_Start_bp",
        "b_end_max": "Genome_B_Block_End_bp",
        "mean_identity_pct": "Mean_BLASTP_Identity_pct",
    }

    op_pair_rename = {
        "operon_pair_id": "Operon_Pair_ID",
        "genome_a": "Genome_A_Name",
        "operon_id_a": "Genome_A_Operon_ID",
        "operon_size_a": "Genome_A_Operon_Size",
        "mapped_genes_a": "Genome_A_Mapped_Genes",
        "coverage_a": "Genome_A_Operon_Coverage",
        "genome_b": "Genome_B_Name",
        "operon_id_b": "Genome_B_Operon_ID",
        "operon_size_b": "Genome_B_Operon_Size",
        "mapped_genes_b": "Genome_B_Mapped_Genes",
        "coverage_b": "Genome_B_Operon_Coverage",
        "mapped_pairs": "Mapped_Gene_Pairs",
        "mean_identity_pct": "Mean_BLASTP_Identity_pct",
        "operon_conservation_degree": "Operon_Conservation_Degree_pct",
    }

    inter_clear = inter.rename(columns={k: v for k, v in pair_rename.items() if k in inter.columns})
    blocks_clear = blocks.rename(columns={k: v for k, v in block_rename.items() if k in blocks.columns})
    operon_pairs_clear = operon_pairs.rename(columns={k: v for k, v in op_pair_rename.items() if k in operon_pairs.columns})
    operon_ready_clear = operon_ready.rename(columns={k: v for k, v in pair_rename.items() if k in operon_ready.columns})
    return inter_clear, blocks_clear, operon_pairs_clear, operon_ready_clear


def load_blast_agg(blast_raw_path: Path) -> pd.DataFrame:
    if not blast_raw_path.exists():
        return pd.DataFrame(columns=["pair_key", "blast_evalue", "blast_identity_pct", "blast_bitscore", "blast_aln_len"])

    blast = pd.read_csv(blast_raw_path, sep="\t", header=None, names=BLAST_COLS)
    if blast.empty:
        return pd.DataFrame(columns=["pair_key", "blast_evalue", "blast_identity_pct", "blast_bitscore", "blast_aln_len"])

    blast["pair_key"] = [pair_key(a, b) for a, b in zip(blast["qseqid"], blast["sseqid"])]
    agg = (
        blast.groupby("pair_key", as_index=False)
        .agg(
            blast_evalue=("evalue", "min"),
            blast_identity_pct=("pident", "max"),
            blast_bitscore=("bitscore", "max"),
            blast_aln_len=("length", "max"),
        )
    )
    return agg


def load_operons(operon_root: Path, genomes: list[str]) -> pd.DataFrame:
    frames = []
    for g in genomes:
        p = operon_root / g / "operons.tsv"
        if not p.exists():
            continue
        try:
            df = pd.read_csv(p, sep="\t", dtype=str)
        except Exception:
            continue

        required = {"feature_id", "operon_id", "operon_size", "operon_position"}
        if not required.issubset(set(df.columns)):
            continue

        if "member_genes" not in df.columns:
            df["member_genes"] = ""
        if "contig" not in df.columns:
            df["contig"] = ""
        if "start" not in df.columns:
            df["start"] = ""
        if "end" not in df.columns:
            df["end"] = ""
        if "strand" not in df.columns:
            df["strand"] = ""

        df["genome"] = g
        frames.append(df)

    if not frames:
        return pd.DataFrame(columns=["feature_id", "operon_id", "operon_size", "operon_position", "member_genes", "contig", "start", "end", "strand", "genome"])

    all_ops = pd.concat(frames, ignore_index=True)
    all_ops = all_ops.drop_duplicates(subset=["feature_id", "genome"], keep="first")
    return all_ops


def add_blast_and_operon_context(inter: pd.DataFrame, blast_agg: pd.DataFrame, operons: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    out = inter.copy()

    # BLAST enrichment
    out["pair_key"] = [pair_key(a, b) for a, b in zip(out["gene_a"], out["gene_b"])]
    out = out.merge(blast_agg, on="pair_key", how="left")

    # Operon enrichment
    op_a = operons.rename(
        columns={
            "feature_id": "gene_a",
            "genome": "genome_a",
            "operon_id": "op_operon_id_a",
            "operon_size": "op_operon_size_a",
            "operon_position": "op_operon_position_a",
            "member_genes": "op_member_genes_a",
            "contig": "op_contig_a",
            "start": "op_start_a",
            "end": "op_end_a",
            "strand": "op_strand_a",
        }
    )
    op_b = operons.rename(
        columns={
            "feature_id": "gene_b",
            "genome": "genome_b",
            "operon_id": "op_operon_id_b",
            "operon_size": "op_operon_size_b",
            "operon_position": "op_operon_position_b",
            "member_genes": "op_member_genes_b",
            "contig": "op_contig_b",
            "start": "op_start_b",
            "end": "op_end_b",
            "strand": "op_strand_b",
        }
    )
    out = out.merge(op_a, on=["gene_a", "genome_a"], how="left")
    out = out.merge(op_b, on=["gene_b", "genome_b"], how="left")

    out["gene_a_in_operon"] = out["op_operon_id_a"].fillna("NA") != "NA"
    out["gene_b_in_operon"] = out["op_operon_id_b"].fillna("NA") != "NA"
    out["both_operonic"] = out["gene_a_in_operon"] & out["gene_b_in_operon"]

    out["operon_pair_id"] = ""
    mask = out["both_operonic"]
    out.loc[mask, "operon_pair_id"] = (
        out.loc[mask, "genome_a"]
        + "::"
            + out.loc[mask, "op_operon_id_a"].astype(str)
        + "<->"
        + out.loc[mask, "genome_b"]
        + "::"
            + out.loc[mask, "op_operon_id_b"].astype(str)
    )

    # Operon pair conservation summary
    both = out[out["both_operonic"]].copy()
    if both.empty:
        operon_pairs = pd.DataFrame(
            columns=[
                "operon_pair_id",
                "genome_a",
                "operon_id_a",
                "operon_size_a",
                "mapped_genes_a",
                "coverage_a",
                "genome_b",
                "operon_id_b",
                "operon_size_b",
                "mapped_genes_b",
                "coverage_b",
                "mapped_pairs",
                "mean_identity_pct",
                "operon_conservation_degree",
            ]
        )
    else:
        both["operon_size_a_num"] = pd.to_numeric(both["op_operon_size_a"], errors="coerce")
        both["operon_size_b_num"] = pd.to_numeric(both["op_operon_size_b"], errors="coerce")

        operon_pairs = (
            both.groupby(["operon_pair_id", "genome_a", "op_operon_id_a", "genome_b", "op_operon_id_b"], as_index=False)
            .agg(
                operon_size_a=("operon_size_a_num", "max"),
                operon_size_b=("operon_size_b_num", "max"),
                mapped_genes_a=("gene_a", "nunique"),
                mapped_genes_b=("gene_b", "nunique"),
                mapped_pairs=("gene_a", "size"),
                mean_identity_pct=("blast_identity_pct", "mean"),
            )
        )
        operon_pairs = operon_pairs.rename(columns={"op_operon_id_a": "operon_id_a", "op_operon_id_b": "operon_id_b"})
        operon_pairs["coverage_a"] = operon_pairs["mapped_genes_a"] / operon_pairs["operon_size_a"].replace(0, pd.NA)
        operon_pairs["coverage_b"] = operon_pairs["mapped_genes_b"] / operon_pairs["operon_size_b"].replace(0, pd.NA)
        operon_pairs["operon_conservation_degree"] = 100.0 * operon_pairs[["coverage_a", "coverage_b"]].min(axis=1)
        operon_pairs = operon_pairs.sort_values(["operon_conservation_degree", "mapped_pairs"], ascending=[False, False])

        out = out.merge(
            operon_pairs[["operon_pair_id", "operon_conservation_degree", "mapped_pairs"]].rename(
                columns={"mapped_pairs": "operon_pair_mapped_pairs"}
            ),
            on="operon_pair_id",
            how="left",
        )

    # Operon stats
    operon_stats = []
    operon_stats.append(("pairs_with_gene_a_in_operon", int(out["gene_a_in_operon"].sum())))
    operon_stats.append(("pairs_with_gene_b_in_operon", int(out["gene_b_in_operon"].sum())))
    operon_stats.append(("pairs_with_both_genes_operonic", int(out["both_operonic"].sum())))
    operon_stats.append(("unique_operon_pairs_with_synteny", int(operon_pairs["operon_pair_id"].nunique()) if not operon_pairs.empty else 0))
    if not operon_pairs.empty:
        operon_stats.append(("mean_operon_conservation_degree", float(operon_pairs["operon_conservation_degree"].mean())))
        operon_stats.append(("max_operon_conservation_degree", float(operon_pairs["operon_conservation_degree"].max())))
    else:
        operon_stats.append(("mean_operon_conservation_degree", 0.0))
        operon_stats.append(("max_operon_conservation_degree", 0.0))

    operon_stats_df = pd.DataFrame(operon_stats, columns=["metric", "value"])
    return out, operon_pairs, operon_stats_df


def summarize_pairs(df: pd.DataFrame, operon_stats_df: pd.DataFrame) -> pd.DataFrame:
    stats = []
    stats.append(("total_pairs", len(df)))
    stats.append(("total_blocks", df["block_id"].nunique()))
    stats.append(("genome_a_count", df["genome_a"].nunique()))
    stats.append(("genome_b_count", df["genome_b"].nunique()))
    stats.append(("unique_genes_a", df["gene_a"].nunique()))
    stats.append(("unique_genes_b", df["gene_b"].nunique()))

    if len(df) > 0:
        stats.append(("median_pairs_per_block", int(df.groupby("block_id").size().median())))
        stats.append(("max_pairs_in_block", int(df.groupby("block_id").size().max())))
        if "blast_identity_pct" in df.columns:
            stats.append(("mean_identity_pct", float(df["blast_identity_pct"].mean())))
            stats.append(("median_identity_pct", float(df["blast_identity_pct"].median())))
    else:
        stats.append(("median_pairs_per_block", 0))
        stats.append(("max_pairs_in_block", 0))
        stats.append(("mean_identity_pct", 0.0))
        stats.append(("median_identity_pct", 0.0))

    stats_df = pd.DataFrame(stats, columns=["metric", "value"])
    return pd.concat([stats_df, operon_stats_df], ignore_index=True)


def block_summary(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=["block_id", "genome_a", "genome_b", "pairs", "a_start_min", "a_end_max", "b_start_min", "b_end_max", "mean_identity_pct"])

    grouped = (
        df.groupby(["block_id", "genome_a", "genome_b"], as_index=False)
        .agg(
            pairs=("gene_a", "size"),
            a_start_min=("start_a", "min"),
            a_end_max=("end_a", "max"),
            b_start_min=("start_b", "min"),
            b_end_max=("end_b", "max"),
            mean_identity_pct=("blast_identity_pct", "mean"),
        )
        .sort_values(["pairs", "block_id"], ascending=[False, True])
    )
    return grouped


def write_text_report(out_path: Path, stats: pd.DataFrame, blocks: pd.DataFrame, operon_pairs: pd.DataFrame, df: pd.DataFrame) -> None:
    top_blocks = blocks.head(10)
    top_operon_pairs = operon_pairs.head(10)
    with out_path.open("w", encoding="utf-8") as f:
        f.write("Synteny + Operon Conservation Report\n")
        f.write("===================================\n\n")
        f.write("Summary statistics\n")
        for _, row in stats.iterrows():
            f.write(f"- {row['metric']}: {row['value']}\n")

        f.write("\nGenome comparisons present\n")
        for _, row in df[["genome_a", "genome_b"]].drop_duplicates().iterrows():
            f.write(f"- {row['genome_a']}  <->  {row['genome_b']}\n")

        f.write("\nTop 10 synteny blocks by pair count\n")
        if top_blocks.empty:
            f.write("- none\n")
        else:
            for _, row in top_blocks.iterrows():
                f.write(
                    f"- block {row['block_id']}: pairs={row['pairs']}, "
                    f"mean_identity={row.get('mean_identity_pct', float('nan')):.2f}%, "
                    f"A[{row['a_start_min']}-{row['a_end_max']}], "
                    f"B[{row['b_start_min']}-{row['b_end_max']}]\n"
                )

        f.write("\nTop 10 operon-pair conservation groups\n")
        if top_operon_pairs.empty:
            f.write("- none (operon outputs unavailable or no operonic synteny)\n")
        else:
            for _, row in top_operon_pairs.iterrows():
                f.write(
                    f"- {row['operon_pair_id']}: degree={row['operon_conservation_degree']:.2f}%, "
                    f"pairs={int(row['mapped_pairs'])}, covA={row['coverage_a']:.2f}, covB={row['coverage_b']:.2f}\n"
                )


def draw_dotplot(df: pd.DataFrame, out_png: Path, genome_a_label: str, genome_b_label: str) -> None:
    fig, ax = plt.subplots(figsize=(9, 7))
    if df.empty:
        ax.text(0.5, 0.5, "No inter-genome synteny pairs", ha="center", va="center")
        ax.axis("off")
    else:
        ax.scatter(df["midpoint_a"], df["midpoint_b"], s=6, alpha=0.35)
        ax.set_xlabel(f"{genome_a_label} gene midpoint (bp)")
        ax.set_ylabel(f"{genome_b_label} gene midpoint (bp)")
        ax.set_title(f"Inter-genome synteny dotplot: {genome_a_label} vs {genome_b_label}")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def draw_linear_map(df: pd.DataFrame, blocks: pd.DataFrame, out_png: Path, genome_a_label: str, genome_b_label: str) -> None:
    fig, ax = plt.subplots(figsize=(12, 6))
    if df.empty or blocks.empty:
        ax.text(0.5, 0.5, "No inter-genome synteny pairs", ha="center", va="center")
        ax.axis("off")
        fig.tight_layout()
        fig.savefig(out_png, dpi=180)
        plt.close(fig)
        return

    top_block_ids = blocks.head(10)["block_id"].tolist()
    sub = df[df["block_id"].isin(top_block_ids)].copy()

    a_min, a_max = sub["midpoint_a"].min(), sub["midpoint_a"].max()
    b_min, b_max = sub["midpoint_b"].min(), sub["midpoint_b"].max()

    ax.hlines(1.0, a_min, a_max, linewidth=4)
    ax.hlines(0.0, b_min, b_max, linewidth=4)

    color_map = {bid: i for i, bid in enumerate(top_block_ids)}
    for _, r in sub.iterrows():
        c = f"C{color_map[r['block_id']] % 10}"
        ax.plot([r["midpoint_a"], r["midpoint_b"]], [1.0, 0.0], color=c, alpha=0.15, linewidth=0.7)

    ax.text(a_min, 1.08, genome_a_label, fontsize=10)
    ax.text(b_min, -0.16, genome_b_label, fontsize=10)
    ax.set_title(f"Synteny genome map (top 10 blocks): {genome_a_label} vs {genome_b_label}")
    ax.set_yticks([])
    ax.set_xlabel("Genomic position (gene midpoint)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate human-readable synteny+operon report package")
    parser.add_argument("--input", required=True, help="synteny_gene_pairs_detailed.tsv")
    parser.add_argument("--outdir", required=True, help="output report directory")
    parser.add_argument("--blast-raw", default="", help="Optional mcscanx_all.blast.raw.tsv for identity enrichment")
    parser.add_argument("--operon-root", default="", help="Optional output/operon root path")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    input_path = Path(args.input)
    df = pd.read_csv(input_path, sep="\t")
    synteny_root = infer_synteny_root(input_path)

    inter = df[df["genome_a"] != df["genome_b"]].copy()
    inter = inter.sort_values(["block_id", "genome_a", "contig_a", "start_a", "start_b"])

    blast_raw = Path(args.blast_raw) if args.blast_raw else pick_existing(
        synteny_root / "native" / "mcscanx_all.blast.raw.tsv",
        synteny_root / "mcscanx_all.blast.raw.tsv",
    )
    blast_agg = load_blast_agg(blast_raw)

    genomes = sorted(set(inter["genome_a"].dropna().tolist()) | set(inter["genome_b"].dropna().tolist()))
    operon_root = Path(args.operon_root) if args.operon_root else synteny_root.parent / "operon"
    operons = load_operons(operon_root, genomes)

    inter, operon_pairs, operon_stats_df = add_blast_and_operon_context(inter, blast_agg, operons)

    inter.to_csv(outdir / "intergenome_synteny_pairs.tsv", sep="\t", index=False)

    stats = summarize_pairs(inter, operon_stats_df)
    stats.to_csv(outdir / "summary_stats.tsv", sep="\t", index=False)

    blocks = block_summary(inter)
    blocks.to_csv(outdir / "block_summary.tsv", sep="\t", index=False)

    operon_ready_cols = [
        "block_id",
        "gene_a",
        "genome_a",
        "contig_a",
        "start_a",
        "end_a",
        "strand_a",
        "gene_b",
        "genome_b",
        "contig_b",
        "start_b",
        "end_b",
        "strand_b",
        "blast_identity_pct",
        "blast_bitscore",
        "op_operon_id_a",
        "op_operon_position_a",
        "op_operon_size_a",
        "op_operon_id_b",
        "op_operon_position_b",
        "op_operon_size_b",
        "operon_conservation_degree",
    ]
    present_cols = [c for c in operon_ready_cols if c in inter.columns]
    operon_ready = inter[present_cols].copy() if not inter.empty else pd.DataFrame(columns=present_cols)
    operon_ready.to_csv(outdir / "operon_join_ready.tsv", sep="\t", index=False)

    operon_pairs.to_csv(outdir / "operon_pair_conservation.tsv", sep="\t", index=False)
    operon_stats_df.to_csv(outdir / "operon_conservation_stats.tsv", sep="\t", index=False)

    inter_clear, blocks_clear, operon_pairs_clear, operon_ready_clear = clear_header_tables(inter, blocks, operon_pairs, operon_ready)
    inter_clear.to_csv(outdir / "intergenome_synteny_pairs.clear_headers.tsv", sep="\t", index=False)
    blocks_clear.to_csv(outdir / "block_summary.clear_headers.tsv", sep="\t", index=False)
    operon_ready_clear.to_csv(outdir / "operon_join_ready.clear_headers.tsv", sep="\t", index=False)
    operon_pairs_clear.to_csv(outdir / "operon_pair_conservation.clear_headers.tsv", sep="\t", index=False)

    summary_clear = stats.rename(columns={"metric": "Metric", "value": "Value"})
    summary_clear.to_csv(outdir / "summary_stats.clear_headers.tsv", sep="\t", index=False)

    genome_a_label, genome_b_label = get_genome_labels(inter)

    write_text_report(outdir / "report.txt", stats, blocks, operon_pairs, inter)
    draw_dotplot(inter, outdir / "synteny_dotplot_intergenome.png", genome_a_label, genome_b_label)
    draw_linear_map(inter, blocks, outdir / "synteny_genome_map_top_blocks.png", genome_a_label, genome_b_label)


if __name__ == "__main__":
    main()
