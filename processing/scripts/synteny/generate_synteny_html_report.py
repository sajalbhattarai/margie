#!/usr/bin/env python3
"""Build a polished standalone HTML report for synteny results with clear headers."""

from __future__ import annotations

import argparse
import base64
from pathlib import Path
from typing import Optional

import pandas as pd


def image_to_data_uri(path: Path) -> str:
    if not path.exists():
        return ""
    mime = "image/png"
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{mime};base64,{data}"


def fmt_int(v: object) -> str:
    try:
        return f"{int(float(v)):,}"
    except Exception:
        return str(v)


def pick_existing(*paths: Path) -> Path:
    for p in paths:
        if p.exists():
            return p
    return paths[0]


def build_summary_cards(summary_df: pd.DataFrame) -> str:
    metric_col = "Metric" if "Metric" in summary_df.columns else "metric"
    value_col = "Value" if "Value" in summary_df.columns else "value"
    values = {str(getattr(r, metric_col)): getattr(r, value_col) for r in summary_df.itertuples(index=False)}
    cards = [
        ("Inter-genome Pairs", values.get("total_pairs", 0)),
        ("Synteny Blocks", values.get("total_blocks", 0)),
        ("Unique Genes (Genome A)", values.get("unique_genes_a", 0)),
        ("Unique Genes (Genome B)", values.get("unique_genes_b", 0)),
        ("Mean % Identity", values.get("mean_identity_pct", 0)),
        ("Median Pairs / Block", values.get("median_pairs_per_block", 0)),
        ("Max Pairs in Block", values.get("max_pairs_in_block", 0)),
        ("Both Genes Operonic", values.get("pairs_with_both_genes_operonic", 0)),
    ]
    out = []
    for title, val in cards:
        if "Identity" in title:
            try:
                txt = f"{float(val):.2f}%"
            except Exception:
                txt = str(val)
        else:
            txt = fmt_int(val)
        out.append(f"<div class='card'><div class='card-title'>{title}</div><div class='card-value'>{txt}</div></div>")
    return "\n".join(out)


def to_html_table(df: pd.DataFrame, table_id: str, max_rows: Optional[int]) -> str:
    if df.empty:
        return "<p class='muted'>No records available.</p>"
    shown = df if not max_rows else df.head(max_rows)
    html = shown.to_html(index=False, classes=f"data-table {table_id}", border=0)
    if max_rows and len(df) > max_rows:
        html += f"<p class='muted'>Showing first {max_rows:,} of {len(df):,} rows.</p>"
    return html


def keep(df: pd.DataFrame, preferred: list[str]) -> pd.DataFrame:
    cols = [c for c in preferred if c in df.columns]
    return df[cols] if cols else df


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate HTML report for synteny outputs")
    parser.add_argument("--report-dir", required=True, help="Directory containing report TSV files")
    parser.add_argument("--output", required=True, help="Path to output HTML file")
    parser.add_argument("--max-rows", type=int, default=0, help="Maximum rows per table; 0 means all rows")
    args = parser.parse_args()

    report_dir = Path(args.report_dir)
    out_path = Path(args.output)
    max_rows: Optional[int] = args.max_rows if args.max_rows > 0 else None

    summary_tsv = pick_existing(report_dir / "summary_stats.clear_headers.tsv", report_dir / "summary_stats.tsv")
    block_tsv = pick_existing(report_dir / "block_summary.clear_headers.tsv", report_dir / "block_summary.tsv")
    inter_tsv = pick_existing(report_dir / "intergenome_synteny_pairs.clear_headers.tsv", report_dir / "intergenome_synteny_pairs.tsv")
    operon_tsv = pick_existing(report_dir / "operon_join_ready.clear_headers.tsv", report_dir / "operon_join_ready.tsv")
    operon_pair_tsv = pick_existing(report_dir / "operon_pair_conservation.clear_headers.tsv", report_dir / "operon_pair_conservation.tsv")
    operon_stats_tsv = report_dir / "operon_conservation_stats.tsv"
    dotplot_png = report_dir / "synteny_dotplot_intergenome.png"
    map_png = report_dir / "synteny_genome_map_top_blocks.png"

    for p in [summary_tsv, block_tsv, inter_tsv]:
        if not p.exists():
            raise FileNotFoundError(f"Missing required input: {p}")

    summary = pd.read_csv(summary_tsv, sep="\t")
    blocks = pd.read_csv(block_tsv, sep="\t")
    pairs = pd.read_csv(inter_tsv, sep="\t")
    operon = pd.read_csv(operon_tsv, sep="\t") if operon_tsv.exists() else pd.DataFrame()
    operon_pairs = pd.read_csv(operon_pair_tsv, sep="\t") if operon_pair_tsv.exists() else pd.DataFrame()
    operon_stats = pd.read_csv(operon_stats_tsv, sep="\t") if operon_stats_tsv.exists() else pd.DataFrame()

    pair_cmp = "n/a"
    if {"Genome_A_Name", "Genome_B_Name"}.issubset(set(pairs.columns)) and not pairs.empty:
        uniq = pairs[["Genome_A_Name", "Genome_B_Name"]].drop_duplicates()
        pair_cmp = f"{uniq.iloc[0]['Genome_A_Name']} <-> {uniq.iloc[0]['Genome_B_Name']}"
    elif {"genome_a", "genome_b"}.issubset(set(pairs.columns)) and not pairs.empty:
        uniq = pairs[["genome_a", "genome_b"]].drop_duplicates()
        pair_cmp = f"{uniq.iloc[0]['genome_a']} <-> {uniq.iloc[0]['genome_b']}"

    dotplot_data = image_to_data_uri(dotplot_png)
    map_data = image_to_data_uri(map_png)

    col_help = "Genome_A_* columns always refer to the left/reference genome; Genome_B_* columns refer to the matched genome."

    blocks_html = to_html_table(
        keep(
            blocks,
            [
                "Synteny_Block_ID",
                "Genome_A_Name",
                "Genome_B_Name",
                "Matched_Gene_Pairs",
                "Genome_A_Block_Start_bp",
                "Genome_A_Block_End_bp",
                "Genome_B_Block_Start_bp",
                "Genome_B_Block_End_bp",
                "Mean_BLASTP_Identity_pct",
            ],
        ),
        "blocks",
        max_rows,
    )

    pairs_html = to_html_table(
        keep(
            pairs,
            [
                "Synteny_Block_ID",
                "Genome_A_Name",
                "Genome_A_Feature_ID",
                "Genome_A_Contig_ID",
                "Genome_A_Start_bp",
                "Genome_A_End_bp",
                "Genome_A_Strand",
                "Genome_B_Name",
                "Genome_B_Feature_ID",
                "Genome_B_Contig_ID",
                "Genome_B_Start_bp",
                "Genome_B_End_bp",
                "Genome_B_Strand",
                "MCScanX_Anchor_Evalue",
                "BLASTP_Percent_Identity",
                "BLASTP_BitScore",
                "BLASTP_Evalue",
                "Genome_A_In_Operon",
                "Genome_B_In_Operon",
                "Both_Genes_In_Operons",
                "Genome_A_Operon_ID",
                "Genome_B_Operon_ID",
                "Operon_Conservation_Degree_pct",
            ],
        ),
        "pairs",
        max_rows,
    )

    operon_html = to_html_table(operon, "operon", max_rows)
    operon_pairs_html = to_html_table(operon_pairs, "operon-pairs", max_rows)
    operon_stats_html = to_html_table(operon_stats, "operon-stats", max_rows)

    html = f"""
<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\" />
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />
  <title>Synteny Report</title>
  <style>
    :root {{
      --bg: #f7f8f3;
      --ink: #12201f;
      --muted: #586665;
      --card: #ffffff;
      --accent: #0f766e;
      --border: #d8dfda;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      color: var(--ink);
      background:
        radial-gradient(1200px 500px at -15% -10%, #dbeee9 0%, transparent 60%),
        radial-gradient(900px 420px at 105% -15%, #f7e7cc 0%, transparent 55%),
        var(--bg);
      font-family: ui-serif, Georgia, Cambria, 'Times New Roman', serif;
      line-height: 1.45;
    }}
    .wrap {{ max-width: 1280px; margin: 0 auto; padding: 28px 22px 48px; }}
    .hero {{
      background: linear-gradient(135deg, #083b37 0%, #145a55 60%, #176f67 100%);
      color: #f7fffd;
      border-radius: 18px;
      padding: 24px 24px;
      box-shadow: 0 16px 34px rgba(8,59,55,0.25);
      margin-bottom: 18px;
    }}
    .hero h1 {{ margin: 0 0 6px; font-size: 32px; letter-spacing: 0.2px; }}
    .hero p {{ margin: 0; color: #d8eeea; }}
    .badge {{
      display: inline-block;
      margin-top: 10px;
      padding: 6px 10px;
      border-radius: 999px;
      background: rgba(255,255,255,0.14);
      color: #fff;
      font-size: 12px;
      font-weight: 600;
    }}
    .cards {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 12px;
      margin: 18px 0 20px;
    }}
    .card {{
      background: var(--card);
      border: 1px solid var(--border);
      border-radius: 14px;
      padding: 12px 14px;
      box-shadow: 0 4px 12px rgba(9, 35, 34, 0.06);
    }}
    .card-title {{ font-size: 12px; color: var(--muted); text-transform: uppercase; letter-spacing: 0.5px; }}
    .card-value {{ font-size: 28px; font-weight: 700; margin-top: 4px; color: var(--accent); }}
    .section {{
      background: #fff;
      border: 1px solid var(--border);
      border-radius: 14px;
      padding: 14px 14px;
      margin-bottom: 14px;
      box-shadow: 0 3px 10px rgba(9,35,34,0.05);
    }}
    h2 {{ margin: 0 0 8px; font-size: 22px; }}
    .muted {{ color: var(--muted); font-size: 13px; }}
    .plot-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(360px, 1fr));
      gap: 12px;
      align-items: start;
    }}
    .plot img {{ width: 100%; border: 1px solid var(--border); border-radius: 10px; background: #fff; }}
    .plot h3 {{ margin: 0 0 8px; font-size: 16px; color: #204341; }}
    .data-table {{
      width: 100%; border-collapse: collapse; font-size: 12px;
      border: 1px solid var(--border); border-radius: 10px; overflow: hidden;
      display: block; overflow-x: auto; max-height: 680px;
    }}
    .data-table thead th {{
      position: sticky; top: 0;
      background: #eff7f5;
      border-bottom: 1px solid var(--border);
      text-align: left;
      padding: 7px 8px;
      white-space: nowrap;
    }}
    .data-table tbody td {{ border-bottom: 1px solid #edf1ee; padding: 6px 8px; white-space: nowrap; }}
    .downloads a {{
      display: inline-block; margin: 6px 10px 0 0; padding: 8px 10px;
      border-radius: 9px; border: 1px solid var(--border);
      text-decoration: none; color: #103432; background: #f4faf8;
      font-size: 13px;
    }}
    .footer-note {{ font-size: 12px; color: var(--muted); margin-top: 10px; }}
  </style>
</head>
<body>
  <div class=\"wrap\">
    <section class=\"hero\">
      <h1>Genome Synteny Report</h1>
      <p>Human-readable inter-genome synteny summary with explicit genome-aware headers and coordinate-level mappings.</p>
      <span class=\"badge\">Comparison: {pair_cmp}</span>
    </section>

    <section class=\"cards\">{build_summary_cards(summary)}</section>

    <section class=\"section\">
      <h2>Genome Map Visualizations</h2>
      <p class=\"muted\">These maps highlight genome-to-genome syntenic structure and block continuity.</p>
      <div class=\"plot-grid\">
        <div class=\"plot\">
          <h3>Inter-genome Dotplot</h3>
          {('<img src="' + dotplot_data + '" alt="dotplot" />') if dotplot_data else '<p class="muted">Dotplot image missing.</p>'}
        </div>
        <div class=\"plot\">
          <h3>Linear Genome Map (Top Blocks)</h3>
          {('<img src="' + map_data + '" alt="linear synteny map" />') if map_data else '<p class="muted">Linear map image missing.</p>'}
        </div>
      </div>
    </section>

    <section class=\"section\">
      <h2>Top Synteny Blocks</h2>
      <p class=\"muted\">{col_help}</p>
      {blocks_html}
    </section>

    <section class=\"section\">
      <h2>Detailed Gene-to-Gene Coordinate Matches</h2>
      <p class=\"muted\">Feature IDs can look similar, so each row includes explicit genome name and A/B-prefixed headers.</p>
      {pairs_html}
    </section>

    <section class=\"section\">
      <h2>Operon Join-Ready Table</h2>
      <p class=\"muted\">Explicit A/B genome headers with operon membership and conservation fields.</p>
      {operon_html}
    </section>

    <section class=\"section\">
      <h2>Operon Pair Conservation</h2>
      {operon_pairs_html}
    </section>

    <section class=\"section\">
      <h2>Operon Conservation Metrics</h2>
      {operon_stats_html}
    </section>

    <section class=\"section\">
      <h2>Downloads</h2>
      <div class=\"downloads\">
        <a href=\"summary_stats.tsv\">summary_stats.tsv</a>
        <a href=\"summary_stats.clear_headers.tsv\">summary_stats.clear_headers.tsv</a>
        <a href=\"block_summary.tsv\">block_summary.tsv</a>
        <a href=\"block_summary.clear_headers.tsv\">block_summary.clear_headers.tsv</a>
        <a href=\"intergenome_synteny_pairs.tsv\">intergenome_synteny_pairs.tsv</a>
        <a href=\"intergenome_synteny_pairs.clear_headers.tsv\">intergenome_synteny_pairs.clear_headers.tsv</a>
        <a href=\"operon_join_ready.tsv\">operon_join_ready.tsv</a>
        <a href=\"operon_join_ready.clear_headers.tsv\">operon_join_ready.clear_headers.tsv</a>
        <a href=\"operon_pair_conservation.tsv\">operon_pair_conservation.tsv</a>
        <a href=\"operon_pair_conservation.clear_headers.tsv\">operon_pair_conservation.clear_headers.tsv</a>
        <a href=\"operon_conservation_stats.tsv\">operon_conservation_stats.tsv</a>
        <a href=\"synteny_dotplot_intergenome.png\">synteny_dotplot_intergenome.png</a>
        <a href=\"synteny_genome_map_top_blocks.png\">synteny_genome_map_top_blocks.png</a>
      </div>
      <p class=\"footer-note\">Use --max-rows N if browser rendering gets slow with very large tables.</p>
    </section>
  </div>
</body>
</html>
"""

    out_path.write_text(html, encoding="utf-8")
    print(f"Wrote HTML report: {out_path}")


if __name__ == "__main__":
    main()
