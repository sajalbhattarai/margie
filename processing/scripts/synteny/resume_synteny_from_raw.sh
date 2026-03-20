#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${1:-/container/output/synteny}"
PREFIX="${2:-mcscanx_all}"
MIN_IDENTITY="${3:-25}"
MIN_ALN_LEN="${4:-40}"
MAX_HITS_PER_QUERY="${5:-30}"

pick_existing() {
  local path
  for path in "$@"; do
    if [[ -e "$path" ]]; then
      printf '%s\n' "$path"
      return 0
    fi
  done
  printf '%s\n' "$1"
  return 1
}

NATIVE_DIR="$OUTDIR/native"
PROCESSED_DIR="$OUTDIR/processed"
REPORT_DIR="$PROCESSED_DIR/human_report"
CONSOLIDATED_DIR="/container/output/consolidation"
ANNOTATION_JOIN="$PROCESSED_DIR/synteny_annotation_comparison.tsv"
mkdir -p "$NATIVE_DIR" "$PROCESSED_DIR" "$REPORT_DIR"

MC_BASE_DIR="$(dirname "$(pick_existing "$NATIVE_DIR/${PREFIX}.gff" "$OUTDIR/${PREFIX}.gff")")"
RAW="$(pick_existing "$MC_BASE_DIR/${PREFIX}.blast.raw.tsv" "$OUTDIR/${PREFIX}.blast.raw.tsv")"
LOCI="$(pick_existing "$MC_BASE_DIR/${PREFIX}.gene_loci.tsv" "$OUTDIR/${PREFIX}.gene_loci.tsv")"
MCBLAST="$MC_BASE_DIR/${PREFIX}.blast"
LOG="$MC_BASE_DIR/${PREFIX}.log"
PAIRS="$OUTDIR/synteny_gene_pairs.tsv"
PAIRSD="$PROCESSED_DIR/synteny_gene_pairs_detailed.tsv"
COL="$MC_BASE_DIR/${PREFIX}.collinearity"

if [[ ! -s "$RAW" ]]; then
  echo "Missing raw blast: $RAW"
  exit 1
fi
if [[ ! -s "$LOCI" ]]; then
  echo "Missing loci: $LOCI"
  exit 1
fi

{
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Resume: filtering existing DIAMOND output"
} >> "$LOG"

TMP1="$MC_BASE_DIR/${PREFIX}.blast.prefilter.tsv"
TMP2="$MC_BASE_DIR/${PREFIX}.blast.sorted.tsv"

awk -F'\t' -v OFS='\t' -v minid="$MIN_IDENTITY" -v minlen="$MIN_ALN_LEN" '
  NR==FNR {
    if (FNR == 1) next
    g[$1] = $2
    next
  }
  {
    if ($1 == $2) next
    if (($3 + 0) < (minid + 0)) next
    if (($4 + 0) < (minlen + 0)) next
    if (g[$1] == "" || g[$2] == "") next
    if (g[$1] == g[$2]) next
    print $0
  }
' "$LOCI" "$RAW" > "$TMP1"

sort -T "$MC_BASE_DIR" -S 1G -t $'\t' -k1,1 -k12,12gr -k11,11g "$TMP1" > "$TMP2"
awk -F'\t' -v OFS='\t' -v topn="$MAX_HITS_PER_QUERY" '
  {
    c[$1]++
    if (c[$1] <= topn) print $0
  }
' "$TMP2" > "$MCBLAST"

rm -f "$TMP1" "$TMP2"

if [[ ! -s "$MCBLAST" ]]; then
  echo "Filtered BLAST file is empty: $MCBLAST"
  exit 1
fi

{
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Resume: running MCScanX"
} >> "$LOG"
MCScanX -s 2 -m 50 -w 10 "$MC_BASE_DIR/$PREFIX" >> "$LOG" 2>&1

python3 /container/scripts/synteny_pairs_from_collinearity.py --collinearity "$COL" --output "$PAIRS" >> "$LOG" 2>&1
python3 /container/scripts/annotate_synteny_pairs_with_loci.py --pairs "$PAIRS" --loci "$LOCI" --output "$PAIRSD" >> "$LOG" 2>&1
if [[ -d "$CONSOLIDATED_DIR" ]]; then
  python3 /container/work/processing/containers/synteny/scripts/join_synteny_with_annotations.py \
    --pairs "$PAIRS" \
    --consolidated-dir "$CONSOLIDATED_DIR" \
    --output "$ANNOTATION_JOIN" >> "$LOG" 2>&1 || true
fi
python3 /container/scripts/render_synteny_plots.py --pairs "$PAIRS" --output-dir "$PROCESSED_DIR" >> "$LOG" 2>&1 || true
python3 /container/workspace/processing/scripts/synteny/generate_human_synteny_report.py --input "$PAIRSD" --outdir "$REPORT_DIR" --blast-raw "$RAW" --operon-root /container/output/operon >> "$LOG" 2>&1 || true
python3 /container/workspace/processing/scripts/synteny/generate_synteny_html_report.py --report-dir "$REPORT_DIR" --output "$REPORT_DIR/report.html" >> "$LOG" 2>&1 || true
python3 /container/workspace/processing/scripts/synteny/generate_synteny_gene_report.py --pairs-detailed "$PAIRSD" --rasttk-root /container/output/rasttk --operon-root /container/output/operon --output-json "$REPORT_DIR/gene_synteny_operon_report.json" --output-txt "$REPORT_DIR/gene_synteny_operon_report.txt" --max-links-per-gene 100 --max-genes-in-txt 300 >> "$LOG" 2>&1 || true

{
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Resume pipeline completed"
} >> "$LOG"
