#!/usr/bin/env bash
# Standalone multi-genome synteny workflow using MCScanX + DIAMOND.

set -euo pipefail

WORKSPACE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$WORKSPACE_ROOT"

RASTTK_DIR="output/rasttk"
CONSOLIDATED_DIR="output/consolidation"
OUTDIR="output/synteny"
THREADS="8"
GENOME_LIST=""
PREFIX="mcscanx_all"

usage() {
  cat <<EOF
Usage: $0 [options]

Standalone synteny workflow for existing RASTtk outputs.

Options:
  --rasttk-dir DIR         RASTtk genome output root (default: output/rasttk)
  --consolidated-dir DIR   Consolidated annotation root (default: output/consolidation)
  --outdir DIR             Output root (default: output/synteny)
  --genomes CSV            Comma-separated genome names to include (default: auto-discover)
  --threads N              Threads for DIAMOND (default: 8)
  --prefix NAME            MCScanX prefix basename (default: mcscanx_all)
  -h, --help               Show this help

Expected per-genome files:
  <rasttk-dir>/<genome>/gene_calls/<genome>.faa
  <rasttk-dir>/<genome>/gene_calls/<genome>.gff

Requirements:
  - diamond
  - MCScanX (binary on PATH)
  - python3

Outputs:
  <outdir>/<prefix>.gff
  <outdir>/<prefix>.blast
  <outdir>/<prefix>.collinearity      (from MCScanX)
  <outdir>/synteny_gene_pairs.tsv
  <outdir>/synteny_annotation_comparison.tsv  (if consolidated files are available)
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --rasttk-dir) RASTTK_DIR="$2"; shift 2 ;;
    --consolidated-dir) CONSOLIDATED_DIR="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --genomes) GENOME_LIST="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

command -v diamond >/dev/null 2>&1 || { echo "ERROR: diamond not found in PATH"; exit 1; }
command -v MCScanX >/dev/null 2>&1 || { echo "ERROR: MCScanX not found in PATH"; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found in PATH"; exit 1; }

mkdir -p "$OUTDIR"

TMP_DIR="$OUTDIR/tmp"
mkdir -p "$TMP_DIR"

COMBINED_FAA="$OUTDIR/${PREFIX}.faa"
MC_GFF="$OUTDIR/${PREFIX}.gff"
RAW_BLAST="$OUTDIR/${PREFIX}.blast.raw.tsv"
MC_BLAST="$OUTDIR/${PREFIX}.blast"
PAIRS_TSV="$OUTDIR/synteny_gene_pairs.tsv"
ANNOTATION_JOIN_TSV="$OUTDIR/synteny_annotation_comparison.tsv"
LOG_FILE="$OUTDIR/${PREFIX}.log"

: > "$LOG_FILE"
: > "$COMBINED_FAA"
: > "$MC_GFF"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

log "Starting MCScanX synteny workflow"
log "Workspace: $WORKSPACE_ROOT"
log "RASTtk dir: $RASTTK_DIR"
log "Consolidated dir: $CONSOLIDATED_DIR"
log "Output dir: $OUTDIR"

if [[ ! -d "$RASTTK_DIR" ]]; then
  log "ERROR: RASTtk directory not found: $RASTTK_DIR"
  exit 1
fi

# Resolve genome list.
declare -a GENOMES
if [[ -n "$GENOME_LIST" ]]; then
  IFS=',' read -r -a GENOMES <<< "$GENOME_LIST"
else
  while IFS= read -r d; do
    GENOMES+=("$(basename "$d")")
  done < <(find "$RASTTK_DIR" -mindepth 1 -maxdepth 1 -type d | sort)
fi

if [[ ${#GENOMES[@]} -lt 2 ]]; then
  log "ERROR: Need at least 2 genomes for synteny analysis"
  exit 1
fi

log "Genomes selected (${#GENOMES[@]}): ${GENOMES[*]}"

VALID_GENOMES=0
for genome in "${GENOMES[@]}"; do
  faa="$RASTTK_DIR/$genome/gene_calls/$genome.faa"
  gff="$RASTTK_DIR/$genome/gene_calls/$genome.gff"

  if [[ ! -f "$faa" ]]; then
    log "WARNING: Missing FAA, skipping genome: $genome"
    continue
  fi
  if [[ ! -f "$gff" ]]; then
    log "WARNING: Missing GFF, skipping genome: $genome"
    continue
  fi

  VALID_GENOMES=$((VALID_GENOMES + 1))

  # Concatenate proteins for all-vs-all homology search.
  cat "$faa" >> "$COMBINED_FAA"

  # Convert GFF3 CDS rows to MCScanX GFF format:
  # gene_id, chromosome, start, end
  awk -F'\t' -v OFS='\t' -v genome="$genome" '
    BEGIN { }
    /^#/ { next }
    NF < 9 { next }
    $3 != "CDS" { next }
    {
      id = ""
      n = split($9, attrs, ";")
      for (i = 1; i <= n; i++) {
        if (attrs[i] ~ /^ID=/) {
          sub(/^ID=/, "", attrs[i])
          id = attrs[i]
          break
        }
      }
      if (id != "") {
        chr = genome "|" $1
        print id, chr, $4, $5
      }
    }
  ' "$gff" >> "$MC_GFF"
done

if [[ $VALID_GENOMES -lt 2 ]]; then
  log "ERROR: Fewer than 2 genomes had both FAA and GFF"
  exit 1
fi

if [[ ! -s "$COMBINED_FAA" || ! -s "$MC_GFF" ]]; then
  log "ERROR: Required MCScanX inputs were not generated"
  exit 1
fi

log "Building DIAMOND database"
diamond makedb --in "$COMBINED_FAA" -d "$OUTDIR/${PREFIX}" >> "$LOG_FILE" 2>&1

log "Running DIAMOND all-vs-all BLASTP"
diamond blastp \
  --query "$COMBINED_FAA" \
  --db "$OUTDIR/${PREFIX}.dmnd" \
  --out "$RAW_BLAST" \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
  --threads "$THREADS" \
  --max-target-seqs 25 \
  --evalue 1e-5 >> "$LOG_FILE" 2>&1

if [[ ! -s "$RAW_BLAST" ]]; then
  log "ERROR: DIAMOND output is empty"
  exit 1
fi

# MCScanX consumes BLAST tabular-like input; keep top 12 columns.
cp "$RAW_BLAST" "$MC_BLAST"

log "Running MCScanX"
MCScanX "$OUTDIR/$PREFIX" >> "$LOG_FILE" 2>&1 || {
  log "ERROR: MCScanX failed"
  exit 1
}

COLLINEARITY="$OUTDIR/${PREFIX}.collinearity"
if [[ ! -f "$COLLINEARITY" ]]; then
  log "ERROR: Missing MCScanX collinearity output: $COLLINEARITY"
  exit 1
fi

log "Extracting gene pairs from collinearity output"
python3 "processing/scripts/synteny/synteny_pairs_from_collinearity.py" \
  --collinearity "$COLLINEARITY" \
  --output "$PAIRS_TSV" >> "$LOG_FILE" 2>&1

if [[ ! -s "$PAIRS_TSV" ]]; then
  log "WARNING: No syntenic pairs extracted"
fi

if [[ -d "$CONSOLIDATED_DIR" ]]; then
  log "Joining synteny pairs with consolidated annotations"
  python3 "processing/scripts/synteny/join_synteny_with_annotations.py" \
    --pairs "$PAIRS_TSV" \
    --consolidated-dir "$CONSOLIDATED_DIR" \
    --output "$ANNOTATION_JOIN_TSV" >> "$LOG_FILE" 2>&1 || {
      log "WARNING: Annotation join step failed; synteny pairs still available"
    }
else
  log "WARNING: Consolidated directory missing, skipping annotation join"
fi

log "Workflow complete"
log "Main outputs:"
log "  $MC_GFF"
log "  $MC_BLAST"
log "  $COLLINEARITY"
log "  $PAIRS_TSV"
[[ -f "$ANNOTATION_JOIN_TSV" ]] && log "  $ANNOTATION_JOIN_TSV"
