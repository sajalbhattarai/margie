#!/usr/bin/env bash
# Run concordance model build (if needed) and score one consolidated genome table.

set -euo pipefail

WORKSPACE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$WORKSPACE_ROOT"

GENOME=""
CONSOLIDATED_FILE=""
CONSOLIDATED_DIR="output/consolidated"
OUTPUT_ROOT="output/scoring"
REBUILD_MODEL=0

usage() {
  cat <<EOF
Usage: $0 [options]

Required:
  --genome NAME                Genome name
  --consolidated FILE          Path to consolidated_annotations.tsv for genome

Optional:
  --consolidated-dir DIR       Root directory of consolidated outputs (default: output/consolidated)
  --output-root DIR            Root output directory for scoring (default: output/scoring)
  --rebuild-model              Force model rebuild before scoring
  -h, --help                   Show help

Outputs:
  <output-root>/model/concordance_reference.json
  <output-root>/<genome>/scoring.tsv
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --genome) GENOME="$2"; shift 2 ;;
    --consolidated) CONSOLIDATED_FILE="$2"; shift 2 ;;
    --consolidated-dir) CONSOLIDATED_DIR="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    --rebuild-model) REBUILD_MODEL=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "$GENOME" || -z "$CONSOLIDATED_FILE" ]]; then
  echo "ERROR: --genome and --consolidated are required"
  usage
  exit 1
fi

if [[ ! -f "$CONSOLIDATED_FILE" ]]; then
  echo "ERROR: Consolidated input not found: $CONSOLIDATED_FILE"
  exit 1
fi

MODEL_FILE="$OUTPUT_ROOT/model/concordance_reference.json"
SCORED_DIR="$OUTPUT_ROOT/$GENOME"
SCORED_FILE="$SCORED_DIR/scoring.tsv"
LOG_DIR="$SCORED_DIR/logs"
LOG_FILE="$LOG_DIR/concordance.log"

mkdir -p "$SCORED_DIR" "$LOG_DIR" "$OUTPUT_ROOT/model"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

log "Concordance scoring started for $GENOME"
log "Consolidated input: $CONSOLIDATED_FILE"
log "Model file: $MODEL_FILE"

if [[ "$REBUILD_MODEL" -eq 1 || ! -f "$MODEL_FILE" ]]; then
  log "Building concordance model from $CONSOLIDATED_DIR"
  python3 processing/scripts/scoring/build_concordance_model.py \
    --consolidated-dir "$CONSOLIDATED_DIR" \
    --output "$MODEL_FILE" >> "$LOG_FILE" 2>&1
  log "Model build complete"
else
  log "Reusing existing concordance model"
fi

log "Applying concordance scoring"
python3 processing/scripts/scoring/apply_concordance_scoring.py \
  --consolidated "$CONSOLIDATED_FILE" \
  --model "$MODEL_FILE" \
  --output "$SCORED_FILE" >> "$LOG_FILE" 2>&1

if [[ ! -f "$SCORED_FILE" ]]; then
  log "ERROR: Scoring output was not created"
  exit 1
fi

log "Concordance scoring complete: $SCORED_FILE"
