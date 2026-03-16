#!/bin/bash
# UniProt annotation using UPIMAPI against Swiss-Prot resources

set -euo pipefail

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_faa> <output_tsv> <swissprot_resource> <threads>"
    exit 1
fi

INPUT_FAA="$1"
OUTPUT_TSV="$2"
SWISSPROT_RESOURCE="$3"
THREADS="$4"

if [ ! -f "$INPUT_FAA" ]; then
    echo "Error: Input FASTA not found: $INPUT_FAA"
    exit 1
fi

if [ ! -f "$SWISSPROT_RESOURCE" ]; then
    echo "Error: Swiss-Prot resource not found: $SWISSPROT_RESOURCE"
    exit 1
fi

OUTPUT_BASE_DIR="$(dirname "$(dirname "$OUTPUT_TSV")")"
LOG_FILE="${OUTPUT_BASE_DIR}/uniprot_upimapi.log"
OUTPUT_DIR="$(dirname "$OUTPUT_TSV")"
RUN_DIR="${OUTPUT_DIR}/upimapi_native"
RESOURCES_DIR="$(dirname "$SWISSPROT_RESOURCE")"
RESULTS_TSV="${RUN_DIR}/UPIMAPI_results.tsv"
INFO_TSV="${RUN_DIR}/uniprotinfo.tsv"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')]  $*" | tee -a "$LOG_FILE"
}

{
    echo "==========================================="
    echo "UniProt UPIMAPI Search"
    echo "==========================================="
    echo "Timestamp: $(date)"
    echo "Host: $(hostname)"
    echo "User: $(whoami)"
    echo ""
    echo "Input:    $INPUT_FAA"
    echo "Output:   $OUTPUT_TSV"
    echo "Resource: $SWISSPROT_RESOURCE"
    echo "DB dir:   $RESOURCES_DIR"
    echo "Threads:  $THREADS"
    echo ""
    echo "==========================================="
    echo "Script: run_upimapi.sh"
    echo "==========================================="
    cat "$0"
    echo ""
    echo "==========================================="
    echo "Execution Output"
    echo "==========================================="
    echo ""
} > "$LOG_FILE"

log "UniProt UPIMAPI Search"
log "==========================================="
log "Input proteins:      $INPUT_FAA"
log "Output file:         $OUTPUT_TSV"
log "Swiss-Prot resource: $SWISSPROT_RESOURCE"
log "UPIMAPI resource dir:${RESOURCES_DIR}"
log "Threads:             $THREADS"
log ""

PROTEIN_COUNT=$(grep -c '^>' "$INPUT_FAA" || echo "0")
log "Total proteins in input: $PROTEIN_COUNT"
log ""

mkdir -p "$OUTPUT_DIR"
mkdir -p "$RUN_DIR"

log "Starting UPIMAPI..."
log "  Database: swissprot"
log "  Resources dir: $RESOURCES_DIR"
log "  Mode: more_sensitive"
log "  E-value threshold: 1e-20"
log "  Max target seqs: 1 (best hit only)"
log ""

START_TIME=$(date +%s)

upimapi \
    -i "$INPUT_FAA" \
    -o "$RUN_DIR" \
    -ot "$INFO_TSV" \
    -rd "$RESOURCES_DIR" \
    -db swissprot \
    --skip-db-check \
    -t "$THREADS" \
    --evalue 1e-20 \
    --max-target-seqs 1 \
    --diamond-mode more_sensitive 2>&1 | tee -a "$LOG_FILE"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

log ""
log "UPIMAPI run completed in ${ELAPSED}s"
log ""

if [ -f "$RESULTS_TSV" ]; then
    cp "$RESULTS_TSV" "$OUTPUT_TSV"
elif [ -f "$INFO_TSV" ]; then
    cp "$INFO_TSV" "$OUTPUT_TSV"
else
    log "ERROR: UPIMAPI did not create an annotation table"
    exit 1
fi

if [ -f "$OUTPUT_TSV" ]; then
    if [ -s "$OUTPUT_TSV" ]; then
        HIT_COUNT=$(tail -n +2 "$OUTPUT_TSV" | wc -l | tr -d ' ')
        PROTEINS_WITH_HITS=$(tail -n +2 "$OUTPUT_TSV" | cut -f1 | sort -u | wc -l | tr -d ' ')
    else
        HIT_COUNT="0"
        PROTEINS_WITH_HITS="0"
    fi

    if [ "$PROTEIN_COUNT" -gt 0 ]; then
        COVERAGE=$(awk "BEGIN {printf \"%.1f\", 100 * $PROTEINS_WITH_HITS / $PROTEIN_COUNT}")
    else
        COVERAGE="0.0"
    fi

    log "Results:"
    log "  Total annotated rows:    $HIT_COUNT"
    log "  Proteins with hits:      $PROTEINS_WITH_HITS / $PROTEIN_COUNT (${COVERAGE}%)"
    log "  Native UPIMAPI output:   $RUN_DIR"
    log ""
else
    log "ERROR: Output file not created!"
    exit 1
fi

log "==========================================="
log "✓ UPIMAPI search complete!"
log "==========================================="
log "Output file: $OUTPUT_TSV"
log "Log file:    $LOG_FILE"
log "==========================================="