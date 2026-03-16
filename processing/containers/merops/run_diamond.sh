#!/bin/bash
# MEROPS Diamond blastp script

set -euo pipefail

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <genome_name> <input_faa> <db_dmnd> <output_dir> <threads>"
    exit 1
fi

GENOME_NAME="$1"
INPUT_FAA="$2"
DB_DMND="$3"
OUTPUT_DIR="$4"
THREADS="$5"

# Setup logging
LOG_FILE="${OUTPUT_DIR}/../merops_diamond.log"
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')]  $*" | tee -a "$LOG_FILE"
}

# Log header with embedded script
{
    echo "==========================================="
    echo "MEROPS Peptidase Annotation - Diamond"
    echo "==========================================="
    echo "Timestamp: $(date)"
    echo "Host: $(hostname)"
    echo "User: $(whoami)"
    echo ""
    echo "Input:    $INPUT_FAA"
    echo "Output:   ${OUTPUT_DIR}/merops_hits.tsv"
    echo "Database: $DB_DMND"
    echo "Threads:  $THREADS"
    echo ""
    echo "==========================================="
    echo "Script: run_diamond.sh"
    echo "==========================================="
    cat "$0"
    echo "==========================================="
    echo ""
} > "$LOG_FILE"

log "MEROPS Peptidase Annotation - Diamond"
log "==========================================="
log "Genome:   $GENOME_NAME"
log "Input:    $INPUT_FAA"
log "Database: $DB_DMND"
log "Threads:  $THREADS"
log ""

PROTEIN_COUNT=$(grep -c "^>" "$INPUT_FAA" || true)
log "Total proteins in input: $PROTEIN_COUNT"
log ""

log "Starting DIAMOND blastp..."
log "  Mode: sensitive"
log "  E-value threshold: 1e-5"
log "  Query coverage: ≥50%"
log "  Subject coverage: ≥50%"
log "  Max target seqs: 1 (best hit only)"
log ""

START_TIME=$(date +%s)

DIAMOND_OUTPUT="$OUTPUT_DIR/merops_hits.tsv"

# Run Diamond blastp with best practice parameters
# Format includes qlen/slen/qcovhsp/scovhsp for proper validity computation
# Keep stitle for metadata parsing in consolidation
diamond blastp \
    --query "$INPUT_FAA" \
    --db "$DB_DMND" \
    --out "$DIAMOND_OUTPUT" \
    --threads "$THREADS" \
    --evalue 1e-5 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp scovhsp stitle \
    --max-target-seqs 1 \
    --sensitive \
    --query-cover 50 \
    --subject-cover 50 2>&1 | tee -a "$LOG_FILE"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

log ""
log "DIAMOND search completed in ${ELAPSED}s"
log ""

if [ -f "$DIAMOND_OUTPUT" ]; then
    HIT_COUNT=$(wc -l < "$DIAMOND_OUTPUT" | tr -d ' ')
    PROTEINS_WITH_HITS=$(awk '{print $1}' "$DIAMOND_OUTPUT" | sort -u | wc -l | tr -d ' ')
    
    if [ "$PROTEIN_COUNT" -gt 0 ]; then
        COVERAGE=$(awk "BEGIN {printf \"%.1f\", 100 * $PROTEINS_WITH_HITS / $PROTEIN_COUNT}")
    else
        COVERAGE="0.0"
    fi
    
    log "Results:"
    log "  Total hits:              $HIT_COUNT"
    log "  Proteins with hits:      $PROTEINS_WITH_HITS / $PROTEIN_COUNT (${COVERAGE}%)"
    log ""
    log "==========================================="
    log "Output file: $DIAMOND_OUTPUT"
    log "Log file:    $LOG_FILE"
    log "==========================================="
else
    log "Warning: No output generated"
fi
