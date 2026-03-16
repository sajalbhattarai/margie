#!/bin/bash
# Run Diamond blastp against EggNOG database for orthology annotation

set -e

# Input parameters
INPUT_FAA="${INPUT_FAA:-/container/input/proteins.faa}"
OUTPUT_DIR="${OUTPUT_DIR:-/container/output}"
THREADS="${THREADS:-8}"
EVALUE="${EVALUE:-1e-5}"

# Database
DIAMOND_DB="/container/db/bacteria.dmnd"

# Output files
NATIVE_DIR="${OUTPUT_DIR}/native"
LOGS_DIR="${OUTPUT_DIR}/logs"
mkdir -p "$NATIVE_DIR" "$LOGS_DIR"

DIAMOND_OUTPUT="${NATIVE_DIR}/eggnog_hits.tsv"
LOG_FILE="${LOGS_DIR}/eggnog_diamond.log"

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')]  $*" | tee -a "$LOG_FILE"
}

# Start logging - embed script for reproducibility
{
    echo "=========================================="
    echo "eggNOG Diamond Search Log"
    echo "=========================================="
    echo "Embedded Script: run_diamond.sh"
    echo "=========================================="
    cat "$0"
    echo ""
    echo "=========================================="
    echo "Diamond Execution Log"
    echo "=========================================="
} > "$LOG_FILE"

START_TIME=$(date +%s)

log "EggNOG Diamond Search"
log "Input:    $INPUT_FAA"
log "Database: $DIAMOND_DB"
log "E-value:  $EVALUE"
log "Threads:  $THREADS"
log "Output:   $DIAMOND_OUTPUT"
log ""

# Check if database exists
if [ ! -f "$DIAMOND_DB" ]; then
    log "ERROR: EggNOG Diamond database not found: $DIAMOND_DB"
    exit 1
fi

# Count input proteins
TOTAL_PROTEINS=$(grep -c "^>" "$INPUT_FAA" || echo "0")
log "Total proteins to search: $TOTAL_PROTEINS"
log ""

log "Running Diamond blastp (mid-sensitive mode with 500MB block size)..."
log "  E-value threshold:    $EVALUE"
log "  Max target seqs:      1 (best hit per protein)"
log "  Block size:           0.5 GB (memory-optimized)"
log "  Mode:                 mid-sensitive"
log ""

# Run Diamond blastp with memory-optimized parameters
# Using --mid-sensitive and --block-size 0.5 to handle the 9.2GB eggNOG database
# Format includes qlen/slen/qcovhsp/scovhsp for proper validity computation
diamond blastp \
    --query "$INPUT_FAA" \
    --db "$DIAMOND_DB" \
    --out "$DIAMOND_OUTPUT" \
    --threads "$THREADS" \
    --evalue "$EVALUE" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp scovhsp \
    --max-target-seqs 1 \
    --mid-sensitive \
    --block-size 0.5

EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
    log "ERROR: Diamond search failed with exit code $EXIT_CODE"
    exit $EXIT_CODE
fi

# Calculate execution time
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
log ""
log "✓ Diamond search complete"
log "  Execution time:  ${ELAPSED}s"
log ""

# Add header
log "Adding TSV header..."
TEMP_OUTPUT="${DIAMOND_OUTPUT}.tmp"
mv "$DIAMOND_OUTPUT" "$TEMP_OUTPUT"
echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tqcovhsp\tscovhsp" > "$DIAMOND_OUTPUT"
cat "$TEMP_OUTPUT" >> "$DIAMOND_OUTPUT"
rm "$TEMP_OUTPUT"
log ""

# Count hits and calculate statistics
TOTAL_HITS=$(($(wc -l < "$DIAMOND_OUTPUT") - 1))
PROTEINS_WITH_HITS=$(cut -f1 "$DIAMOND_OUTPUT" | tail -n +2 | sort -u | wc -l | tr -d ' ')

log "Results summary:"
log "  Total hits:              $TOTAL_HITS"
log "  Proteins with hits:      $PROTEINS_WITH_HITS / $TOTAL_PROTEINS"

if [ "$TOTAL_PROTEINS" -gt 0 ]; then
    COVERAGE=$(awk "BEGIN {printf \"%.1f\", 100*$PROTEINS_WITH_HITS/$TOTAL_PROTEINS}")
    log "  Coverage:                ${COVERAGE}%"
fi

log ""
log "Output files:"
log "  eggNOG hits:  $DIAMOND_OUTPUT"
log "  Diamond log:  $LOG_FILE"
log ""
log "eggNOG Diamond search completed successfully!"
