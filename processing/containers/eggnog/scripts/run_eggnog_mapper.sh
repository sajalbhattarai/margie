#!/bin/bash
# Run official eggnog-mapper for orthology annotation

set -euo pipefail

# Input parameters
INPUT_FAA="${INPUT_FAA:-/container/input/proteins.faa}"
OUTPUT_DIR="${OUTPUT_DIR:-/container/output}"
THREADS="${THREADS:-8}"
ORGANISM_NAME="${ORGANISM_NAME:-genome}"
EGGNOG_TAX_SCOPE="${EGGNOG_TAX_SCOPE:-2,2157}"
EGGNOG_TAX_SCOPE_MODE="${EGGNOG_TAX_SCOPE_MODE:-inner_narrowest}"
EGGNOG_TARGET_TAXA="${EGGNOG_TARGET_TAXA:-}"

# Database directory
DB_DIR="/container/db"

# Output files
NATIVE_DIR="${OUTPUT_DIR}/native"
mkdir -p "$NATIVE_DIR"

# Log file
LOG_FILE="${OUTPUT_DIR}/eggnog_mapper.log"

# Function to log with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

# Start logging
echo "========================================" | tee "$LOG_FILE"
echo "EggNOG Mapper Annotation" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
log "Timestamp: $(date)"
log "Host: $(hostname)"
log "User: $(whoami)"
log ""
log "Input FASTA: $INPUT_FAA"
log "Database directory: $DB_DIR"
log "Output directory: $OUTPUT_DIR"
log "Threads: $THREADS"
log "Organism name: $ORGANISM_NAME"
log "Tax scope: $EGGNOG_TAX_SCOPE (mode: $EGGNOG_TAX_SCOPE_MODE)"
if [ -n "$EGGNOG_TARGET_TAXA" ]; then
    log "Target taxa: $EGGNOG_TARGET_TAXA"
fi
log ""

# Embed this script in the log file
log "========================================" 
log "Script: run_eggnog_mapper.sh"
log "========================================" 
cat /container/scripts/run_eggnog_mapper.sh >> "$LOG_FILE"
echo "" >> "$LOG_FILE"
log "========================================" 
log ""

# Check if input exists
if [ ! -f "$INPUT_FAA" ]; then
    log "ERROR: Input FASTA not found: $INPUT_FAA"
    exit 1
fi

if [ ! -f "$DB_DIR/eggnog.db" ] || [ ! -f "$DB_DIR/eggnog_proteins.dmnd" ]; then
    log "ERROR: EggNOG database is incomplete in $DB_DIR"
    log "       Required files: eggnog.db and eggnog_proteins.dmnd"
    exit 1
fi

# Count input proteins
TOTAL_PROTEINS=$(grep -c "^>" "$INPUT_FAA" || echo "0")
log "Total proteins to annotate: $TOTAL_PROTEINS"
log ""

# Check if eggnog-mapper is installed
if ! command -v emapper.py &> /dev/null; then
    log "ERROR: emapper.py not found. Is eggnog-mapper installed?"
    exit 1
fi

log "Running eggnog-mapper..."
log "Command: emapper.py -i $INPUT_FAA --output ${ORGANISM_NAME} --output_dir $NATIVE_DIR --data_dir $DB_DIR --cpu $THREADS -m diamond --tax_scope $EGGNOG_TAX_SCOPE --tax_scope_mode $EGGNOG_TAX_SCOPE_MODE --override"
log ""

SCOPE_ARGS=(--tax_scope "$EGGNOG_TAX_SCOPE" --tax_scope_mode "$EGGNOG_TAX_SCOPE_MODE")
if [ -n "$EGGNOG_TARGET_TAXA" ]; then
    SCOPE_ARGS+=(--target_taxa "$EGGNOG_TARGET_TAXA")
fi

# Run eggnog-mapper
# -i: input FASTA file
# --output: output file prefix  
# --output_dir: output directory
# --data_dir: database directory
# --cpu: number of threads
# -m: search mode (diamond is default and recommended)
# --override: overwrite existing files
emapper.py \
    -i "$INPUT_FAA" \
    --output "${ORGANISM_NAME}" \
    --output_dir "$NATIVE_DIR" \
    --data_dir "$DB_DIR" \
    --cpu "$THREADS" \
    -m diamond \
    "${SCOPE_ARGS[@]}" \
    --override 2>&1 | tee -a "$LOG_FILE"

log ""
log "✓ EggNOG mapper complete"

# Check output files
ANNOTATIONS_FILE="${NATIVE_DIR}/${ORGANISM_NAME}.emapper.annotations"
if [ -f "$ANNOTATIONS_FILE" ]; then
    # Count annotations (excluding header lines starting with #)
    TOTAL_ANNOTATIONS=$(grep -v "^#" "$ANNOTATIONS_FILE" | wc -l | tr -d ' ')
    log "Total annotations: $TOTAL_ANNOTATIONS"
    
    if [ "$TOTAL_PROTEINS" -gt 0 ]; then
        PERCENT=$(awk "BEGIN {printf \"%.1f\", 100*$TOTAL_ANNOTATIONS/$TOTAL_PROTEINS}")
        log "Coverage: ${PERCENT}%"
    fi
    
    log ""
    log "Annotations saved to: $ANNOTATIONS_FILE"
else
    log "ERROR: Annotations file not found: $ANNOTATIONS_FILE"
    exit 1
fi

log ""
log "========================================"
log "EggNOG Mapper Complete"
log "========================================"
log "Output directory: $NATIVE_DIR"
log "Log file: $LOG_FILE"
log ""

# Generate consolidated eggnog.tsv for downstream aggregation
if python3 /container/scripts/consolidate_eggnog.py \
        "$ANNOTATIONS_FILE" "$INPUT_FAA" "$ORGANISM_NAME" "$OUTPUT_DIR/eggnog.tsv" \
        >> "$LOG_FILE" 2>&1; then
    log "Created eggnog.tsv: $(tail -n +2 "$OUTPUT_DIR/eggnog.tsv" 2>/dev/null | wc -l | tr -d ' ') entries"
else
    log "Warning: eggnog TSV consolidation failed; raw annotations still available at $ANNOTATIONS_FILE"
fi
