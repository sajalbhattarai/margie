#!/usr/bin/env bash
# run_kofamscan.sh
# Runs KofamScan inside container with mounted database
# This script runs inside the container

set -euo pipefail

# Default parameters
INPUT_FASTA="/container/input/proteins.faa"
OUTPUT_FILE="/container/output/kofamscan_output.txt"
DB_DIR="/container/db"
THREADS=8
FORMAT="detail"

usage() {
    cat <<EOF
Usage: run_kofamscan.sh [OPTIONS]

Optional:
  --input FILE      Input protein FASTA file (default: /container/input/proteins.faa)
  --output FILE     Output annotation file (default: /container/output/kofamscan_output.txt)
  --db DIR          Database directory (default: /container/db)
  --threads NUM     Number of CPU threads (default: 4)
  --format FORMAT   Output format: detail|mapper|detail-tsv (default: detail)
  -h, --help        Show this help message

Example (run inside container):
  /container/scripts/run_kofamscan.sh --input /container/input/test.faa --threads 8
EOF
    exit 0
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input) INPUT_FASTA="$2"; shift 2 ;;
        --output) OUTPUT_FILE="$2"; shift 2 ;;
        --db) DB_DIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --format) FORMAT="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# Validate inputs
if [ ! -f "$INPUT_FASTA" ]; then
    echo "ERROR: Input file not found: $INPUT_FASTA"
    exit 1
fi

if [ ! -d "$DB_DIR" ]; then
    echo "ERROR: Database directory not found: $DB_DIR"
    exit 1
fi

# Check database files (support DB root at /container/db or /container/db/kegg)
if [ -d "${DB_DIR}/profiles" ] && [ -f "${DB_DIR}/ko_list" ]; then
    :
elif [ -d "${DB_DIR}/kegg/profiles" ] && [ -f "${DB_DIR}/kegg/ko_list" ]; then
    :
else
    echo "ERROR: KOfam DB files not found under ${DB_DIR}"
    echo "Expected either:"
    echo "  - ${DB_DIR}/profiles and ${DB_DIR}/ko_list"
    echo "  - ${DB_DIR}/kegg/profiles and ${DB_DIR}/kegg/ko_list"
    echo "Database may not be properly set up"
    exit 1
fi

# Generate config.yml dynamically
CONFIG_FILE="/tmp/kegg_config_$$.yml"
/container/scripts/generate_config.sh "$DB_DIR" "$CONFIG_FILE" "$THREADS"

# Create output directory
mkdir -p "$(dirname "$OUTPUT_FILE")"

# Setup logging
OUTPUT_BASE_DIR="$(dirname "$OUTPUT_FILE")"
LOG_FILE="${OUTPUT_BASE_DIR}/kegg_kofamscan.log"

# Function to log with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')]  $*" | tee -a "$LOG_FILE"
}

# Start logging
echo "=========================================" | tee "$LOG_FILE"
echo "KofamScan KEGG Annotation" | tee -a "$LOG_FILE"
echo "=========================================" | tee -a "$LOG_FILE"
log "Timestamp: $(date)"
log "Host: $(hostname)"
log "User: $(whoami)"
log ""
log "Input:    $INPUT_FASTA"
log "Output:   $OUTPUT_FILE"
log "Database: $DB_DIR"
log "Threads:  $THREADS"
log "Format:   $FORMAT"
log ""
log "========================================="
log "Script: run_kofamscan.sh"
log "========================================="
cat /container/scripts/run_kofamscan.sh >> "$LOG_FILE"
echo "" >> "$LOG_FILE"
log "========================================="
log ""

# Count input sequences
SEQ_COUNT=$(grep -c "^>" "$INPUT_FASTA" || echo 0)
log "Total proteins to search: $SEQ_COUNT"

# Build KofamScan command
KOFAMSCAN_CMD="/opt/kofam_scan/exec_annotation"
TMP_DIR="/tmp/kofamscan_$$"
mkdir -p "$TMP_DIR"

ARGS=(
    -o "$OUTPUT_FILE"
    --config="$CONFIG_FILE"
    --cpu="$THREADS"
    --tmp-dir="$TMP_DIR"
    --format="$FORMAT"
    "$INPUT_FASTA"
)

# Run KofamScan
log "Running KofamScan..."
log "Command: $KOFAMSCAN_CMD ${ARGS[*]}"
log ""
START_TIME=$(date +%s)
"$KOFAMSCAN_CMD" "${ARGS[@]}" 2>&1 | tee -a "$LOG_FILE"
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

# Cleanup
rm -rf "$TMP_DIR" "$CONFIG_FILE"

# Count results
if [ -f "$OUTPUT_FILE" ]; then
    if [ "$FORMAT" = "detail" ]; then
        HIT_COUNT=$(grep -c "^\*" "$OUTPUT_FILE" || echo 0)
        BELOW_THRESHOLD=$(grep -c "^[^*#]" "$OUTPUT_FILE" | tail -1 || echo 0)
        TOTAL_COUNT=$((HIT_COUNT + BELOW_THRESHOLD))
        
        log ""
        log "========================================"
        log "✓ KofamScan annotation complete!"
        log "========================================"
        log "Time:              ${ELAPSED}s"
        log "Above threshold:   ${HIT_COUNT} hits"
        log "Below threshold:   ${BELOW_THRESHOLD} hits"
        log "Total hits:        ${TOTAL_COUNT}"
        log "Output:            $OUTPUT_FILE"
        log "Log:               $LOG_FILE"
        log "========================================"
    else
        RESULT_COUNT=$(grep -c "^[^#]" "$OUTPUT_FILE" || echo 0)
        log ""
        log "========================================"
        log "✓ KofamScan annotation complete!"
        log "========================================"
        log "Time:    ${ELAPSED}s"
        log "Results: ${RESULT_COUNT} lines"
        log "Output:  $OUTPUT_FILE"
        log "Log:     $LOG_FILE"
        log "========================================"
    fi
else
    log "ERROR: Output file was not created!"
    exit 1
fi
