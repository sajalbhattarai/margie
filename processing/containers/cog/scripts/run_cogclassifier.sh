#!/bin/bash
# Run official COGclassifier workflow for COG functional classification

set -euo pipefail

INPUT_FAA="${INPUT_FAA:-/container/input/proteins.faa}"
OUTPUT_DIR="${OUTPUT_DIR:-/container/output}"
THREADS="${THREADS:-4}"
EVALUE="${EVALUE:-0.01}"
DB_DIR="${DB_DIR:-/container/db}"
COG_QUIET="${COG_QUIET:-0}"

NATIVE_DIR="${OUTPUT_DIR}/native"
mkdir -p "$NATIVE_DIR"

LOG_FILE="${OUTPUT_DIR}/cogclassifier.log"

{
    echo "=========================================="
    echo "COGclassifier Annotation"
    echo "=========================================="
    echo "Date: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Input: $INPUT_FAA"
    echo "Output: $OUTPUT_DIR"
    echo "Database: $DB_DIR"
    echo "Threads: $THREADS"
    echo "E-value: $EVALUE"
    echo ""
} | tee "$LOG_FILE"

if [ ! -f "$INPUT_FAA" ]; then
    echo "ERROR: Input protein file not found: $INPUT_FAA" | tee -a "$LOG_FILE"
    exit 1
fi

TOTAL_PROTEINS=$(grep -c "^>" "$INPUT_FAA" || echo "0")
echo "Total proteins to classify: $TOTAL_PROTEINS" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

echo "Running COGclassifier..." | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

COG_ARGS=(
    -i "$INPUT_FAA"
    -o "$NATIVE_DIR"
    -d "$DB_DIR"
    -t "$THREADS"
    -e "$EVALUE"
)

# Default is live output; set COG_QUIET=1 to suppress verbose tool chatter.
if [[ "$COG_QUIET" == "1" ]]; then
    COG_ARGS+=(--quiet)
fi

if command -v stdbuf >/dev/null 2>&1; then
    COG_EXIT=0
    stdbuf -oL -eL COGclassifier "${COG_ARGS[@]}" 2>&1 | tee -a "$LOG_FILE" || COG_EXIT=$?
else
    COG_EXIT=0
    COGclassifier "${COG_ARGS[@]}" 2>&1 | tee -a "$LOG_FILE" || COG_EXIT=$?
fi

if [[ "${COG_EXIT:-0}" -ne 0 ]]; then
    if grep -q "KeyError:" "$LOG_FILE"; then
        echo "" | tee -a "$LOG_FILE"
        echo "[warn] COGclassifier crashed on unmapped CDD ID (known upstream edge case)." | tee -a "$LOG_FILE"
        echo "[run]  Falling back to tolerant COG pipeline..." | tee -a "$LOG_FILE"
        python3 /container/scripts/run_cogclassifier_safe.py \
            --input "$INPUT_FAA" \
            --outdir "$NATIVE_DIR" \
            --db "$DB_DIR" \
            --threads "$THREADS" \
            --evalue "$EVALUE" 2>&1 | tee -a "$LOG_FILE"
        COG_EXIT=0
    else
        echo "ERROR: COGclassifier failed. See log: $LOG_FILE" | tee -a "$LOG_FILE"
        exit "$COG_EXIT"
    fi
fi

echo "" | tee -a "$LOG_FILE"
echo "✓ COGclassifier complete" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

if [ -f "${NATIVE_DIR}/cog_classify.tsv" ]; then
    TOTAL_HITS=$(($(wc -l < "${NATIVE_DIR}/cog_classify.tsv") - 1))
    echo "Total COG hits: $TOTAL_HITS" | tee -a "$LOG_FILE"

    if [ "$TOTAL_PROTEINS" -gt 0 ]; then
        PERCENT=$(awk "BEGIN {printf \"%.1f\", 100*$TOTAL_HITS/$TOTAL_PROTEINS}")
        echo "Coverage: ${PERCENT}%" | tee -a "$LOG_FILE"
    fi
else
    echo "WARNING: COG classification file not found" | tee -a "$LOG_FILE"
fi

# Explicitly generate publication-ready charts when plot CLIs are available.
if [ -f "${NATIVE_DIR}/cog_count.tsv" ]; then
    if command -v plot_cog_count_barchart >/dev/null 2>&1; then
        echo "Generating COG bar charts..." | tee -a "$LOG_FILE"
        plot_cog_count_barchart \
            -i "${NATIVE_DIR}/cog_count.tsv" \
            -o "${OUTPUT_DIR}/cog_count_barchart.png" 2>&1 | tee -a "$LOG_FILE" || true
        plot_cog_count_barchart \
            -i "${NATIVE_DIR}/cog_count.tsv" \
            -o "${OUTPUT_DIR}/cog_count_barchart.html" 2>&1 | tee -a "$LOG_FILE" || true
    else
        echo "WARNING: plot_cog_count_barchart CLI not found; skipping chart generation" | tee -a "$LOG_FILE"
    fi

    if command -v plot_cog_count_piechart >/dev/null 2>&1; then
        echo "Generating COG pie charts..." | tee -a "$LOG_FILE"
        plot_cog_count_piechart \
            -i "${NATIVE_DIR}/cog_count.tsv" \
            -o "${OUTPUT_DIR}/cog_count_piechart.png" 2>&1 | tee -a "$LOG_FILE" || true
        plot_cog_count_piechart \
            -i "${NATIVE_DIR}/cog_count.tsv" \
            -o "${OUTPUT_DIR}/cog_count_piechart.html" 2>&1 | tee -a "$LOG_FILE" || true
    else
        echo "WARNING: plot_cog_count_piechart CLI not found; skipping chart generation" | tee -a "$LOG_FILE"
    fi
fi

echo "" | tee -a "$LOG_FILE"
echo "Output files:" | tee -a "$LOG_FILE"
echo "  RPS-BLAST hits:      ${NATIVE_DIR}/rpsblast.tsv" | tee -a "$LOG_FILE"
echo "  COG classification:  ${NATIVE_DIR}/cog_classify.tsv" | tee -a "$LOG_FILE"
echo "  COG counts:          ${NATIVE_DIR}/cog_count.tsv" | tee -a "$LOG_FILE"
echo "  COG bar chart PNG:   ${OUTPUT_DIR}/cog_count_barchart.png" | tee -a "$LOG_FILE"
echo "  COG pie chart PNG:   ${OUTPUT_DIR}/cog_count_piechart.png" | tee -a "$LOG_FILE"
echo "  COG bar chart HTML:  ${OUTPUT_DIR}/cog_count_barchart.html" | tee -a "$LOG_FILE"
echo "  COG pie chart HTML:  ${OUTPUT_DIR}/cog_count_piechart.html" | tee -a "$LOG_FILE"
echo "  Log file:            $LOG_FILE" | tee -a "$LOG_FILE"
