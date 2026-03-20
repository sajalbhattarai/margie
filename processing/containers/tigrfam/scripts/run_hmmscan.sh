#!/bin/bash
# TIGRfam HMMER annotation script with TC/NC threshold decision logic
# Mirrors Pfam workflow with strictness assessment

set -euo pipefail

# Default parameters
INPUT_FAA="/container/input/proteins.faa"
OUTPUT_DIR="/container/output"
DATABASE="/container/db/TIGRFAMs_15.0_HMM.LIB"
THREADS=8

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input)
            INPUT_FAA="$2"
            shift 2
            ;;
        --output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --database)
            DATABASE="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Validate inputs
if [[ ! -f "$INPUT_FAA" ]]; then
    echo "ERROR: Input file not found: $INPUT_FAA"
    exit 1
fi

if [[ ! -f "$DATABASE" ]]; then
    echo "ERROR: Database not found: $DATABASE"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Extract genome name
GENOME_NAME=$(basename "$INPUT_FAA" .faa)

# Setup logging
LOG_FILE="${OUTPUT_DIR}/tigrfam_hmmscan.log"

# Function to log with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')]  $*" | tee -a "$LOG_FILE"
}

# Start logging
echo "=========================================" | tee "$LOG_FILE"
echo "TIGRfam HMMER hmmscan Annotation" | tee -a "$LOG_FILE"
echo "=========================================" | tee -a "$LOG_FILE"
log "Timestamp: $(date)"
log "Host: $(hostname)"
log "User: $(whoami)"
log ""
log "Input:    $INPUT_FAA"
log "Output:   $OUTPUT_DIR"
log "Database: $DATABASE"
log "Threads:  $THREADS"
log ""
log "========================================="
log "Script: run_hmmscan.sh"
log "========================================="
cat /container/scripts/run_hmmscan.sh >> "$LOG_FILE"
log ""
log "========================================="
log "Execution Output"
log "========================================="
log ""

# Count input sequences
PROTEIN_COUNT=$(grep -c '^>' "$INPUT_FAA" || echo "0")
log "Total proteins in proteome: $PROTEIN_COUNT"
log ""

# Check if database is pressed
if [[ ! -f "${DATABASE}.h3i" ]]; then
    log "WARNING: HMM database not pressed. Performance will be slower."
    log "Press the database with: hmmpress $DATABASE"
fi

#################################################
# Run hmmscan with --cut_tc (trusted cutoff)
#################################################
log "========================================="
log "Running hmmscan with --cut_tc"
log "========================================="
log "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
log "Command: hmmscan --cpu $THREADS --cut_tc --domtblout $OUTPUT_DIR/tigrfam.domtblout --noali -o $OUTPUT_DIR/tigrfam.out $DATABASE $INPUT_FAA"
log ""

START_TIME=$(date +%s)
hmmscan \
    --cpu "$THREADS" \
    --cut_tc \
    --domtblout "$OUTPUT_DIR/tigrfam.domtblout" \
    --noali \
    -o "$OUTPUT_DIR/tigrfam.out" \
    "$DATABASE" \
    "$INPUT_FAA" 2>&1 | tee -a "$LOG_FILE"
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

log ""
log "Scan completed in ${ELAPSED}s"

#################################################
# Post-processing domain hits
#################################################
log ""
log "========================================="
log "Post-processing domain hits"
log "========================================="

# a. Extract all hits (no filtering)
log "Extracting all domain hits..."
grep -v '^#' "$OUTPUT_DIR/tigrfam.domtblout" > "$OUTPUT_DIR/tigrfam.all.domains"
ALL_COUNT=$(wc -l < "$OUTPUT_DIR/tigrfam.all.domains" | tr -d ' ')
log "  All domain hits: $ALL_COUNT"

# b. Extract best hit per protein
log "Extracting best hit per protein..."
grep -v '^#' "$OUTPUT_DIR/tigrfam.domtblout" | \
    awk '{q=$4; sc=$14; if(!(q in best) || sc>bestsc[q]){best[q]=$0; bestsc[q]=sc}} END{for(q in best) print best[q]}' \
    > "$OUTPUT_DIR/tigrfam.best1.domains"
BEST1_COUNT=$(wc -l < "$OUTPUT_DIR/tigrfam.best1.domains" | tr -d ' ')
log "  Best-hit domains: $BEST1_COUNT"

# c. Extract non-overlapping domains (for final annotation)
log "Extracting non-overlapping domains..."
grep -v '^#' "$OUTPUT_DIR/tigrfam.domtblout" | \
    sort -k4,4 -k14,14nr | \
    awk '{q=$4; s=$20; e=$21; if(q!=pq){delete a; delete b; n=0; pq=q} ok=1; for(i=1;i<=n;i++){ if(!(e<a[i] || s>b[i])){ok=0; break} } if(ok){n++; a[n]=s; b[n]=e; print}}' \
    > "$OUTPUT_DIR/tigrfam.nonoverlap.domains"
NONOVERLAP_COUNT=$(wc -l < "$OUTPUT_DIR/tigrfam.nonoverlap.domains" | tr -d ' ')
log "  Non-overlapping domains: $NONOVERLAP_COUNT"

# Count unique proteins
UNIQUE_PROTEINS=$(cut -f4 "$OUTPUT_DIR/tigrfam.all.domains" | sort -u | wc -l | tr -d ' ')

log ""
log "========================================="
log "✓ HMMER hmmscan complete!"
log "========================================="
log "Scan time:               ${ELAPSED}s"
log ""
log "Results (--cut_tc trusted cutoff):"
log "  All domain hits:       $ALL_COUNT"
log "  Best-hit domains:      $BEST1_COUNT"
log "  Non-overlapping:       $NONOVERLAP_COUNT"
log "  Proteins with domains: $UNIQUE_PROTEINS / $PROTEIN_COUNT"
log ""
log "Output files:"
log "  - Raw output:              $OUTPUT_DIR/tigrfam.out"
log "  - Domain table:            $OUTPUT_DIR/tigrfam.domtblout"
log "  - All domains:             $OUTPUT_DIR/tigrfam.all.domains"
log "  - Best-hit domains:        $OUTPUT_DIR/tigrfam.best1.domains"
log "  - Non-overlapping:         $OUTPUT_DIR/tigrfam.nonoverlap.domains"
log "  - Log file:                $LOG_FILE"
log "========================================="

# Generate structured TSV for downstream consolidation and interpro mapping
awk -v OFS='\t' '
  BEGIN { print "feature_id", "TIGRFAM_id", "TIGRFAM_description" }
  !/^#/ { print $4, $1, $1 }
' "$OUTPUT_DIR/tigrfam.nonoverlap.domains" > "$OUTPUT_DIR/tigrfam.tsv" || true
log "Created tigrfam.tsv: $(tail -n +2 "$OUTPUT_DIR/tigrfam.tsv" 2>/dev/null | wc -l | tr -d ' ') entries"

exit 0