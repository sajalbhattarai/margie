#!/bin/bash
set -euo pipefail

# Pfam HMMER Annotation Script
# Runs hmmscan against Pfam-A database for protein family annotation
# Input: Protein FASTA file
# Output: HMMER domtblout format

# Default values
INPUT_FAA=""
OUTPUT_DIR="/container/output"
DATABASE="/container/db/Pfam-A.hmm"
THREADS=8
EVALUE_CUTOFF=0.01
DOMEVAL_CUTOFF=0.01

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
        --evalue)
            EVALUE_CUTOFF="$2"
            shift 2
            ;;
        --domevalue)
            DOMEVAL_CUTOFF="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Validate inputs
if [[ -z "$INPUT_FAA" ]]; then
    echo "ERROR: --input required"
    exit 1
fi

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
LOG_FILE="${OUTPUT_DIR}/pfam_hmmscan.log"

# Function to log with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')]  $*" | tee -a "$LOG_FILE"
}

# Start logging
echo "=========================================" | tee "$LOG_FILE"
echo "Pfam HMMER hmmscan Annotation" | tee -a "$LOG_FILE"
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
log "Input proteins: $PROTEIN_COUNT"
log ""

# Check if database is pressed
if [[ ! -f "${DATABASE}.h3i" ]]; then
    log "WARNING: HMM database not pressed. Performance will be slower."
    log "Press the database with: hmmpress $DATABASE"
fi

# Run hmmscan with curated gathering thresholds and aggressive acceleration
log "Running hmmscan (using curated --cut_ga thresholds with aggressive speedup)..."
log "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
log "Command: hmmscan --cpu $THREADS --cut_ga --nobias --F1 0.005 --F2 1e-5 --domtblout $OUTPUT_DIR/pfam.domtblout --noali -o $OUTPUT_DIR/pfam.out $DATABASE $INPUT_FAA"
log "Optimizations: --nobias (no composition bias), --F1 0.005, --F2 1e-5 (aggressive filtering, ~40-50% faster)"
log ""

START_TIME=$(date +%s)
hmmscan \
    --cpu "$THREADS" \
    --cut_ga \
    --nobias \
    --F1 0.005 \
    --F2 1e-5 \
    --domtblout "$OUTPUT_DIR/pfam.domtblout" \
    --noali \
    -o "$OUTPUT_DIR/pfam.out" \
    "$DATABASE" \
    "$INPUT_FAA" 2>&1 | tee -a "$LOG_FILE"
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

log ""
log "End time: $(date '+%Y-%m-%d %H:%M:%S')"
log "Elapsed time: ${ELAPSED}s"
log ""

# Check output
if [[ ! -f "$OUTPUT_DIR/pfam.domtblout" ]]; then
    log "ERROR: HMMER output not created"
    exit 1
fi

# Post-process domain hits
log "Post-processing domain hits..."
log ""

# a. Extract all hits (no filtering)
log "Extracting all domain hits..."
grep -v '^#' "$OUTPUT_DIR/pfam.domtblout" > "$OUTPUT_DIR/pfam.all.domains"
ALL_COUNT=$(wc -l < "$OUTPUT_DIR/pfam.all.domains" | tr -d ' ')
log "  All domain hits: $ALL_COUNT"

# b. Extract best hit per protein
log "Extracting best hit per protein..."
grep -v '^#' "$OUTPUT_DIR/pfam.domtblout" | \
    awk '{q=$4; sc=$14; if(!(q in best) || sc>bestsc[q]){best[q]=$0; bestsc[q]=sc}} END{for(q in best) print best[q]}' \
    > "$OUTPUT_DIR/pfam.best1.domains"
BEST1_COUNT=$(wc -l < "$OUTPUT_DIR/pfam.best1.domains" | tr -d ' ')
log "  Best-hit domains: $BEST1_COUNT"

# c. Extract non-overlapping domains (for final annotation)
log "Extracting non-overlapping domains..."
grep -v '^#' "$OUTPUT_DIR/pfam.domtblout" | \
    sort -k4,4 -k14,14nr | \
    awk '{q=$4; s=$20; e=$21; if(q!=pq){delete a; delete b; n=0; pq=q} ok=1; for(i=1;i<=n;i++){ if(!(e<a[i] || s>b[i])){ok=0; break} } if(ok){n++; a[n]=s; b[n]=e; print}}' \
    > "$OUTPUT_DIR/pfam.nonoverlap.domains"
NONOVERLAP_COUNT=$(wc -l < "$OUTPUT_DIR/pfam.nonoverlap.domains" | tr -d ' ')
log "  Non-overlapping domains: $NONOVERLAP_COUNT"

# Count unique proteins
UNIQUE_PROTEINS=$(cut -f4 "$OUTPUT_DIR/pfam.nonoverlap.domains" | sort -u | wc -l | tr -d ' ')

log ""
log "========================================="
log "✓ HMMER hmmscan complete!"
log "========================================="
log "Time:                    ${ELAPSED}s"
log "All domain hits:         $ALL_COUNT"
log "Best-hit domains:        $BEST1_COUNT"
log "Non-overlapping domains: $NONOVERLAP_COUNT"
log "Proteins with domains:   $UNIQUE_PROTEINS / $PROTEIN_COUNT"
log "Output files:"
log "  - Raw output:           $OUTPUT_DIR/pfam.out"
log "  - Domain table:         $OUTPUT_DIR/pfam.domtblout"
log "  - All domains:          $OUTPUT_DIR/pfam.all.domains"
log "  - Best-hit domains:     $OUTPUT_DIR/pfam.best1.domains"
log "  - Non-overlapping:      $OUTPUT_DIR/pfam.nonoverlap.domains"
log "  - Log file:             $LOG_FILE"
log "========================================="
exit 0