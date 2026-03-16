#!/bin/bash
# Run NCBI BLASTP against TCDB database for official-style transporter annotation

set -euo pipefail

# Input parameters (passed as environment variables or command line)
INPUT_FAA="${INPUT_FAA:-/container/input/proteins.faa}"
OUTPUT_DIR="${OUTPUT_DIR:-/container/output}"
THREADS="${THREADS:-8}"
EVALUE="${EVALUE:-1e-5}"

# Database prefix
BLAST_DB_PREFIX="/container/db/tcdb_blast"

# Output files
NATIVE_DIR="${OUTPUT_DIR}/native"
mkdir -p "$NATIVE_DIR"

BLAST_RAW_OUTPUT="${NATIVE_DIR}/tcdb_hits_raw.tsv"
BLAST_OUTPUT="${NATIVE_DIR}/tcdb_hits.tsv"
LOG_FILE="${OUTPUT_DIR}/tcdb_blast.log"

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')]  $*" | tee -a "$LOG_FILE"
}

# Write script header to log file
{
    echo "=============================================================================="
    echo "TCDB BLASTP Search Log"
    echo "=============================================================================="
    echo ""
    echo "Run Details:"
    echo "  Timestamp:         $(date '+%Y-%m-%d %H:%M:%S')"
    echo "  Input FASTA:       ${INPUT_FAA}"
    echo "  Output Directory:  ${OUTPUT_DIR}"
    echo "  Database Prefix:   ${BLAST_DB_PREFIX}"
    echo "  E-value:           ${EVALUE}"
    echo "  Threads:           ${THREADS}"
    echo ""
    echo "------------------------------------------------------------------------------"
    echo "Embedded Script: run_blast.sh"
    echo "------------------------------------------------------------------------------"
    cat "$0"
    echo ""
    echo "=============================================================================="
    echo "BLAST Execution Log"
    echo "=============================================================================="
    echo ""
} > "$LOG_FILE"

log "Starting TCDB BLASTP search..."
log ""

# Check if BLAST database exists
if [ ! -f "${BLAST_DB_PREFIX}.pin" ] && [ ! -f "${BLAST_DB_PREFIX}.phr" ]; then
    FASTA_CANDIDATE=""
    if [ -f "/container/db/tcdb.fasta" ]; then
        FASTA_CANDIDATE="/container/db/tcdb.fasta"
    elif [ -f "/container/db/tcdb.fa" ]; then
        FASTA_CANDIDATE="/container/db/tcdb.fa"
    fi

    if [ -n "$FASTA_CANDIDATE" ] && command -v makeblastdb >/dev/null 2>&1; then
        log "TCDB BLAST DB missing; building index from $FASTA_CANDIDATE"
        makeblastdb -in "$FASTA_CANDIDATE" -dbtype prot -out "$BLAST_DB_PREFIX" >> "$LOG_FILE" 2>&1 || true
    fi

    if [ ! -f "${BLAST_DB_PREFIX}.pin" ] && [ ! -f "${BLAST_DB_PREFIX}.phr" ]; then
        log "ERROR: TCDB BLAST database not found: ${BLAST_DB_PREFIX}"
        exit 1
    fi
fi

# Count input proteins
TOTAL_PROTEINS=$(grep -c "^>" "$INPUT_FAA" || echo "0")
log "Input validation:"
log "  Total proteins:    $TOTAL_PROTEINS"
log ""

# Record start time
START_TIME=$(date +%s)

log "Running BLASTP..."
log "  E-value threshold:    $EVALUE"
log "  Coverage thresholds:  50% (query) and 50% (subject)"
log "  Max target seqs:      1 (best hit per protein)"
log ""

# Output schema matches consolidate_tcdb.py expectations after post-processing:
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp scovhsp stitle
blastp \
    -query "$INPUT_FAA" \
    -db "$BLAST_DB_PREFIX" \
    -out "$BLAST_RAW_OUTPUT" \
    -num_threads "$THREADS" \
    -evalue "$EVALUE" \
    -max_target_seqs 1 \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp stitle'

EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
    log "ERROR: BLASTP search failed with exit code $EXIT_CODE"
    exit $EXIT_CODE
fi

# Add derived subject coverage and apply the same coverage thresholds used by the DIAMOND path.
# Subject coverage proxy is aligned_length / subject_length * 100.
awk -F '\t' -v OFS='\t' '
    {
        qcov = ($15 == "" ? 0 : $15)
        slen = ($14 == "" ? 0 : $14)
        scov = (slen > 0 ? (100.0 * $4 / slen) : 0)
        if (qcov >= 50 && scov >= 50) {
            print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,qcov,scov,$16
        }
    }
' "$BLAST_RAW_OUTPUT" > "$BLAST_OUTPUT"

# Calculate execution time
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
log ""
log "✓ BLASTP search complete"
log "  Execution time:  ${ELAPSED}s"
log ""

# Add header expected by consolidation step
log "Adding TSV header..."
TEMP_OUTPUT="${BLAST_OUTPUT}.tmp"
mv "$BLAST_OUTPUT" "$TEMP_OUTPUT"
echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tqcovhsp\tscovhsp\tstitle" > "$BLAST_OUTPUT"
cat "$TEMP_OUTPUT" >> "$BLAST_OUTPUT"
rm -f "$TEMP_OUTPUT" "$BLAST_RAW_OUTPUT"
log ""

# Count hits and calculate statistics
TOTAL_HITS=$(($(wc -l < "$BLAST_OUTPUT") - 1))
PROTEINS_WITH_HITS=$(cut -f1 "$BLAST_OUTPUT" | tail -n +2 | sort -u | wc -l | tr -d ' ')

log "Results summary:"
log "  Total hits:              $TOTAL_HITS"
log "  Proteins with hits:      $PROTEINS_WITH_HITS / $TOTAL_PROTEINS"

if [ "$TOTAL_PROTEINS" -gt 0 ]; then
    COVERAGE=$(awk "BEGIN {printf \"%.1f\", 100*$PROTEINS_WITH_HITS/$TOTAL_PROTEINS}")
    log "  Coverage:                ${COVERAGE}%"
fi

log ""
log "Output files:"
log "  TCDB hits:   $BLAST_OUTPUT"
log "  BLAST log:   $LOG_FILE"
log ""
log "TCDB BLASTP search completed successfully!"
