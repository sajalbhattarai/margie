#!/bin/bash
set -euo pipefail

# dbCAN HMMER Annotation Script
# Runs hmmscan against dbCAN-HMMdb for CAZyme annotation
# Input: Protein FASTA file
# Output: HMMER domtblout format

# Default values
INPUT_FAA=""
OUTPUT_DIR="/container/output"
DATABASE="/container/db/dbCAN-HMMdb-V12.txt"
THREADS=8
EVALUE_CUTOFF=1e-15
COVERAGE_CUTOFF=0.35

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
        --coverage)
            COVERAGE_CUTOFF="$2"
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

echo "========================================"
echo "dbCAN CAZyme Annotation - HMMER"
echo "========================================"
echo "Genome: $GENOME_NAME"
echo "Input: $INPUT_FAA"
echo "Database: $DATABASE"
echo "Threads: $THREADS"
echo "E-value cutoff: $EVALUE_CUTOFF"
echo "Coverage cutoff: $COVERAGE_CUTOFF"
echo "Output: $OUTPUT_DIR"
echo "========================================"
echo ""

# Count input sequences
PROTEIN_COUNT=$(grep -c '^>' "$INPUT_FAA" || echo "0")
echo "Input proteins: $PROTEIN_COUNT"
echo ""

# Check if database is pressed
if [[ ! -f "${DATABASE}.h3i" ]]; then
    echo "WARNING: HMM database not pressed. This will be slow."
    echo "Press the database with: hmmpress $DATABASE"
fi

# Run hmmscan with E-value threshold (not all models have GA scores)
echo "Running hmmscan (E-value threshold: $EVALUE_CUTOFF)..."
echo "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

hmmscan \
    --cpu "$THREADS" \
    -E "$EVALUE_CUTOFF" \
    --domE "$EVALUE_CUTOFF" \
    --domtblout "$OUTPUT_DIR/dbcan_domains.tsv" \
    --tblout "$OUTPUT_DIR/dbcan_hits.tsv" \
    --notextw \
    -o "$OUTPUT_DIR/hmmscan_output.txt" \
    "$DATABASE" \
    "$INPUT_FAA"

echo ""
echo "End time: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Check output
if [[ ! -f "$OUTPUT_DIR/dbcan_domains.tsv" ]]; then
    echo "ERROR: HMMER output not created"
    exit 1
fi

# Count domain hits
DOMAIN_COUNT=$(grep -v '^#' "$OUTPUT_DIR/dbcan_domains.tsv" | wc -l | tr -d ' ')
UNIQUE_PROTEINS=$(grep -v '^#' "$OUTPUT_DIR/dbcan_domains.tsv" | cut -f1 | sort -u | wc -l | tr -d ' ')

echo "========================================"
echo "HMMER Scan Complete"
echo "========================================"
echo "Domain hits: $DOMAIN_COUNT"
echo "Proteins with domains: $UNIQUE_PROTEINS / $PROTEIN_COUNT"
echo "Output files:"
echo "  - Domain table: $OUTPUT_DIR/dbcan_domains.tsv"
echo "  - Full output: $OUTPUT_DIR/hmmscan_output.txt"
echo "========================================"
