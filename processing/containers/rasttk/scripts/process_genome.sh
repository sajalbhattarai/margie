#!/bin/bash

# Simple wrapper for the complex pipeline script
# Usage: run_genome <name> <fasta> [scientific_name]

if [ $# -lt 2 ]; then
    echo "Usage: run_genome <genome_name> <fasta_path> [scientific_name]"
    echo "Example: run_genome my_strain data/contigs.fasta \"Backteroides thetaiotaomicron\""
    exit 1
fi

GENOME_NAME=$1
INPUT_FASTA=$2
SCIENTIFIC_NAME=${3:-$GENOME_NAME}

# Default output directory in the container
OUTPUT_DIR="output/rasttk/$GENOME_NAME"

echo "Starting annotation for $GENOME_NAME..."
echo "Input: $INPUT_FASTA"
echo "Output: $OUTPUT_DIR"

/scripts/run_rasttk_incremental.sh \
    "$GENOME_NAME" \
    "$INPUT_FASTA" \
    "$OUTPUT_DIR" \
    "$SCIENTIFIC_NAME" \
    11 \
    "Bacteria"

# Run finalization (md files)
python3 /scripts/finalize_structure.py "$OUTPUT_DIR" "$GENOME_NAME" "$SCIENTIFIC_NAME"

echo ""
echo "Done! Results are in $OUTPUT_DIR"
