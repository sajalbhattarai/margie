#!/bin/bash
#
# Batch Genome Annotation Script
# Processes multiple genomes sequentially with BV-BRC RASTtk
#
# Usage: ./run_batch_genomes.sh <batch_file>
#   Or:  run_batch_genomes <batch_file>
#
# Batch file format (TSV with 3 columns):
#   genome_name<TAB>fasta_file<TAB>scientific_name
#
# Example batch file (genomes.txt):
#   ref_colik12	input/ref_colik12/ref_colik12.fna	Escherichia coli K12
#   ref_longum	input/ref_longum/ref_longum.fna	Bifidobacterium longum
#   ref_subtilis168	input/ref_subtilis168/ref_subtilis168.fna	Bacillus subtilis 168
#   ref_theta	input/ref_theta/ref_theta.fna	Bacteroides thetaiotaomicron
#

set -e

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${GREEN}============================================${NC}"
echo -e "${GREEN}  BV-BRC Batch Genome Annotation           ${NC}"
echo -e "${GREEN}============================================${NC}"
echo ""

# Check if batch file provided
if [ $# -eq 0 ]; then
    echo -e "${RED}Error: No batch file specified${NC}"
    echo ""
    echo "Usage: run_batch_genomes <batch_file>"
    echo ""
    echo "Batch file format (TSV):"
    echo "  genome_name<TAB>fasta_file<TAB>scientific_name"
    echo ""
    echo "Example:"
    echo "  ref_colik12	input/ref_colik12/ref_colik12.fna	Escherichia coli K12"
    echo "  ref_longum	input/ref_longum/ref_longum.fna	Bifidobacterium longum"
    echo ""
    exit 1
fi

BATCH_FILE="$1"

# Check if batch file exists
if [ ! -f "$BATCH_FILE" ]; then
    echo -e "${RED}Error: Batch file not found: $BATCH_FILE${NC}"
    exit 1
fi

# Check if logged in
if ! p3-whoami &>/dev/null; then
    echo -e "${RED}Error: Not logged in to BV-BRC${NC}"
    echo ""
    echo "Please login first:"
    echo "  p3-login <username>"
    echo ""
    exit 1
fi

USER=$(p3-whoami 2>/dev/null | grep -o 'user .*' | awk '{print $2}')
echo -e "${GREEN}✓ Logged in as: $USER${NC}"
echo ""

# Count total genomes
TOTAL=$(grep -v '^#' "$BATCH_FILE" | grep -v '^[[:space:]]*$' | wc -l | tr -d ' ')
echo -e "${BLUE}Found $TOTAL genome(s) to process${NC}"
echo ""

# Process each genome
CURRENT=0
SUCCESS=0
FAILED=0
START_TIME=$(date +%s)

while IFS=$'\t' read -r genome_name fasta_file sci_name || [ -n "$genome_name" ]; do
    # Skip comments and empty lines
    [[ "$genome_name" =~ ^#.*$ ]] && continue
    [[ -z "$genome_name" ]] && continue
    
    CURRENT=$((CURRENT + 1))
    
    echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${YELLOW}[$CURRENT/$TOTAL] Processing: $genome_name${NC}"
    echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo "  FASTA: $fasta_file"
    echo "  Name: $sci_name"
    echo ""
    
    # Check if FASTA file exists
    if [ ! -f "$fasta_file" ]; then
        echo -e "${RED}✗ Error: FASTA file not found: $fasta_file${NC}"
        echo -e "${RED}  Skipping this genome...${NC}"
        echo ""
        FAILED=$((FAILED + 1))
        continue
    fi
    
    # Check if already processed
    if [ -f "output/rasttk/$genome_name/rast.tsv" ]; then
        echo -e "${YELLOW}⚠ Warning: output/rasttk/$genome_name/rast.tsv already exists${NC}"
        read -p "  Overwrite? (y/N): " -n 1 -r
        echo ""
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "  Skipping..."
            echo ""
            continue
        fi
        echo "  Removing existing output..."
        rm -rf "output/rasttk/$genome_name"
    fi
    
    # Run annotation
    GENOME_START=$(date +%s)
    
    if /scripts/process_genome.sh "$genome_name" "$fasta_file" "$sci_name"; then
        GENOME_END=$(date +%s)
        GENOME_TIME=$((GENOME_END - GENOME_START))
        echo ""
        echo -e "${GREEN}✓ Success: $genome_name completed in ${GENOME_TIME}s${NC}"
        SUCCESS=$((SUCCESS + 1))
    else
        GENOME_END=$(date +%s)
        GENOME_TIME=$((GENOME_END - GENOME_START))
        echo ""
        echo -e "${RED}✗ Failed: $genome_name (took ${GENOME_TIME}s)${NC}"
        FAILED=$((FAILED + 1))
    fi
    
    echo ""
    
done < "$BATCH_FILE"

END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))
MINUTES=$((TOTAL_TIME / 60))
SECONDS=$((TOTAL_TIME % 60))

echo -e "${GREEN}============================================${NC}"
echo -e "${GREEN}  Batch Processing Complete                ${NC}"
echo -e "${GREEN}============================================${NC}"
echo ""
echo "Results:"
echo "  Total genomes: $TOTAL"
echo -e "  ${GREEN}Successful: $SUCCESS${NC}"
if [ $FAILED -gt 0 ]; then
    echo -e "  ${RED}Failed: $FAILED${NC}"
fi
echo "  Total time: ${MINUTES}m ${SECONDS}s"
echo ""

if [ $SUCCESS -gt 0 ]; then
    echo "Output files:"
    for output in output/rasttk/*/rast.tsv; do
        if [ -f "$output" ]; then
            genome=$(echo "$output" | sed 's|output/rasttk/||; s|/rast.tsv||')
            features=$(tail -n +2 "$output" | wc -l | tr -d ' ')
            echo "  ✓ $genome: $features features"
        fi
    done
    echo ""
fi

if [ $FAILED -gt 0 ]; then
    echo -e "${YELLOW}Some genomes failed. Check the output above for errors.${NC}"
    echo ""
fi

echo "Done!"
