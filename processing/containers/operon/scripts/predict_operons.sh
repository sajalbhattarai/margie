#!/bin/bash
#
# Operon Prediction Wrapper Script
# Uses UniOP (Universal Operon Predictor) for machine learning-based operon prediction
# Supports local Docker images or Apptainer/Singularity SIFs built by setup_containers.sh
#
# Usage: predict_operons.sh --gff <gff_file> --output <output_dir>
#

set -euo pipefail

# Default values
GFF_FILE=""
OUTPUT_DIR=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --gff)
            GFF_FILE="$2"
            shift 2
            ;;
        --output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Validate arguments
if [ -z "$GFF_FILE" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 --gff <gff_file> --output <output_dir>"
    exit 1
fi

# Verify GFF file exists
if [ ! -f "$GFF_FILE" ]; then
    echo "ERROR: GFF file not found: $GFF_FILE"
    exit 1
fi

# Get absolute paths
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
WORKSPACE_ROOT="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
GFF_FILE_ABS="$(cd "$(dirname "$GFF_FILE")" && pwd)/$(basename "$GFF_FILE")"

# Create parent directory for OUTPUT_DIR if it doesn't exist (needed to get absolute path)
mkdir -p "$(dirname "$OUTPUT_DIR")"
OUTPUT_DIR_ABS="$(cd "$(dirname "$OUTPUT_DIR")" && pwd)/$(basename "$OUTPUT_DIR")"

# Determine genome name and protein FASTA from GFF path
# GFF path format: output/rasttk/<genome_name>/gene_calls/<genome_name>.gff
GENOME_NAME=$(basename "$(dirname "$(dirname "$GFF_FILE")")")
PROTEIN_FASTA="$(dirname "$GFF_FILE")/${GENOME_NAME}.faa"

if [ ! -f "$PROTEIN_FASTA" ]; then
    echo "ERROR: Protein FASTA not found: $PROTEIN_FASTA"
    exit 1
fi

# Container and script configuration
LOCAL_CONTAINER="operon-annotation:1.0"
FALLBACK_CONTAINER="ghcr.io/sajalbhattarai/for_oprn:latest"
LOCAL_SIF="$WORKSPACE_ROOT/processing/containers/operon/operon.sif"
CONTAINER_NAME=""
CONTAINER_RUNTIME=""
RUNTIME_BIN=""
CONVERT_SCRIPT="$SCRIPT_DIR/convert_rast_to_uniop.py"
REFORMAT_SCRIPT="$SCRIPT_DIR/reformat_uniop_output.py"
COMPREHENSIVE_SCRIPT="$SCRIPT_DIR/create_comprehensive_operon_table.py"

# Check runtime availability (local Docker image first, then fallback image, then Apptainer/Singularity SIF)
if command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1 && docker image inspect "$LOCAL_CONTAINER" &> /dev/null; then
    CONTAINER_NAME="$LOCAL_CONTAINER"
    CONTAINER_RUNTIME="docker"
elif command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1 && docker image inspect "$FALLBACK_CONTAINER" &> /dev/null; then
    CONTAINER_NAME="$FALLBACK_CONTAINER"
    CONTAINER_RUNTIME="docker"
elif command -v apptainer >/dev/null 2>&1 && [ -f "$LOCAL_SIF" ]; then
    CONTAINER_NAME="$LOCAL_SIF"
    CONTAINER_RUNTIME="apptainer"
    RUNTIME_BIN="apptainer"
elif command -v singularity >/dev/null 2>&1 && [ -f "$LOCAL_SIF" ]; then
    CONTAINER_NAME="$LOCAL_SIF"
    CONTAINER_RUNTIME="apptainer"
    RUNTIME_BIN="singularity"
else
    echo "ERROR: Operon container not found"
    echo "  Expected local image: $LOCAL_CONTAINER"
    echo "  Or fallback image: $FALLBACK_CONTAINER"
    echo "  Or local SIF: $LOCAL_SIF"
    echo "  Build all containers with: ./setup_containers.sh --docker"
    echo "  Or on HPC: ./setup_containers.sh --apptainer"
    exit 1
fi

# Create output directories
mkdir -p "$OUTPUT_DIR_ABS"
mkdir -p "$OUTPUT_DIR_ABS/native"
mkdir -p "$OUTPUT_DIR_ABS/logs"
mkdir -p "$OUTPUT_DIR_ABS/scripts"

# Setup logging
LOG_FILE="$OUTPUT_DIR_ABS/logs/operon_prediction.log"
START_TIME=$(date +%s)

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')]  $*" | tee -a "$LOG_FILE"
}

# Start log file with embedded script
{
    echo "=========================================="
    echo "Operon Prediction Log"
    echo "=========================================="
    echo "Embedded Script: predict_operons.sh"
    echo "=========================================="
    cat "$0"
    echo ""
    echo "=========================================="
    echo "Execution Log"
    echo "=========================================="
} > "$LOG_FILE"

# Save script to output for reproducibility
cp "$0" "$OUTPUT_DIR_ABS/scripts/"
cp "$CONVERT_SCRIPT" "$OUTPUT_DIR_ABS/scripts/" 2>/dev/null || true
cp "$REFORMAT_SCRIPT" "$OUTPUT_DIR_ABS/scripts/" 2>/dev/null || true
cp "$COMPREHENSIVE_SCRIPT" "$OUTPUT_DIR_ABS/scripts/" 2>/dev/null || true

log "======================================================================"
log "Operon Prediction - UniOP"
log "======================================================================"
log ""
log "Configuration:"
log "  Genome:          $GENOME_NAME"
log "  Input GFF:       $GFF_FILE"
log "  Input proteins:  $PROTEIN_FASTA"
log "  Output dir:      $OUTPUT_DIR_ABS"
log "  Runtime:         $CONTAINER_RUNTIME"
log "  Container:       $CONTAINER_NAME"
log ""

# Count proteins
TOTAL_PROTEINS=$(grep -c "^>" "$PROTEIN_FASTA" || echo "0")
log "  Total proteins:  $TOTAL_PROTEINS"
log ""

# Step 1: Convert RAST format to Prodigal-style format for UniOP
PRODIGAL_FAA="$OUTPUT_DIR_ABS/${GENOME_NAME}_prodigal_style.faa"

log "[Step 1/4] Converting RAST format to UniOP format..."
log ""

python3 "$CONVERT_SCRIPT" \
    "$GFF_FILE_ABS" \
    "$PROTEIN_FASTA" \
    "$PRODIGAL_FAA" \
    2>&1 | tee -a "$LOG_FILE"

if [ ${PIPESTATUS[0]} -ne 0 ]; then
    log "ERROR: Format conversion failed!"
    exit 1
fi

log ""
log "✓ Format conversion complete"

# Step 2: Run UniOP (containerized)
log ""
log "[Step 2/4] Running operon prediction..."
log ""

if [ "$CONTAINER_RUNTIME" = "docker" ]; then
    docker run --rm \
        -v "$PRODIGAL_FAA:/input/proteins.faa:ro" \
        -v "$OUTPUT_DIR_ABS/native:/output" \
        "$CONTAINER_NAME" \
        -i /input/proteins.faa -o /output --genome "$GENOME_NAME" \
        2>&1 | tee -a "$LOG_FILE"
else
    "$RUNTIME_BIN" exec \
        --bind "$PRODIGAL_FAA:/input/proteins.faa" \
        --bind "$OUTPUT_DIR_ABS/native:/output" \
        "$CONTAINER_NAME" \
        operon_exec -i /input/proteins.faa -o /output --genome "$GENOME_NAME" \
        2>&1 | tee -a "$LOG_FILE"
fi

if [ ${PIPESTATUS[0]} -ne 0 ]; then
    log "ERROR: UniOP prediction failed!"
    exit 1
fi

log ""
log "✓ Operon prediction complete"

FINAL_TSV="$OUTPUT_DIR_ABS/operons.tsv"
if [ -f "$OUTPUT_DIR_ABS/native/$GENOME_NAME/operons.tsv" ] && [ ! -f "$FINAL_TSV" ]; then
    cp "$OUTPUT_DIR_ABS/native/$GENOME_NAME/operons.tsv" "$FINAL_TSV"
fi

if [ -f "$FINAL_TSV" ]; then
    log ""
    log "Detected modern operon_exec output: $FINAL_TSV"
    LINE_COUNT=$(wc -l < "$FINAL_TSV" | tr -d ' ')
    GENE_COUNT=$((LINE_COUNT - 1))
    OPERON_COUNT=$(tail -n +2 "$FINAL_TSV" | awk -F'\t' '$2 != "NA" {print $2}' | sort -u | wc -l | tr -d ' ')
    GENES_IN_OPS=$(tail -n +2 "$FINAL_TSV" | awk -F'\t' '$2 != "NA"' | wc -l | tr -d ' ')
    if [ $TOTAL_PROTEINS -gt 0 ]; then
        COVERAGE=$(awk "BEGIN {printf \"%.1f\", 100.0 * $GENES_IN_OPS / $TOTAL_PROTEINS}")
    else
        COVERAGE="0.0"
    fi

    log ""
    log "=========================================================================="
    log "                   Operon Prediction Complete!                      "
    log "=========================================================================="
    log ""
    log "Summary:"
    log "  • Total genes analyzed: $TOTAL_PROTEINS"
    log "  • Genes in operons: $GENES_IN_OPS ($COVERAGE%)"
    log "  • Total operons: $OPERON_COUNT"
    log ""
    log "Main output file:"
    log "  • $FINAL_TSV"
    log ""
    log "Columns: feature_id, operon_id, operon_size, operon_position, member_genes, ..."
    log ""
    exit 0
fi

# Legacy UniOP path (older images)
UNIOP_PRED="$OUTPUT_DIR_ABS/native/uniop.pred"
UNIOP_OPERON="$OUTPUT_DIR_ABS/native/uniop.operon"

if [ ! -f "$UNIOP_PRED" ] || [ ! -f "$UNIOP_OPERON" ]; then
    log "ERROR: Operon output files not found (neither modern operons.tsv nor legacy uniop.* files)"
    exit 1
fi

# Step 3: Reformat UniOP output
log ""
log "[Step 3/4] Reformatting operon predictions..."
log ""

python3 "$REFORMAT_SCRIPT" \
    "$GENOME_NAME" \
    "$PROTEIN_FASTA" \
    "$UNIOP_PRED" \
    "$UNIOP_OPERON" \
    "$OUTPUT_DIR_ABS" \
    2>&1 | tee -a "$LOG_FILE"

if [ ${PIPESTATUS[0]} -ne 0 ]; then
    log "ERROR: Reformatting failed!"
    exit 1
fi

# Step 4: Create comprehensive operon annotation table (legacy UniOP mode)
log ""
log "[Step 4/4] Creating comprehensive operon annotation table..."
log ""

# First create simple table for native_processed
REFORMATTED_OPERON="$OUTPUT_DIR_ABS/processed_uniop_${GENOME_NAME}.operon"
mkdir -p "$OUTPUT_DIR_ABS/native_processed"
SIMPLE_TSV="$OUTPUT_DIR_ABS/native_processed/operons_simple.tsv"

# Parse the reformatted operon file to create simple TSV
python3 - <<EOF 2>&1 | tee -a "$LOG_FILE"
import sys

operon_file = "$REFORMATTED_OPERON"
output_file = "$SIMPLE_TSV"

# Read operon assignments
operon_data = []
with open(operon_file, 'r') as f:
    for line in f:
        if line.startswith('#') or line.startswith('Organism'):
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 4:
            organism = parts[0]
            operon_id = parts[1]
            num_genes = int(parts[2])
            gene_ids_str = parts[3]
            
            # Parse gene IDs
            gene_ids = [gid.strip() for gid in gene_ids_str.split(',')]
            
            # Create entry for each gene
            for position, gene_id in enumerate(gene_ids, start=1):
                operon_data.append({
                    'feature_id': gene_id,
                    'operon_id': operon_id,
                    'operon_size': num_genes,
                    'operon_position': position
                })

# Write TSV
with open(output_file, 'w') as f:
    # Write header
    f.write('feature_id\toperon_id\toperon_size\toperon_position\n')
    
    # Write data
    for entry in operon_data:
        f.write(f"{entry['feature_id']}\t{entry['operon_id']}\t{entry['operon_size']}\t{entry['operon_position']}\n")

# Count operons
num_operons = len(set(entry['operon_id'] for entry in operon_data))
num_genes = len(operon_data)

print(f"  Simple table saved to native_processed/")
EOF

if [ ${PIPESTATUS[0]} -ne 0 ]; then
    log "ERROR: Simple TSV creation failed!"
    exit 1
fi

# Now create comprehensive table
FINAL_TSV="$OUTPUT_DIR_ABS/operons.tsv"

log ""
log "  Creating comprehensive table with prediction scores and distances..."

python3 "$COMPREHENSIVE_SCRIPT" \
    "$GENOME_NAME" \
    "$PRODIGAL_FAA" \
    "$UNIOP_PRED" \
    "$REFORMATTED_OPERON" \
    "$FINAL_TSV" \
    2>&1 | tee -a "$LOG_FILE"

if [ ${PIPESTATUS[0]} -ne 0 ]; then
    log "ERROR: Comprehensive table creation failed!"
    exit 1
fi

# Verify output
if [ -f "$FINAL_TSV" ]; then
    LINE_COUNT=$(wc -l < "$FINAL_TSV" | tr -d ' ')
    GENE_COUNT=$((LINE_COUNT - 1))
    
    # Count operons (unique operon_ids excluding NA)
    OPERON_COUNT=$(tail -n +2 "$FINAL_TSV" | awk -F'\t' '$2 != "NA" {print $2}' | sort -u | wc -l | tr -d ' ')
    GENES_IN_OPS=$(tail -n +2 "$FINAL_TSV" | awk -F'\t' '$2 != "NA"' | wc -l | tr -d ' ')
    
    if [ $TOTAL_PROTEINS -gt 0 ]; then
        COVERAGE=$(awk "BEGIN {printf \"%.1f\", 100.0 * $GENES_IN_OPS / $TOTAL_PROTEINS}")
    else
        COVERAGE="0.0"
    fi
    
    log ""
    log "=========================================================================="
    log "                   Operon Prediction Complete!                      "
    log "=========================================================================="
    log ""
    log "Comprehensive operon table created"
    log ""
    log "Summary:"
    log "  • Total genes analyzed: $TOTAL_PROTEINS"
    log "  • Genes in operons: $GENES_IN_OPS ($COVERAGE%)"
    log "  • Total operons: $OPERON_COUNT"
    log ""
    log "Main output file:"
    log "  • $FINAL_TSV"
    log ""
    log "    Columns: feature_id, operon_id, operon_size, operon_position,"
    log "             member_genes, upstream_gene, downstream_gene,"
    log "             upstream_operon_score, downstream_operon_score,"
    log "             intergenic_distance_upstream, intergenic_distance_downstream,"
    log "             contig, start, end, strand"
    log ""
    log "Native UniOP files (for reference):"
    log "  • $OUTPUT_DIR_ABS/native/uniop.pred (pairwise predictions)"
    log "  • $OUTPUT_DIR_ABS/native/uniop.operon (operon clusters)"
    log ""
    log "Processed files (for PUL detection):"
    log "  • $OUTPUT_DIR_ABS/native_processed/operons_simple.tsv (basic operon assignments)"
    log "  • $OUTPUT_DIR_ABS/processed_uniop_${GENOME_NAME}.pred (reformatted predictions)"
    log "  • $OUTPUT_DIR_ABS/processed_uniop_${GENOME_NAME}.operon (reformatted clusters)"
    log ""
else
    log "ERROR: Output file not created"
    exit 1
fi
