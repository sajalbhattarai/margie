#!/bin/bash
# Detect Polysaccharide Utilization Loci (PULs) from consolidated annotations
# 
# PULs are genomic regions encoding systems for degrading complex carbohydrates:
# - Multiple CAZymes clustered together (≥3 within 10kb)
# - SusC/SusD outer membrane transport systems
# - Co-localized with regulatory genes

set -euo pipefail

# Default values
CONSOLIDATED_FILE=""
OUTPUT_DIR=""
GENOME_NAME=""
GENERATE_HTML=false

usage() {
    cat <<EOF
Usage: $0 [OPTIONS]

Required:
  -i, --input FILE        Consolidated annotations TSV file
  -o, --output DIR        Output directory
  -n, --name NAME         Genome name (for output files)

Optional:
  --html                  Generate HTML visualization report
  -h, --help             Show this help message

Example:
  $0 -i output/consolidated/test/consolidated_annotations.tsv \\
     -o output/puls -n test --html

What this detects:
  - CAZyme clustering: ≥3 CAZymes within 10kb
  - SusC/SusD systems: Outer membrane carbohydrate uptake
  - PUL boundaries: Gene clusters with carbohydrate degradation machinery

Output:
  output_dir/
    ├── genome_name/
    │   ├── puls.tsv                    # Detected PUL regions
    │   ├── puls_report.html            # Visualization (if --html)
    │   └── logs/
    │       └── pul_detection.log       # Execution log
EOF
    exit 0
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input) CONSOLIDATED_FILE="$2"; shift 2 ;;
        -o|--output) OUTPUT_DIR="$2"; shift 2 ;;
        -n|--name) GENOME_NAME="$2"; shift 2 ;;
        --html) GENERATE_HTML=true; shift ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# Validate required arguments
if [ -z "$CONSOLIDATED_FILE" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$GENOME_NAME" ]; then
    echo "ERROR: Missing required arguments"
    usage
fi

# Validate input file
if [ ! -f "$CONSOLIDATED_FILE" ]; then
    echo "ERROR: Input file not found: $CONSOLIDATED_FILE"
    exit 1
fi

# Get absolute paths
CONSOLIDATED_FILE=$(cd "$(dirname "$CONSOLIDATED_FILE")" && pwd)/$(basename "$CONSOLIDATED_FILE")
OUTPUT_DIR=$(cd "$(dirname "$OUTPUT_DIR")" 2>/dev/null && pwd)/$(basename "$OUTPUT_DIR") || OUTPUT_DIR=$(pwd)/$OUTPUT_DIR

# Create output directory structure
GENOME_OUTPUT="${OUTPUT_DIR}/${GENOME_NAME}"
mkdir -p "${GENOME_OUTPUT}/logs"
mkdir -p "${GENOME_OUTPUT}/scripts"

# Setup logging
LOG_FILE="${GENOME_OUTPUT}/logs/pul_detection.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')]  $*" | tee -a "$LOG_FILE"
}

# Embed script in log for reproducibility
{
    echo "=========================================="
    echo "PUL Detection Log"
    echo "Embedded Script: run_pul_detection.sh"
    cat "$0"
    echo "=========================================="
    echo "Execution Log"
    echo "=========================================="
} > "$LOG_FILE"

# Copy script
cp "$0" "${GENOME_OUTPUT}/scripts/run_pul_detection.sh"

log "======================================================================"
log "Polysaccharide Utilization Loci (PUL) Detection"
log "======================================================================"
log ""
log "Configuration:"
log "  Genome:      $GENOME_NAME"
log "  Input:       $CONSOLIDATED_FILE"
log "  Output:      $GENOME_OUTPUT"
log "  HTML report: $GENERATE_HTML"
log ""

# Count input proteins
TOTAL_PROTEINS=$(tail -n +2 "$CONSOLIDATED_FILE" | wc -l | tr -d ' ')
log "  Total proteins: $TOTAL_PROTEINS"
log ""

# Count proteins with CAZymes (DBCAN_family_count > 0)
CAZYME_COUNT=$(awk -F'\t' '
    NR==1 {
        for(i=1; i<=NF; i++) {
            if($i == "DBCAN_family_count") col=i
        }
    }
    NR>1 && col>0 && $col > 0 {
        count++
    }
    END {print count+0}
' "$CONSOLIDATED_FILE")

if [ -z "$CAZYME_COUNT" ] || [ "$CAZYME_COUNT" -eq 0 ]; then
    CAZYME_COUNT=0
fi

log "  Proteins with CAZymes: $CAZYME_COUNT"
log ""

if [ "$CAZYME_COUNT" -lt 3 ]; then
    log "WARNING: Only $CAZYME_COUNT CAZyme(s) found"
    log "         PUL detection requires ≥3 CAZymes clustered within 10kb"
    log "         Results may be limited"
    log ""
fi

# Get workspace root
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
WORKSPACE_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# Detection script
DETECTION_SCRIPT="${SCRIPT_DIR}/detect_pul.py"

if [ ! -f "$DETECTION_SCRIPT" ]; then
    log "ERROR: PUL detection script not found: $DETECTION_SCRIPT"
    exit 1
fi

log "Validating detection script..."
log "✓ Script found: detect_pul.py"
log ""

# Run PUL detection
log "Running PUL detection..."
log ""
log "Detection criteria:"
log "  - CAZyme clustering: ≥3 CAZymes within 10kb"
log "  - SusC/SusD systems: TonB-dependent transporters + lipoproteins"
log "  - Gene co-localization: Sequential genes on same contig"
log ""

# Build command
DETECTION_CMD=(
    python3 "$DETECTION_SCRIPT"
    --annotations "$CONSOLIDATED_FILE"
    --output "${GENOME_OUTPUT}/puls.tsv"
)

if [ "$GENERATE_HTML" = true ]; then
    DETECTION_CMD+=(--html "${GENOME_OUTPUT}/puls_report.html")
fi

"${DETECTION_CMD[@]}" 2>&1 | tee -a "$LOG_FILE"

if [ ${PIPESTATUS[0]} -ne 0 ]; then
    log ""
    log "ERROR: PUL detection failed"
    exit 1
fi

log ""
log "✓ PUL detection complete"
log ""

# Verify and summarize output
if [ -f "${GENOME_OUTPUT}/puls.tsv" ]; then
    LINE_COUNT=$(wc -l < "${GENOME_OUTPUT}/puls.tsv" | tr -d ' ')
    PUL_COUNT=$((LINE_COUNT - 1))
    
    log "======================================================================"
    log "PUL Detection Complete!"
    log "======================================================================"
    log ""
    log "Summary:"
    log "  Total proteins analyzed:   $TOTAL_PROTEINS"
    log "  Proteins with CAZymes:     $CAZYME_COUNT"
    log "  PULs detected:             $PUL_COUNT"
    log ""
    
    if [ $PUL_COUNT -gt 0 ]; then
        log "Detected PUL details:"
        log ""
        
        # Parse PUL summary from TSV
        tail -n +2 "${GENOME_OUTPUT}/puls.tsv" | while IFS=$'\t' read -r pul_id contig start end num_genes num_cazymes has_suscd genes; do
            log "  PUL: $pul_id"
            log "    Location:   ${contig}:${start}-${end} ($(( (end - start) / 1000 ))kb)"
            log "    Genes:      $num_genes total, $num_cazymes CAZymes"
            log "    SusC/SusD:  $has_suscd"
            log ""
        done
    else
        log "  No PULs detected"
        log ""
        if [ "$CAZYME_COUNT" -lt 3 ]; then
            log "  Reason: Insufficient CAZymes (need ≥3 clustered within 10kb)"
        else
            log "  CAZymes may be dispersed across genome (not clustered)"
        fi
        log ""
    fi
    
    log "Output files:"
    log "  PUL regions: ${GENOME_OUTPUT}/puls.tsv"
    
    if [ "$GENERATE_HTML" = true ] && [ -f "${GENOME_OUTPUT}/puls_report.html" ]; then
        log "  HTML report: ${GENOME_OUTPUT}/puls_report.html"
    fi
    
    log "  Execution log: ${LOG_FILE}"
    log ""
else
    log "ERROR: Output file not created"
    exit 1
fi

log "======================================================================"
log ""
