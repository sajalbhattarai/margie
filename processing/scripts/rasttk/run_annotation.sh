#!/usr/bin/env bash
# RASTtk gene-calling + annotation wrapper for Margie
set -euo pipefail

usage() {
    cat <<EOF
Usage: $0 --input <fasta> --output <dir> [OPTIONS]

Required:
  --input, -i         Input genome FASTA
  --output, -o        Output directory

Optional:
  --name, -n          Genome name (default: basename of input)
  --scientific        Scientific name (default: Unknown organism)
  --genetic-code      Genetic code (default: 11)
  --domain            Domain (default: Bacteria)
  --threads, -t       Threads (default: 8)
EOF
    exit 1
}

SOURCE="${BASH_SOURCE[0]}"
while [[ -L "$SOURCE" ]]; do
    SRC_DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"
    SOURCE="$(readlink "$SOURCE")"
    [[ "$SOURCE" != /* ]] && SOURCE="$SRC_DIR/$SOURCE"
done
SCRIPT_DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
source "$ROOT_DIR/processing/scripts/lib/runtime_db_resolver.sh"

IMAGE_NAME="rasttk-annotation:1.0"
THREADS=8
SCIENTIFIC_NAME=""
GENETIC_CODE=11
DOMAIN="Bacteria"
INPUT_FASTA=""
OUTPUT_DIR=""
GENOME_NAME=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input|-i) INPUT_FASTA="$2"; shift 2 ;;
        --output|-o) OUTPUT_DIR="$2"; shift 2 ;;
        --name|-n) GENOME_NAME="$2"; shift 2 ;;
        --scientific) SCIENTIFIC_NAME="$2"; shift 2 ;;
        --genetic-code) GENETIC_CODE="$2"; shift 2 ;;
        --domain) DOMAIN="$2"; shift 2 ;;
        --threads|-t) THREADS="$2"; shift 2 ;;
        --help|-h) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

if [[ -z "$INPUT_FASTA" || -z "$OUTPUT_DIR" ]]; then
    usage
fi
if [[ ! -f "$INPUT_FASTA" ]]; then
    echo "Error: Input file not found: $INPUT_FASTA"
    exit 1
fi

INPUT_FASTA="$(cd "$(dirname "$INPUT_FASTA")" && pwd)/$(basename "$INPUT_FASTA")"
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"

if [[ -z "$GENOME_NAME" ]]; then
    GENOME_NAME="$(basename "$INPUT_FASTA")"
    GENOME_NAME="${GENOME_NAME%.fasta}"
    GENOME_NAME="${GENOME_NAME%.fna}"
    GENOME_NAME="${GENOME_NAME%.fa}"
fi
if [[ -z "$SCIENTIFIC_NAME" ]]; then
    SCIENTIFIC_NAME="Unknown organism"
fi

GENOME_OUTPUT="$OUTPUT_DIR/$GENOME_NAME"
mkdir -p "$GENOME_OUTPUT"

DB_ROOT="$(margie_resolve_runtime_db_root "$ROOT_DIR")"
RAST_DB="${DB_ROOT}/rasttk"

echo "=========================================="
echo "BV-BRC RASTtk Gene Calling"
echo "=========================================="
echo "Genome:          $GENOME_NAME"
echo "Input:           $INPUT_FASTA"
echo "Output:          $GENOME_OUTPUT"
echo "Scientific name: $SCIENTIFIC_NAME"
echo "Genetic code:    $GENETIC_CODE"
echo "Domain:          $DOMAIN"
echo "Threads:         $THREADS"
echo "Runtime DB root: $DB_ROOT"
echo ""

if ! command -v docker >/dev/null 2>&1 || ! docker info >/dev/null 2>&1; then
    echo "Error: Docker daemon is required for run_annotation.sh"
    exit 1
fi

if ! docker images --format "{{.Repository}}:{{.Tag}}" | grep -q "^${IMAGE_NAME}$"; then
    echo "Error: Container image not found: ${IMAGE_NAME}"
    echo "Build with: ./setup_containers.sh --docker --skip-db"
    exit 1
fi

TMP_INPUT="$(mktemp -d)"
cp "$INPUT_FASTA" "$TMP_INPUT/genome.fasta"

docker run --rm \
    --platform linux/amd64 \
    -v "$TMP_INPUT:/container/input:ro" \
    -v "$GENOME_OUTPUT:/container/output:rw" \
    -v "$RAST_DB:/container/db:ro" \
    -w /container \
    "$IMAGE_NAME" \
    /container/scripts/run_rasttk_incremental.sh \
        "$GENOME_NAME" \
        "/container/input/genome.fasta" \
        "/container/output" \
        "$SCIENTIFIC_NAME" \
        "$GENETIC_CODE" \
        "$DOMAIN"

rm -rf "$TMP_INPUT"

docker run --rm \
    --platform linux/amd64 \
    -v "$GENOME_OUTPUT:/container/output:rw" \
    -v "$RAST_DB:/container/db:ro" \
    -w /container \
    "$IMAGE_NAME" \
    python3 /container/scripts/generate_rast_tsv.py \
        "/container/output" \
        "$GENOME_NAME"

echo ""
echo "✓ RASTtk annotation complete"
echo "Output files:"
echo "  Gene calls:       $GENOME_OUTPUT/gene_calls/"
echo "  Native outputs:   $GENOME_OUTPUT/native/"
echo "  Annotated TSV:    $GENOME_OUTPUT/rast.tsv"
