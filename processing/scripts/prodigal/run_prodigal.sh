#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $0 --input <genome.fna> --output <dir> --name <genome_name> [--mode single|meta]
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
SIF_PATH="$ROOT_DIR/processing/containers/prodigal/prodigal.sif"

INPUT_FILE=""
OUTPUT_DIR=""
GENOME_NAME=""
MODE="single"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input|-i) INPUT_FILE="$2"; shift 2 ;;
        --output|-o) OUTPUT_DIR="$2"; shift 2 ;;
        --name|-n) GENOME_NAME="$2"; shift 2 ;;
        --mode) MODE="$2"; shift 2 ;;
        *) usage ;;
    esac
done

if [[ -z "$INPUT_FILE" || -z "$OUTPUT_DIR" || -z "$GENOME_NAME" ]]; then
    usage
fi

bash "$SCRIPT_DIR/preflight_check_prodigal.sh" --input "$INPUT_FILE" --output "$OUTPUT_DIR"

INPUT_FILE="$(cd "$(dirname "$INPUT_FILE")" && pwd)/$(basename "$INPUT_FILE")"
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"

OUT_PREFIX="$OUTPUT_DIR/$GENOME_NAME"

if command -v prodigal >/dev/null 2>&1; then
    prodigal \
        -i "$INPUT_FILE" \
        -a "${OUT_PREFIX}.faa" \
        -d "${OUT_PREFIX}.fna" \
        -f gff \
        -o "${OUT_PREFIX}.gff" \
        -p "$MODE"
elif command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
    docker run --rm \
        --platform linux/amd64 \
        -v "$(dirname "$INPUT_FILE"):/container/input:ro" \
        -v "$OUTPUT_DIR:/container/output:rw" \
        prodigal-annotation:1.0 \
        prodigal \
            -i "/container/input/$(basename "$INPUT_FILE")" \
            -a "/container/output/${GENOME_NAME}.faa" \
            -d "/container/output/${GENOME_NAME}.fna" \
            -f gff \
            -o "/container/output/${GENOME_NAME}.gff" \
            -p "$MODE"
else
    RUNTIME_BIN=""
    if command -v apptainer >/dev/null 2>&1; then
        RUNTIME_BIN="apptainer"
    elif command -v singularity >/dev/null 2>&1; then
        RUNTIME_BIN="singularity"
    fi

    if [[ -z "$RUNTIME_BIN" || ! -f "$SIF_PATH" ]]; then
        echo "Error: no supported runtime path available for Prodigal"
        exit 1
    fi

    "$RUNTIME_BIN" exec \
        --bind "$(dirname "$INPUT_FILE"):/container/input" \
        --bind "$OUTPUT_DIR:/container/output" \
        "$SIF_PATH" \
        prodigal \
            -i "/container/input/$(basename "$INPUT_FILE")" \
            -a "/container/output/${GENOME_NAME}.faa" \
            -d "/container/output/${GENOME_NAME}.fna" \
            -f gff \
            -o "/container/output/${GENOME_NAME}.gff" \
            -p "$MODE"
fi

echo "Prodigal output: ${OUT_PREFIX}.faa, ${OUT_PREFIX}.fna, ${OUT_PREFIX}.gff"
