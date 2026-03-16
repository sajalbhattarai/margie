#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $0 --input <genome.fna> [--output <dir>]
EOF
    exit 1
}

INPUT_FILE=""
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input|-i) INPUT_FILE="$2"; shift 2 ;;
        --output|-o) OUTPUT_DIR="$2"; shift 2 ;;
        *) usage ;;
    esac
done

if [[ -z "$INPUT_FILE" ]]; then
    usage
fi

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "[FAIL] Input genome not found: $INPUT_FILE"
    exit 1
fi

if [[ -n "$OUTPUT_DIR" ]]; then
    mkdir -p "$OUTPUT_DIR"
    if [[ ! -w "$OUTPUT_DIR" ]]; then
        echo "[FAIL] Output directory not writable: $OUTPUT_DIR"
        exit 1
    fi
fi

SOURCE="${BASH_SOURCE[0]}"
while [[ -L "$SOURCE" ]]; do
    SRC_DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"
    SOURCE="$(readlink "$SOURCE")"
    [[ "$SOURCE" != /* ]] && SOURCE="$SRC_DIR/$SOURCE"
done
SCRIPT_DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
SIF_PATH="$ROOT_DIR/processing/containers/prodigal/prodigal.sif"

if command -v prodigal >/dev/null 2>&1; then
    echo "[PASS] Host prodigal available"
    exit 0
fi

if command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
    if docker image inspect prodigal-annotation:1.0 >/dev/null 2>&1; then
        echo "[PASS] Docker runtime + prodigal-annotation:1.0 image available"
        exit 0
    fi
    echo "[WARN] Docker available but prodigal image missing: prodigal-annotation:1.0"
fi

if command -v apptainer >/dev/null 2>&1 || command -v singularity >/dev/null 2>&1; then
    if [[ -f "$SIF_PATH" ]]; then
        echo "[PASS] Apptainer/Singularity runtime + SIF available"
        exit 0
    fi
    echo "[WARN] Apptainer/Singularity found but SIF missing: $SIF_PATH"
fi

echo "[FAIL] No executable path found for Prodigal"
echo "  Options:"
echo "  1) Install host prodigal"
echo "  2) Build docker image: ./setup_containers.sh --docker --skip-db"
echo "  3) Build/pull apptainer image: ./setup_containers.sh --apptainer --skip-db"
exit 1
