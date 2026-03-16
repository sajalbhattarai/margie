#!/bin/bash
# Run BV-BRC RASTtk Container
# Mounts input, output, and database directories for annotation

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
source "$ROOT_DIR/processing/scripts/lib/runtime_db_resolver.sh"

IMAGE_NAME="rasttk-annotation:1.0"

# Default directories (relative to workspace root)
INPUT_DIR="${ROOT_DIR}/input"
OUTPUT_DIR="${ROOT_DIR}/output"
DB_DIR="$(margie_resolve_runtime_db_root "$ROOT_DIR")/rasttk"

# Parse arguments
COMMAND=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --input|-i)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output|-o)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --db|-d)
            DB_DIR="$2"
            shift 2
            ;;
        *)
            # Everything else is the command to run in container
            COMMAND="$@"
            break
            ;;
    esac
done

# Validate directories exist
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory not found: $INPUT_DIR"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

if [ ! -d "$DB_DIR" ]; then
    echo "Warning: Database directory not found: $DB_DIR"
    echo "RASTtk subsystem enrichment will not be available"
fi

echo "Runtime DB source: $(dirname "$DB_DIR")"

# Run container with mounted volumes
docker run --rm \
    --platform linux/amd64 \
    -v "$INPUT_DIR:/container/input:ro" \
    -v "$OUTPUT_DIR:/container/output:rw" \
    -v "$DB_DIR:/container/db:ro" \
    -w /container \
    ${IMAGE_NAME} \
    ${COMMAND}
