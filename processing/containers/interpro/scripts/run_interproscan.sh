#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_interproscan.sh --input proteins.faa --output OUT_DIR --db DB_DIR [--threads N]

Required:
  --input      Protein FASTA input (.faa/.fasta/.fa)
  --output     Output directory for InterProScan results
  --db         InterPro database root (contains interproscan-*/interproscan.sh)

Optional:
  --threads    CPU count (default: 1)

Environment:
  MARGIE_INTERPRO_DISABLE_LOOKUP=1   Pass -dp to disable remote lookup service
EOF
}

INPUT=""
OUTPUT_DIR=""
DB_DIR=""
THREADS="1"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input) INPUT="${2:-}"; shift 2 ;;
        --output) OUTPUT_DIR="${2:-}"; shift 2 ;;
        --db) DB_DIR="${2:-}"; shift 2 ;;
        --threads) THREADS="${2:-1}"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

if [[ -z "$INPUT" || -z "$OUTPUT_DIR" || -z "$DB_DIR" ]]; then
    usage
    exit 1
fi

if [[ ! -f "$INPUT" ]]; then
    echo "Input FASTA not found: $INPUT" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

INTERPRO_DIR=""
if compgen -G "$DB_DIR/interproscan-*/interproscan.sh" >/dev/null; then
    INTERPRO_DIR="$(dirname "$(ls -1 "$DB_DIR"/interproscan-*/interproscan.sh | head -n1)")"
fi

if [[ -z "$INTERPRO_DIR" || ! -x "$INTERPRO_DIR/interproscan.sh" ]]; then
    echo "InterProScan launcher not found or not executable under: $DB_DIR" >&2
    exit 1
fi

OUT_TSV="$OUTPUT_DIR/interproscan.tsv"
TEMP_DIR="$OUTPUT_DIR/interproscan_tmp"
mkdir -p "$TEMP_DIR"

CMD=(
    "$INTERPRO_DIR/interproscan.sh"
    -i "$INPUT"
    -f tsv
    -cpu "$THREADS"
    -T "$TEMP_DIR"
    -o "$OUT_TSV"
)

if [[ "${MARGIE_INTERPRO_DISABLE_LOOKUP:-0}" == "1" ]]; then
    CMD+=( -dp )
fi

"${CMD[@]}"

if [[ ! -s "$OUT_TSV" ]]; then
    echo "InterProScan did not produce output: $OUT_TSV" >&2
    exit 1
fi

echo "InterProScan output written to: $OUT_TSV"
