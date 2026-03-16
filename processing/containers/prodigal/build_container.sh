#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

docker build --platform linux/amd64 -t prodigal-annotation:1.0 .

echo "Built prodigal-annotation:1.0"
