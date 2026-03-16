#!/bin/bash
# Backward-compatible wrapper for legacy callers.
# The active COG workflow is now run_cogclassifier.sh.

set -euo pipefail

echo "[WARN] run_diamond.sh is deprecated for COG. Redirecting to run_cogclassifier.sh" >&2
exec /container/scripts/run_cogclassifier.sh "$@"
