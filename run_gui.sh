#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v streamlit >/dev/null 2>&1; then
    echo "Streamlit is not installed."
    echo "Install with: python3 -m pip install -r $SCRIPT_DIR/requirements-gui.txt"
    exit 1
fi

exec streamlit run "$SCRIPT_DIR/app/gui/streamlit_app.py"
