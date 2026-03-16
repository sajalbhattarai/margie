#!/usr/bin/env bash
# GUI container entrypoint.
# Verifies the Streamlit app is reachable and starts it.
set -euo pipefail

APP_PATH="${MARGIE_APP_PATH:-/margie/app/gui/streamlit_app.py}"
PORT="${STREAMLIT_SERVER_PORT:-8501}"

if [[ ! -f "$APP_PATH" ]]; then
    echo "ERROR: Streamlit app not found at $APP_PATH"
    echo ""
    echo "The Margie repository must be mounted into the container."
    echo "Example:"
    echo "  docker run --rm -it -p ${PORT}:${PORT} \\"
    echo "      -v /path/to/margie:/margie:ro \\"
    echo "      -v /path/to/margie/input:/margie/input:rw \\"
    echo "      -v /path/to/margie/output:/margie/output:rw \\"
    echo "      margie-gui:1.0"
    exit 1
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Margie GUI"
echo "  App : $APP_PATH"
echo "  Port: $PORT"
echo "  Open http://localhost:${PORT} in your browser"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

exec streamlit run "$APP_PATH" \
    --server.port "$PORT" \
    --server.headless true \
    --browser.gatherUsageStats false \
    "$@"
