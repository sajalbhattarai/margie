#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TARGET_DIR="${MARGIE_BIN_DIR:-$HOME/bin}"
QUIET=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --quiet)
            QUIET=1
            shift
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

mkdir_cmd() {
    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "[DRY-RUN] mkdir -p $TARGET_DIR"
    else
        mkdir -p "$TARGET_DIR"
    fi
}

link_cmd() {
    local source_path="$1"
    local link_name="$2"
    local link_path="$TARGET_DIR/$link_name"
    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "[DRY-RUN] ln -sfn $source_path $link_path"
    else
        ln -sfn "$source_path" "$link_path"
    fi
}

mkdir_cmd

link_cmd "$SCRIPT_DIR/setup_containers" "setup_containers"
link_cmd "$SCRIPT_DIR/setup_databases" "setup_databases"
link_cmd "$SCRIPT_DIR/annotate-genome" "annotate-genome"

if [[ "$QUIET" -eq 0 ]]; then
    echo "Installed Margie commands to: $TARGET_DIR"
    echo "Available commands: setup_containers, setup_databases, annotate-genome"
fi
