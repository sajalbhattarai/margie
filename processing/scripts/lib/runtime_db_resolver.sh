#!/usr/bin/env bash
# Resolve runtime database source for execution-time pipeline steps.
# Order:
# 1) Shared DB root (default: /depot/lindems/data/margie/db) when present and non-empty
# 2) Local repo DB root (<root-folder>/db)

set -euo pipefail

margie_dir_has_content() {
    local d="$1"
    [[ -d "$d" ]] && find "$d" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null | grep -q .
}

margie_runtime_shared_db_root() {
    echo "${MARGIE_SHARED_DB_ROOT:-/depot/lindems/data/margie/db}"
}

margie_resolve_runtime_db_root() {
    local root_dir="$1"
    local shared_db
    shared_db="$(margie_runtime_shared_db_root)"

    if margie_dir_has_content "$shared_db"; then
        echo "$shared_db"
    else
        echo "$root_dir/db"
    fi
}
