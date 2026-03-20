#!/usr/bin/env bash
# generate_config.sh
# Dynamically generates config.yml for KofamScan using mounted database

set -euo pipefail

DB_DIR="${1:-/container/db}"
OUTPUT_FILE="${2:-/tmp/kegg_config.yml}"
CPU="${3:-8}"

resolve_kofam_db_dir() {
    local base_dir="$1"
    if [ -d "${base_dir}/profiles" ] && [ -f "${base_dir}/ko_list" ]; then
        echo "$base_dir"
        return 0
    fi

    if [ -d "${base_dir}/kegg/profiles" ] && [ -f "${base_dir}/kegg/ko_list" ]; then
        echo "${base_dir}/kegg"
        return 0
    fi

    return 1
}

# Validate database directory
if [ ! -d "$DB_DIR" ]; then
    echo "ERROR: Database directory not found: $DB_DIR"
    exit 1
fi

if ! DB_DIR_RESOLVED="$(resolve_kofam_db_dir "$DB_DIR")"; then
    echo "ERROR: Could not find KOfam database under ${DB_DIR}"
    echo "Expected either:"
    echo "  - ${DB_DIR}/profiles and ${DB_DIR}/ko_list"
    echo "  - ${DB_DIR}/kegg/profiles and ${DB_DIR}/kegg/ko_list"
    echo "Database may not be properly set up"
    exit 1
fi

# Generate config.yml from official template if present.
TEMPLATE_PATH="/opt/kofam_scan/config-template.yml"
if [ -f "$TEMPLATE_PATH" ]; then
    cp "$TEMPLATE_PATH" "$OUTPUT_FILE"
else
    cat > "$OUTPUT_FILE" <<EOF
profile: ${DB_DIR_RESOLVED}/profiles
ko_list: ${DB_DIR_RESOLVED}/ko_list
hmmsearch: /usr/local/bin/hmmsearch
parallel: /usr/bin/parallel
cpu: ${CPU}
EOF
fi

# Replace keys while preserving official template structure where possible.
# The official kofamscan config-template.yml ships with all entries commented out
# (e.g. "# ko_list: /path/to/ko_list"), so we must match both commented and
# uncommented forms.  Use [# ]* to cover both cases.
sed -i.bak -E "s|^[# ]*profile:.*$|profile: ${DB_DIR_RESOLVED}/profiles|" "$OUTPUT_FILE" || true
sed -i.bak -E "s|^[# ]*ko_list:.*$|ko_list: ${DB_DIR_RESOLVED}/ko_list|" "$OUTPUT_FILE" || true
sed -i.bak -E "s|^[# ]*hmmsearch:.*$|hmmsearch: /usr/local/bin/hmmsearch|" "$OUTPUT_FILE" || true
sed -i.bak -E "s|^[# ]*parallel:.*$|parallel: /usr/bin/parallel|" "$OUTPUT_FILE" || true
if grep -qE '^[# ]*cpu:' "$OUTPUT_FILE"; then
    sed -i.bak -E "s|^[# ]*cpu:.*$|cpu: ${CPU}|" "$OUTPUT_FILE" || true
else
    printf '\ncpu: %s\n' "$CPU" >> "$OUTPUT_FILE"
fi
rm -f "$OUTPUT_FILE.bak"
# Safety check: ensure required keys are actually present (un-commented)
if ! grep -qE '^ko_list:' "$OUTPUT_FILE"; then
    printf '\nko_list: %s/ko_list\n' "${DB_DIR_RESOLVED}" >> "$OUTPUT_FILE"
fi
if ! grep -qE '^profile:' "$OUTPUT_FILE"; then
    printf '\nprofile: %s/profiles\n' "${DB_DIR_RESOLVED}" >> "$OUTPUT_FILE"
fi

echo "Generated KofamScan config at ${OUTPUT_FILE}"
echo "Resolved KOfam DB directory: ${DB_DIR_RESOLVED}"
