#!/bin/bash

# Batch Annotation Pipeline with Proper Formatting
# Uses complete annotation pipelines that produce formatted, mergeable outputs
# Date: February 11, 2026

set -e

echo "=============================================================================="
echo "  Batch Genome Annotation Pipeline (Formatted Outputs)"
echo "=============================================================================="
echo "Start time: $(date)"
echo ""

# Configuration  
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${ROOT_DIR}/processing/scripts/lib/runtime_db_resolver.sh"

# Directories
RASTTK_DIR="${ROOT_DIR}/output/rasttk"
MARGIE_DIR="${ROOT_DIR}"
MARGIE_INPUT="${MARGIE_DIR}/input"
MARGIE_OUTPUT="${MARGIE_DIR}/output"
CONTAINER_NAME="margie-annotation-run"

# Container paths
CONTAINER_INPUT="/margie/input"
CONTAINER_OUTPUT="/margie/output"
CONTAINER_SCRIPTS="/margie/scripts"
HOST_DB_ROOT="$(margie_resolve_runtime_db_root "$ROOT_DIR")"

# Databases to run
DATABASES=("cog" "pfam" "tigrfam" "kegg" "cazy" "dbcan" "merops" "tcdb" "eggnog" "operon")

# Threads
THREADS=8

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Prepare genome input
prepare_genome_input() {
    local genome="$1"
    local input_faa="${RASTTK_DIR}/${genome}/gene_calls/${genome}.faa"
    local genome_input_dir="${MARGIE_INPUT}/${genome}"
    
    mkdir -p "${genome_input_dir}"
    
    if [ ! -f "${genome_input_dir}/${genome}.faa" ]; then
        if [ -f "$input_faa" ]; then
            cp "${input_faa}" "${genome_input_dir}/${genome}.faa"
            echo "  Copied ${genome}.faa"
        else
            echo -e "${YELLOW}⚠ ${input_faa} not found${NC}"
            return 1
        fi
    fi
    return 0
}

# Run annotation with proper formatting
run_annotation() {
    local genome="$1"
    local db="$2"
    local container_input="${CONTAINER_INPUT}/${genome}/${genome}.faa"
    local host_output="${MARGIE_OUTPUT}/${db}/${genome}"
    local container_output="${CONTAINER_OUTPUT}/${db}/${genome}"
    
    echo ""
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${GREEN}  ${genome} → ${db}${NC}"
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    
    mkdir -p "${host_output}"
    
    # Check if formatted output exists
    local parsed_file="${host_output}/${genome}_${db}_parsed.tsv"
    if [ -f "${parsed_file}" ]; then
        echo -e "${YELLOW}✓ Parsed output exists, skipping${NC}"
        return 0
    fi
    
    # Run complete pipeline based on database
    case "$db" in
        cog)
            docker exec "${CONTAINER_NAME}" bash -c "
                ${CONTAINER_SCRIPTS}/cog/one_stop_cog_pipeline.sh \
                    -i ${container_input} \
                    -o ${container_output} \
                    -t ${THREADS}
            "
            ;;
            
        pfam)
            docker exec "${CONTAINER_NAME}" bash -c "
                ${CONTAINER_SCRIPTS}/pfam/one_stop_pfam_pipeline.sh \
                    -i ${container_input} \
                    -o ${container_output} \
                    -t ${THREADS}
            "
            ;;
            
        tigrfam)
            docker exec "${CONTAINER_NAME}" bash -c "
                ${CONTAINER_SCRIPTS}/tigrfam/one_stop_tigr_pipeline.sh \
                    -i ${container_input} \
                    -o ${container_output} \
                    -t ${THREADS}
            "
            ;;
            
        kegg)
            docker exec "${CONTAINER_NAME}" bash -c "
                ${CONTAINER_SCRIPTS}/kegg/one_stop_kegg_pipeline.sh \
                    -i ${container_input} \
                    -o ${container_output} \
                    -t ${THREADS}
            "
            ;;
            
        cazy)
            docker exec "${CONTAINER_NAME}" bash -c "
                ${CONTAINER_SCRIPTS}/cazy/one_stop_cazy_pipeline.sh \
                    -i ${container_input} \
                    -o ${container_output} \
                    -t ${THREADS}
            "
            ;;
            
        dbcan)
            docker exec "${CONTAINER_NAME}" bash -c "
                ${CONTAINER_SCRIPTS}/dbcan/one_stop_dbcan_pipeline.sh \
                    -i ${container_input} \
                    -o ${container_output} \
                    -t ${THREADS}
            "
            ;;
            
        eggnog)
            docker exec "${CONTAINER_NAME}" bash -c "
                ${CONTAINER_SCRIPTS}/eggnog/one_stop_eggnog_pipeline.sh \
                    -i ${container_input} \
                    -o ${container_output} \
                    -t ${THREADS}
            "
            ;;

        operon)
            local gff_input="${RASTTK_DIR}/${genome}/gene_calls/${genome}.gff"
            if [ ! -f "${gff_input}" ]; then
                echo -e "${YELLOW}⚠ GFF not found for operon calling: ${gff_input}${NC}"
                return 1
            fi

            local operon_script="${ROOT_DIR}/processing/containers/operon/scripts/predict_operons.sh"
            if [ ! -x "${operon_script}" ]; then
                chmod +x "${operon_script}" 2>/dev/null || true
            fi

            bash "${operon_script}" \
                --gff "${gff_input}" \
                --output "${MARGIE_OUTPUT}/operon/${genome}"
            ;;
            
        merops|tcdb)
            # These may need simple formatter scripts - check if they exist
            echo -e "${YELLOW}⚠ ${db} formatter not yet implemented, skipping${NC}"
            return 1
            ;;
    esac
    
    local exit_code=$?
    if [ $exit_code -eq 0 ]; then
        echo -e "${GREEN}✓ Complete (formatted output generated)${NC}"
        return 0
    else
        echo -e "${YELLOW}⚠ Failed (exit ${exit_code})${NC}"
        return 1
    fi
}

# Main
GENOMES=($(ls -1 "${RASTTK_DIR}"))

echo "Genomes: ${GENOMES[*]}"
echo "Databases: ${DATABASES[*]}"
echo "Threads: ${THREADS}"
echo "Host runtime DB source: ${HOST_DB_ROOT}"
echo ""

# Check container
if ! docker ps | grep -q "${CONTAINER_NAME}"; then
    echo -e "${YELLOW}⚠ Container not running${NC}"
    exit 1
fi

echo "DB source policy: prefer /depot/lindems/data/margie/db, fallback to ${MARGIE_DIR}/db"
echo "Resolved host DB source: ${HOST_DB_ROOT}"

# Prepare inputs
echo "Preparing inputs..."
for genome in "${GENOMES[@]}"; do
    prepare_genome_input "$genome" && echo -e "${GREEN}✓ ${genome}${NC}"
done
echo ""

# Run annotations
total=$((${#GENOMES[@]} * ${#DATABASES[@]}))
completed=0
failed=0

for genome in "${GENOMES[@]}"; do
    for db in "${DATABASES[@]}"; do
        if run_annotation "$genome" "$db"; then
            ((completed++))
        else
            ((failed++))
        fi
    done
done

echo ""
echo "=============================================================================="
echo "Complete: ${completed}/${total}, Failed: ${failed}"
echo "Formatted results: ${MARGIE_OUTPUT}/"
echo "End: $(date)"
echo "=============================================================================="
