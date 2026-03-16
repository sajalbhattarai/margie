#!/bin/bash

# Batch Annotation Pipeline for BV-BRC Genomes
# Date: February 11, 2026

set -e

echo "=============================================================================="
echo "  Batch Genome Annotation Pipeline"
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
CONTAINER_DB=""
HOST_DB_ROOT="$(margie_resolve_runtime_db_root "$ROOT_DIR")"

# Databases and threads
DATABASES=("cog" "pfam" "tigrfam" "kegg" "dbcan" "cazy" "merops" "tcdb")
THREADS=8

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

resolve_container_db_root() {
    docker exec "${CONTAINER_NAME}" bash -lc '
        if [ -d /depot/lindems/data/margie/db ] && [ "$(ls -A /depot/lindems/data/margie/db 2>/dev/null)" ]; then
            echo /depot/lindems/data/margie/db
        elif [ -d /margie/db ] && [ "$(ls -A /margie/db 2>/dev/null)" ]; then
            echo /margie/db
        elif [ -d /margie/databases ] && [ "$(ls -A /margie/databases 2>/dev/null)" ]; then
            echo /margie/databases
        else
            exit 1
        fi
    '
}

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

# Run annotation
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
    
    # Check if exists
    local outfile="${host_output}/${genome}_${db}_raw.tsv"
    [[ "$db" =~ ^(pfam|tigrfam|dbcan)$ ]] && outfile="${host_output}/${genome}_${db}.tbl"
    
    if [ -f "${outfile}" ]; then
        echo -e "${YELLOW}✓ Output exists, skipping${NC}"
        return 0
    fi
    
    # Run based on DB type
    case "$db" in
        cog|kegg|cazy|merops)
            local db_path="${CONTAINER_DB}/${db}"
            [[ "$db" = "cog" ]] && db_path="${db_path}/COG.dmnd"
            [[ "$db" = "kegg" ]] && db_path="${db_path}/kegg.dmnd"
            [[ "$db" = "cazy" ]] && db_path="${db_path}/CAZy.dmnd"
            [[ "$db" = "merops" ]] && db_path="${db_path}/merops.dmnd"
            
            local max_targets=1
            local evalue="1e-5"
            [[ "$db" = "cazy" ]] && max_targets=5 && evalue="1e-10"
            
            docker exec "${CONTAINER_NAME}" bash -c "mkdir -p ${container_output} && cd ${container_output} && diamond blastp --query ${container_input} --db ${db_path} --out ${genome}_${db}_raw.tsv --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle --threads ${THREADS} --max-target-seqs ${max_targets} --evalue ${evalue} && echo '${db} complete'"
            ;;

        tcdb)
            local db_path="${CONTAINER_DB}/tcdb/tcdb_blast"
            docker exec "${CONTAINER_NAME}" bash -c "mkdir -p ${container_output} && cd ${container_output} && blastp -query ${container_input} -db ${db_path} -out ${genome}_${db}_raw.tsv -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -num_threads ${THREADS} -max_target_seqs 1 -evalue 1e-5 && echo 'tcdb complete'"
            ;;
            
        pfam)
            docker exec "${CONTAINER_NAME}" bash -c "mkdir -p ${container_output} && cd ${container_output} && hmmscan --tblout ${genome}_pfam.tbl --domtblout ${genome}_pfam_dom.tbl --cut_ga --cpu ${THREADS} ${CONTAINER_DB}/pfam/Pfam-A.hmm ${container_input} > ${genome}_pfam_full.txt && echo 'Pfam complete'"
            ;;
            
        tigrfam)
            docker exec "${CONTAINER_NAME}" bash -c "mkdir -p ${container_output} && cd ${container_output} && hmmscan --tblout ${genome}_tigrfam.tbl --domtblout ${genome}_tigrfam_dom.tbl --cut_tc --cpu ${THREADS} ${CONTAINER_DB}/tigrfam/TIGRFAMs_15.0_HMM.LIB ${container_input} > ${genome}_tigrfam_full.txt && echo 'TIGRFam complete'"
            ;;
            
        dbcan)
            docker exec "${CONTAINER_NAME}" bash -c "mkdir -p ${container_output} && cd ${container_output} && hmmscan --tblout ${genome}_dbcan.tbl --domtblout ${genome}_dbcan_dom.tbl --cpu ${THREADS} ${CONTAINER_DB}/dbcan/dbCAN_sub.hmm ${container_input} > ${genome}_dbcan_full.txt && echo 'dbCAN complete'"
            ;;
    esac
    
    local exit_code=$?
    if [ $exit_code -eq 0 ]; then
        echo -e "${GREEN}✓ Complete${NC}"
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

if ! CONTAINER_DB="$(resolve_container_db_root)"; then
    echo -e "${YELLOW}⚠ No database mount found in container. Expected /depot/lindems/data/margie/db or /margie/db${NC}"
    exit 1
fi
echo "Container DB source: ${CONTAINER_DB}"

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
echo "Results: ${MARGIE_OUTPUT}/"
echo "End: $(date)"
echo "=============================================================================="
