#!/usr/bin/env bash
# setup_databases.sh
# Downloads all required annotation databases.
# Run once after cloning the repository.
# Usage:
#   ./setup_databases.sh
#   ./setup_databases.sh --db cog,pfam,tigrfam
#   ./setup_databases.sh --cog --pfam --tigrfam
#   ./setup_databases.sh --all | -all
#   ./setup_databases.sh --runtime docker
#   ./setup_databases.sh --platform auto|linux/amd64|linux/arm64
#   ./setup_databases.sh --dry-run
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DB_DIR="$SCRIPT_DIR/db"
STATE_DIR="$SCRIPT_DIR/.margie_state"
STATE_DATABASES_FILE="$STATE_DIR/databases.ready"

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; BOLD='\033[1m'; NC='\033[0m'

DRY_RUN=0
SELECTED="all"
MODE="auto"
DOCKER_PLATFORM="auto"
FORCE_ALL_DBS=0
REQUESTED_DBS=()

append_requested_db() {
    local db_name="$1"
    local existing
    for existing in "${REQUESTED_DBS[@]-}"; do
        if [[ "$existing" == "$db_name" ]]; then
            return 0
        fi
    done
    REQUESTED_DBS+=("$db_name")
}

append_requested_db_list() {
    local db_list="$1"
    local db_name
    IFS=',' read -ra db_names <<< "$db_list"
    for db_name in "${db_names[@]}"; do
        db_name="${db_name//[[:space:]]/}"
        if [[ -n "$db_name" ]]; then
            append_requested_db "$db_name"
        fi
    done
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --db)
            SELECTED="${2:-all}"
            if [[ "$SELECTED" == "all" ]]; then
                FORCE_ALL_DBS=1
            else
                append_requested_db_list "$SELECTED"
            fi
            shift 2
            ;;
        --all|-all|all)
            FORCE_ALL_DBS=1
            shift
            ;;
        --cog)
            append_requested_db "cog"
            shift
            ;;
        --pfam)
            append_requested_db "pfam"
            shift
            ;;
        --tigrfam)
            append_requested_db "tigrfam"
            shift
            ;;
        --dbcan)
            append_requested_db "dbcan"
            shift
            ;;
        --kegg)
            append_requested_db "kegg"
            shift
            ;;
        --eggnog)
            append_requested_db "eggnog"
            shift
            ;;
        --merops)
            append_requested_db "merops"
            shift
            ;;
        --tcdb)
            append_requested_db "tcdb"
            shift
            ;;
        --uniprot)
            append_requested_db "uniprot"
            shift
            ;;
        --interpro)
            append_requested_db "interpro"
            shift
            ;;
        --rasttk)
            append_requested_db "rasttk"
            shift
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --runtime)
            MODE="${2:-auto}"
            shift 2
            ;;
        --platform)
            DOCKER_PLATFORM="${2:-auto}"
            shift 2
            ;;
        all)
            SELECTED="all"
            shift
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            exit 1
            ;;
    esac
done

if [[ "$FORCE_ALL_DBS" -eq 1 ]]; then
    SELECTED="all"
elif [[ "${#REQUESTED_DBS[@]}" -gt 0 ]]; then
    IFS=',' SELECTED="${REQUESTED_DBS[*]}"
fi

# Mirror links provided in "links to databases.rtf"
COG_BASE="https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data"
CDD_BASE="https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd"
KOFAM_BASE="https://www.genome.jp/ftp/db/kofam"
KOFAM_FALLBACK_BASE="https://ftp.genome.jp/pub/db/kofam"
EGGNOG_BASE="http://eggnog5.embl.de/download/eggnog_5.0"
EGGNOG_EMAPPER_BASE="http://eggnog5.embl.de/download/emapperdb-5.0.2"
MEROPS_BASE="https://ftp.ebi.ac.uk/pub/databases/merops/current_release"
TCDB_URL="https://www.tcdb.org/download.php"
PFAM_BASE="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"
TIGRFAM_BASE="https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0"
PGAP_BASE="https://ftp.ncbi.nlm.nih.gov/hmm/current"
UNIPROT_SPROT_URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
INTERPRO_BASE="https://ftp.ebi.ac.uk/pub/databases/interpro"
INTERPROSCAN_VERSION="5.77-108.0"
INTERPROSCAN_ARCHIVE="interproscan-${INTERPROSCAN_VERSION}-64-bit.tar.gz"
INTERPROSCAN_MD5_FILE="${INTERPROSCAN_ARCHIVE}.md5"
INTERPROSCAN_SOFTWARE_BASE="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${INTERPROSCAN_VERSION}"
OFFICIAL_DBCAN_IMAGE="haidyi/run_dbcan:latest"
OFFICIAL_EGGNOG_IMAGE="quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0"

is_relaxed_tls_url() {
    local url="$1"
    [[ "$url" == *"bcb.unl.edu"* || "$url" == *"pro.unl.edu"* ]]
}

detect_runtime() {
    if [[ "$MODE" != "auto" ]]; then
        echo "$MODE"
        return 0
    fi
    if command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
        echo "docker"
    elif command -v apptainer >/dev/null 2>&1; then
        echo "apptainer"
    elif command -v singularity >/dev/null 2>&1; then
        echo "singularity"
    else
        echo "none"
    fi
}

RUNTIME="$(detect_runtime)"

resolve_docker_platform() {
    if [[ "$DOCKER_PLATFORM" != "auto" ]]; then
        echo "$DOCKER_PLATFORM"
        return 0
    fi

    local arch
    arch="$(uname -m 2>/dev/null || echo unknown)"
    case "$arch" in
        arm64|aarch64)
            echo "linux/arm64"
            ;;
        x86_64|amd64)
            echo "linux/amd64"
            ;;
        *)
            echo "linux/amd64"
            ;;
    esac
}

DOCKER_PLATFORM_RESOLVED="$(resolve_docker_platform)"

docker_platform_for_db() {
    case "$1" in
        eggnog) echo "linux/amd64" ;;
        dbcan) echo "linux/amd64" ;;
        *) echo "$DOCKER_PLATFORM_RESOLVED" ;;
    esac
}

echo -e "${BOLD}=== Margie Annotation Pipeline — Database Setup ===${NC}"
echo "Databases will be downloaded to: $DB_DIR"
echo "Runtime mode: $RUNTIME"
if [[ "$RUNTIME" == "docker" ]]; then
    echo "Docker platform: $DOCKER_PLATFORM_RESOLVED"
fi
if [[ "$DRY_RUN" -eq 1 ]]; then
    echo -e "${YELLOW}DRY-RUN: commands will be printed but not executed${NC}"
fi
echo ""

need_cmd() {
    if [[ "$DRY_RUN" -eq 1 ]]; then
        return
    fi
    if ! command -v "$1" >/dev/null 2>&1; then
        echo -e "${RED}Missing required command: $1${NC}"
        exit 1
    fi
}

have_cmd() {
    command -v "$1" >/dev/null 2>&1
}

python_download() {
    local url="$1"
    local out="${2:-}"
    local optional="${3:-0}"
    local relaxed_tls="${4:-0}"

    if [[ "$DRY_RUN" -eq 1 ]]; then
        if [[ -n "$out" ]]; then
            cmd python3 -c "download" "$url" "$out"
        else
            cmd python3 -c "download" "$url"
        fi
        return 0
    fi

    if ! have_cmd python3; then
        echo -e "${RED}Missing required command: python3${NC}"
        return 1
    fi

    python3 - "$url" "$out" "$relaxed_tls" <<'PY'
import os
import ssl
import sys
import urllib.request

url = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else ""
relaxed_tls = bool(int(sys.argv[3])) if len(sys.argv) > 3 else False

if out:
    target = out
else:
    target = url.rsplit("/", 1)[-1] or "downloaded_file"

tmp = target + ".part"
start = 0
if os.path.exists(tmp):
    start = os.path.getsize(tmp)

req = urllib.request.Request(url)
if start > 0:
    req.add_header("Range", f"bytes={start}-")

ctx = ssl._create_unverified_context() if relaxed_tls else None

with urllib.request.urlopen(req, context=ctx) as resp:
    length = resp.headers.get("Content-Length")
    total = int(length) + start if length is not None else None
    mode = "ab" if start > 0 else "wb"
    downloaded = start

    with open(tmp, mode) as f:
        while True:
            chunk = resp.read(1024 * 256)
            if not chunk:
                break
            f.write(chunk)
            downloaded += len(chunk)
            if total:
                pct = (downloaded / total) * 100
                sys.stdout.write(f"\rDownloading {target}: {pct:6.2f}% ({downloaded}/{total} bytes)")
            else:
                sys.stdout.write(f"\rDownloading {target}: {downloaded} bytes")
            sys.stdout.flush()

os.replace(tmp, target)
sys.stdout.write("\n")
PY
}

download_with_wget() {
    local url="$1"
    local out="${2:-}"
    local relaxed_tls="${3:-0}"

    local show_progress=0
    if wget --help 2>&1 | grep -q -- '--show-progress'; then
        show_progress=1
    fi

    local -a args
    args=(-c)
    if [[ "$show_progress" -eq 1 ]]; then
        args+=(--show-progress --progress=bar:force:noscroll)
    fi
    if [[ "$relaxed_tls" -eq 1 ]]; then
        args+=(--no-check-certificate)
    fi

    if [[ -n "$out" ]]; then
        wget "${args[@]}" "$url" -O "$out"
    else
        wget "${args[@]}" "$url"
    fi
}

warn_missing_indexer() {
    local tool="$1"
    local db_name="$2"
    echo -e "${YELLOW}  ! ${tool} not found. Download completed for ${db_name}, but local index creation was skipped.${NC}"
    echo -e "${YELLOW}    Install ${tool} or build containers and rerun: ./setup_databases.sh --db ${db_name}${NC}"
}

image_for_db() {
    case "$1" in
        cog) echo "cog-annotation:1.0" ;;
        eggnog) echo "$OFFICIAL_EGGNOG_IMAGE" ;;
        merops) echo "merops-annotation:1.0" ;;
        tcdb) echo "tcdb-annotation:1.0" ;;
        uniprot) echo "uniprot-annotation:1.0" ;;
        pfam) echo "pfam-annotation:1.0" ;;
        tigrfam) echo "tigrfam-annotation:1.0" ;;
        dbcan) echo "haidyi/run_dbcan:latest" ;;
        *) echo "" ;;
    esac
}

sif_for_db() {
    local db_name="$1"
    echo "$SCRIPT_DIR/processing/containers/${db_name}/${db_name}.sif"
}

run_indexer() {
    local db_name="$1"
    shift

    # Keep pre-indexed HMMER sidecars by default so copied DBs are reusable.
    # Set MARGIE_FORCE_HMMPRESS=1 to force regeneration.
    if [[ "$1" == "hmmpress" && "$#" -ge 2 ]]; then
        local hmm_target
        hmm_target="${@: -1}"

        if [[ "${MARGIE_FORCE_HMMPRESS:-0}" != "1" ]] \
            && [[ -s "$hmm_target" ]] \
            && [[ -s "${hmm_target}.h3f" ]] \
            && [[ -s "${hmm_target}.h3i" ]] \
            && [[ -s "${hmm_target}.h3m" ]] \
            && [[ -s "${hmm_target}.h3p" ]]; then
            echo -e "${YELLOW}  ! Existing HMMER index detected for ${hmm_target}; skipping hmmpress (set MARGIE_FORCE_HMMPRESS=1 to rebuild).${NC}"
            return 0
        fi

        # Make forced reruns idempotent by removing stale sidecar indices first.
        if [[ "$DRY_RUN" -eq 1 ]]; then
            cmd rm -f "${hmm_target}.h3f" "${hmm_target}.h3i" "${hmm_target}.h3m" "${hmm_target}.h3p"
        else
            rm -f "${hmm_target}.h3f" "${hmm_target}.h3i" "${hmm_target}.h3m" "${hmm_target}.h3p"
        fi
    fi

    if [[ "$DRY_RUN" -eq 1 ]]; then
        cmd "$@"
        return 0
    fi

    if have_cmd "$1"; then
        cmd "$@"
        return 0
    fi

    local image sif
    local docker_platform
    image="$(image_for_db "$db_name")"
    sif="$(sif_for_db "$db_name")"
    docker_platform="$(docker_platform_for_db "$db_name")"

    if [[ "$RUNTIME" == "docker" ]]; then
        if [[ -n "$image" ]] && docker image inspect "$image" >/dev/null 2>&1; then
            echo -e "${YELLOW}  ! Host $1 missing, using Docker image ${image} for indexing.${NC}"
            if ! docker run --rm --platform "$docker_platform" "$image" /bin/sh -lc "command -v '$1' >/dev/null 2>&1"; then
                warn_missing_indexer "$1" "$db_name"
                return 0
            fi
            if docker run --rm --platform "$docker_platform" -v "$DB_DIR:/db" -w "/db/${db_name}" "$image" "$@"; then
                return 0
            fi
            return 1
        fi
    elif [[ "$RUNTIME" == "apptainer" || "$RUNTIME" == "singularity" ]]; then
        if [[ -f "$sif" ]]; then
            echo -e "${YELLOW}  ! Host $1 missing, using ${RUNTIME} image ${sif} for indexing.${NC}"
            if ! "$RUNTIME" exec --bind "$DB_DIR:/db" --pwd "/db/${db_name}" "$sif" /bin/sh -lc "command -v '$1' >/dev/null 2>&1"; then
                warn_missing_indexer "$1" "$db_name"
                return 0
            fi
            if "$RUNTIME" exec --bind "$DB_DIR:/db" --pwd "/db/${db_name}" "$sif" "$@"; then
                return 0
            fi
            return 1
        fi
    fi

    warn_missing_indexer "$1" "$db_name"
    return 0
}

has_container_indexer_for_db() {
    local db_name="$1"
    local image sif
    image="$(image_for_db "$db_name")"
    sif="$(sif_for_db "$db_name")"

    if [[ "$RUNTIME" == "docker" ]]; then
        [[ -n "$image" ]] && docker image inspect "$image" >/dev/null 2>&1
        return $?
    fi

    if [[ "$RUNTIME" == "apptainer" || "$RUNTIME" == "singularity" ]]; then
        [[ -f "$sif" ]]
        return $?
    fi

    return 1
}

first_time_guidance() {
    if [[ "$DRY_RUN" -eq 1 ]]; then
        return 0
    fi

    local needs_diamond=(eggnog merops uniprot)
    local needs_hmm=(pfam tigrfam)
    local needs_help=0
    local db

    if ! have_cmd diamond; then
        for db in "${needs_diamond[@]}"; do
            if [[ "$SELECTED" == "all" || ",$SELECTED," == *",$db,"* ]]; then
                if ! has_container_indexer_for_db "$db"; then
                    needs_help=1
                    break
                fi
            fi
        done
    fi

    if [[ "$needs_help" -eq 0 ]] && ! have_cmd hmmpress; then
        for db in "${needs_hmm[@]}"; do
            if [[ "$SELECTED" == "all" || ",$SELECTED," == *",$db,"* ]]; then
                if ! has_container_indexer_for_db "$db"; then
                    needs_help=1
                    break
                fi
            fi
        done
    fi

    if [[ "$needs_help" -eq 1 ]]; then
        echo -e "${YELLOW}First-time setup note:${NC}"
        echo "  Some local indexers are missing and matching container images are not built yet."
        echo "  Recommended first-time command: ./setup_containers.sh"
        echo "  This builds containers first, then runs database setup automatically."
        echo ""
    fi
}

cmd() {
    if [[ "$DRY_RUN" -eq 1 ]]; then
        printf '[DRY-RUN]'
        printf ' %q' "$@"
        printf '\n'
        return 0
    fi
    "$@"
}

ensure_state_dir() {
    mkdir -p "$STATE_DIR"
}

write_database_state() {
    ensure_state_dir
    {
        printf 'timestamp=%s\n' "$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
        printf 'runtime=%s\n' "$RUNTIME"
        printf 'selected=%s\n' "$SELECTED"
        printf 'docker_platform=%s\n' "$DOCKER_PLATFORM_RESOLVED"
    } > "$STATE_DATABASES_FILE"
}

download() {
    local url="$1"
    local out="${2:-}"
    local relaxed_tls=0
    if is_relaxed_tls_url "$url"; then
        relaxed_tls=1
    fi

    if [[ "$DRY_RUN" -eq 1 ]]; then
        if [[ -n "$out" ]]; then
            cmd wget -c --show-progress --progress=bar:force:noscroll "$url" -O "$out"
        else
            cmd wget -c --show-progress --progress=bar:force:noscroll "$url"
        fi
        return 0
    fi

    if have_cmd wget; then
        if download_with_wget "$url" "$out" 0; then
            return 0
        fi

        if [[ "$relaxed_tls" -eq 1 ]]; then
            echo -e "${YELLOW}  ! TLS verification failed for $url; retrying with relaxed TLS.${NC}"
            if download_with_wget "$url" "$out" 1; then
                return 0
            fi
        fi

        echo -e "${YELLOW}  ! wget download failed; using built-in Python downloader with progress.${NC}"
    else
        echo -e "${YELLOW}  ! wget not found; using built-in Python downloader with progress.${NC}"
    fi

    if [[ -n "$out" ]]; then
        python_download "$url" "$out" 0 "$relaxed_tls"
    else
        python_download "$url" "" 0 "$relaxed_tls"
    fi
}

download_optional() {
    local url="$1"
    local out="${2:-}"
    local relaxed_tls=0
    if is_relaxed_tls_url "$url"; then
        relaxed_tls=1
    fi

    if [[ "$DRY_RUN" -eq 1 ]]; then
        if [[ -n "$out" ]]; then
            cmd wget -c --show-progress --progress=bar:force:noscroll "$url" -O "$out"
        else
            cmd wget -c --show-progress --progress=bar:force:noscroll "$url"
        fi
        return 0
    fi

    if have_cmd wget; then
        if download_with_wget "$url" "$out" 0; then
            return 0
        fi

        if [[ "$relaxed_tls" -eq 1 ]]; then
            if download_with_wget "$url" "$out" 1; then
                return 0
            fi
        fi

        echo -e "${YELLOW}  Optional wget download failed; trying Python downloader: $url${NC}"
    else
        echo -e "${YELLOW}  ! wget not found; using built-in Python downloader with progress for optional files.${NC}"
    fi

    if [[ -n "$out" ]]; then
        python_download "$url" "$out" 1 "$relaxed_tls" || echo -e "${YELLOW}  Optional file unavailable: $url${NC}"
    else
        python_download "$url" "" 1 "$relaxed_tls" || echo -e "${YELLOW}  Optional file unavailable: $url${NC}"
    fi
}

download_from_mirrors() {
    local out="$1"
    shift

    local url
    local tried=0
    for url in "$@"; do
        tried=1
        echo "  Trying source: $url"
        if download "$url" "$out"; then
            return 0
        fi
        echo -e "${YELLOW}  ! Source failed: $url${NC}"
    done

    if [[ "$tried" -eq 0 ]]; then
        echo -e "${RED}  No download mirrors were provided for $out${NC}"
    else
        echo -e "${RED}  Failed to download $out from all configured mirrors.${NC}"
    fi
    return 1
}

dbcan_build_official() {
    local cpus="${1:-4}"
    local image sif
    local docker_platform
    image="$(image_for_db "dbcan")"
    sif="$(sif_for_db "dbcan")"
    docker_platform="$(docker_platform_for_db "dbcan")"

    if [[ "$DRY_RUN" -eq 1 ]]; then
        if [[ "$RUNTIME" == "docker" ]]; then
            cmd docker run --rm --platform "$docker_platform" --entrypoint dbcan_build -v "$DB_DIR:/db" -w "/db/dbcan" "$image" --cpus "$cpus" --db-dir "/db/dbcan" --clean
            return 0
        fi
        if [[ "$RUNTIME" == "apptainer" || "$RUNTIME" == "singularity" ]]; then
            cmd "$RUNTIME" exec --bind "$DB_DIR:/db" --pwd "/db/dbcan" "$sif" dbcan_build --cpus "$cpus" --db-dir "/db/dbcan" --clean
            return 0
        fi
        return 1
    fi

    if [[ "$RUNTIME" == "docker" ]]; then
        if docker image inspect "$image" >/dev/null 2>&1; then
            echo -e "${YELLOW}  Using official dbCAN runner: $image${NC}"
            docker run --rm --platform "$docker_platform" --entrypoint dbcan_build -v "$DB_DIR:/db" -w "/db/dbcan" "$image" --cpus "$cpus" --db-dir "/db/dbcan" --clean
            return $?
        fi
        return 1
    fi

    if [[ "$RUNTIME" == "apptainer" || "$RUNTIME" == "singularity" ]]; then
        if [[ -f "$sif" ]]; then
            echo -e "${YELLOW}  Using official dbCAN runner SIF: $sif${NC}"
            "$RUNTIME" exec --bind "$DB_DIR:/db" --pwd "/db/dbcan" "$sif" dbcan_build --cpus "$cpus" --db-dir "/db/dbcan" --clean
            return $?
        fi
        return 1
    fi

    return 1
}

dbcan_run_tool() {
    local tool="$1"
    shift

    if [[ "$DRY_RUN" -eq 1 ]]; then
        if [[ "$RUNTIME" == "docker" ]]; then
            cmd docker run --rm --platform "$(docker_platform_for_db "dbcan")" --entrypoint "$tool" -v "$DB_DIR:/db" -w "/db/dbcan" "$OFFICIAL_DBCAN_IMAGE" "$@"
            return 0
        fi
        cmd "$tool" "$@"
        return 0
    fi

    if have_cmd "$tool"; then
        "$tool" "$@"
        return $?
    fi

    if [[ "$RUNTIME" == "docker" ]] && docker image inspect "$OFFICIAL_DBCAN_IMAGE" >/dev/null 2>&1; then
        echo -e "${YELLOW}  ! Host $tool missing, using Docker image $OFFICIAL_DBCAN_IMAGE.${NC}"
        docker run --rm --platform "$(docker_platform_for_db "dbcan")" --entrypoint "$tool" -v "$DB_DIR:/db" -w "/db/dbcan" "$OFFICIAL_DBCAN_IMAGE" "$@"
        return $?
    fi

    warn_missing_indexer "$tool" "dbcan"
    return 0
}

dbcan_build_fallback() {
    echo -e "${YELLOW}  ! Official dbCAN build failed. Falling back to direct mirrored downloads.${NC}"

    # dbcan_build --clean can remove /db/dbcan; recreate and ensure cwd exists.
    cmd mkdir -p "$DB_DIR/dbcan"
    [[ "$DRY_RUN" -eq 0 ]] && cd "$DB_DIR/dbcan"

    # Latest dbCAN resources from browse_download API on pro.unl.edu.
    download "https://pro.unl.edu/dbCAN2/download_file.php?file=Databases/V14/dbCAN-HMMdb-V14.txt" "dbCAN-HMMdb-V14.txt"
    download_optional "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/fam-substrate-mapping.tsv" "fam-substrate-mapping.tsv"
    download_optional "https://pro.unl.edu/dbCAN2/download_file.php?file=Databases/fam-substrate-mapping-08262025.tsv" "fam-substrate-mapping-08262025.tsv"
    download_optional "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/CAZyDB.07242025.fa" "CAZyDB.07242025.fa"
    download_optional "https://pro.unl.edu/dbCAN2/download_file.php?file=Databases/dbCAN_sub.hmm" "dbCAN_sub.hmm"
    download_optional "https://pro.unl.edu/dbCAN2/download_file.php?file=Databases/tcdb.fa" "tcdb.fa"

    if [[ "$DRY_RUN" -eq 1 ]]; then
        cmd cp dbCAN-HMMdb-V14.txt dbCAN.txt
        # Preserve historical filename expected by existing wrappers.
        cmd cp dbCAN-HMMdb-V14.txt dbCAN-HMMdb-V12.txt
    else
        cp dbCAN-HMMdb-V14.txt dbCAN.txt
        cp dbCAN-HMMdb-V14.txt dbCAN-HMMdb-V12.txt
    fi

    dbcan_run_tool hmmpress dbCAN-HMMdb-V14.txt

    if [[ -f dbCAN_sub.hmm ]]; then
        dbcan_run_tool hmmpress -f dbCAN_sub.hmm
    fi

    if [[ -f fam-substrate-mapping.tsv || -f fam-substrate-mapping-08262025.tsv ]]; then
        local src_map=""
        if [[ -f fam-substrate-mapping.tsv ]]; then
            src_map="fam-substrate-mapping.tsv"
        else
            src_map="fam-substrate-mapping-08262025.tsv"
        fi
        if [[ "$DRY_RUN" -eq 1 ]]; then
            cmd cp "$src_map" fam-substrate-mapping.tsv
            if [[ ! -f substrate-mappings.tsv ]]; then
                cmd cp fam-substrate-mapping.tsv substrate-mappings.tsv
            fi
        else
            cp "$src_map" fam-substrate-mapping.tsv
            if [[ ! -f substrate-mappings.tsv ]]; then
                cp fam-substrate-mapping.tsv substrate-mappings.tsv
            fi
        fi
    fi

    if [[ -f CAZyDB.07242025.fa ]]; then
        if [[ "$DRY_RUN" -eq 1 ]]; then
            cmd cp CAZyDB.07242025.fa CAZyDB.fa
        else
            cp CAZyDB.07242025.fa CAZyDB.fa
        fi
    fi

    if [[ -f tcdb.fa ]]; then
        dbcan_run_tool diamond makedb --in tcdb.fa -d tcdb || true
    fi
}

dbcan_pick_hmm_db() {
    if [[ -f dbCAN-HMMdb-V14.txt ]]; then
        echo "dbCAN-HMMdb-V14.txt"
        return 0
    fi
    if [[ -f dbCAN.txt ]]; then
        echo "dbCAN.txt"
        return 0
    fi
    if [[ -f dbCAN-HMMdb-V12.txt ]]; then
        echo "dbCAN-HMMdb-V12.txt"
        return 0
    fi
    return 1
}

dbcan_assets_ready() {
    local hmm_db=""
    if ! hmm_db="$(dbcan_pick_hmm_db)"; then
        return 1
    fi

    [[ -f "$hmm_db" ]] || return 1

    # Keep checks intentionally permissive: if these core assets exist,
    # skip re-downloading and reuse local dbCAN/CAZy files.
    if [[ ! -f CAZyDB.fa && ! -f CAZyDB.07242025.fa ]]; then
        return 1
    fi
    if [[ ! -f fam-substrate-mapping.tsv && ! -f fam-substrate-mapping-08262025.tsv && ! -f substrate-mappings.tsv ]]; then
        return 1
    fi

    return 0
}

dbcan_try_press_if_missing() {
    local hmm_db="$1"

    if [[ -n "$hmm_db" && ( ! -f "${hmm_db}.h3f" || ! -f "${hmm_db}.h3i" || ! -f "${hmm_db}.h3m" || ! -f "${hmm_db}.h3p" ) ]]; then
        echo -e "${YELLOW}  ! Existing dbCAN HMM database found but index is incomplete; attempting hmmpress.${NC}"
        dbcan_run_tool hmmpress "$hmm_db" || true
    fi

    if [[ -f dbCAN_sub.hmm && ( ! -f dbCAN_sub.hmm.h3f || ! -f dbCAN_sub.hmm.h3i || ! -f dbCAN_sub.hmm.h3m || ! -f dbCAN_sub.hmm.h3p ) ]]; then
        echo -e "${YELLOW}  ! Existing dbCAN_sub.hmm found but index is incomplete; attempting hmmpress.${NC}"
        dbcan_run_tool hmmpress -f dbCAN_sub.hmm || true
    fi
}

gunzip_if_present() {
    local gz="$1"
    if [[ "$DRY_RUN" -eq 1 ]]; then
        cmd gunzip -kf "$gz"
        return 0
    fi
    if [[ -f "$gz" ]]; then
        cmd gunzip -kf "$gz"
    fi
}

file_size_bytes() {
    local file="$1"
    wc -c < "$file" | tr -d ' '
}

require_nonempty_file() {
    local file="$1"
    local label="$2"
    if [[ ! -f "$file" ]]; then
        echo -e "${RED}  ✗ Missing required file: $label ($file).${NC}"
        return 1
    fi
    if [[ ! -s "$file" ]]; then
        echo -e "${RED}  ✗ File is empty/incomplete: $label ($file).${NC}"
        return 1
    fi
    return 0
}

require_min_size_bytes() {
    local file="$1"
    local min_bytes="$2"
    local label="$3"
    local size

    require_nonempty_file "$file" "$label" || return 1
    size="$(file_size_bytes "$file")"
    if [[ "$size" -lt "$min_bytes" ]]; then
        echo -e "${RED}  ✗ File too small (likely partial): $label ($file).${NC}"
        echo "    expected at least ${min_bytes} bytes, found ${size} bytes"
        return 1
    fi
    return 0
}

validate_gzip_file() {
    local file="$1"
    local label="$2"

    require_nonempty_file "$file" "$label" || return 1
    if [[ "$DRY_RUN" -eq 0 ]]; then
        if ! gzip -t "$file" >/dev/null 2>&1; then
            echo -e "${RED}  ✗ Corrupt/incomplete gzip archive: $label ($file).${NC}"
            return 1
        fi
    fi
    return 0
}

verify_md5_checksum_file() {
    local archive_file="$1"
    local md5_file="$2"

    if [[ "$DRY_RUN" -eq 1 ]]; then
        if have_cmd md5sum; then
            cmd md5sum -c "$md5_file"
        else
            cmd sh -c "md5 -q '$archive_file'"
        fi
        return 0
    fi

    require_nonempty_file "$archive_file" "$archive_file" || return 1
    require_nonempty_file "$md5_file" "$md5_file" || return 1

    if have_cmd md5sum; then
        if md5sum -c "$md5_file"; then
            return 0
        fi
        echo -e "${RED}  ✗ Checksum verification failed for $archive_file using md5sum.${NC}"
        return 1
    fi

    if have_cmd md5; then
        local expected actual
        expected="$(awk '{print $1}' "$md5_file" | head -n1)"
        actual="$(md5 -q "$archive_file")"
        if [[ -n "$expected" && "$expected" == "$actual" ]]; then
            echo "${archive_file}: OK"
            return 0
        fi
        echo -e "${RED}  ✗ Checksum verification failed for $archive_file using md5.${NC}"
        echo "    expected: $expected"
        echo "    actual:   $actual"
        return 1
    fi

    echo -e "${RED}  ✗ Neither md5sum nor md5 is available; cannot validate $archive_file.${NC}"
    return 1
}

file_nonempty() {
    local file="$1"
    [[ -f "$file" && -s "$file" ]]
}

diamond_index_ready() {
    local db_prefix="$1"
    [[ -f "${db_prefix}.dmnd" ]]
}

blast_index_ready() {
    local db_prefix="$1"
    [[ -f "${db_prefix}.pin" || -f "${db_prefix}.phr" || -f "${db_prefix}.psq" ]]
}

is_fasta_like() {
    local fasta_file="$1"
    [[ -s "$fasta_file" ]] || return 1
    head -n 1 "$fasta_file" 2>/dev/null | grep -q '^>'
}

sanitize_fasta_sequences() {
    local input="$1"
    local output="$2"
    local label="$3"
    local stats_file="${output}.sanitize.stats"

    require_nonempty_file "$input" "$label" || return 1

    if [[ "$DRY_RUN" -eq 1 ]]; then
        cmd awk 'BEGIN{removed=0; changed=0} /^>/{print; next} {orig=$0; gsub(/[^A-Za-z]/, "", $0); if ($0 != orig) changed++; removed += (length(orig)-length($0)); if (length($0)>0) print} END{printf("%d\n%d\n", removed, changed)}' "$input" ">" "$output" "2>" "$stats_file"
        return 0
    fi

    awk -v stats="$stats_file" '
        BEGIN { removed=0; changed=0 }
        /^>/ { print; next }
        {
            orig=$0
            gsub(/[^A-Za-z]/, "", $0)
            if ($0 != orig) changed++
            removed += (length(orig)-length($0))
            if (length($0) > 0) print
        }
        END {
            printf("%d\n%d\n", removed, changed) > stats
        }
    ' "$input" > "$output"

    if [[ ! -s "$output" ]]; then
        echo -e "${RED}  ✗ Failed to sanitize FASTA for ${label}: output is empty (${output}).${NC}"
        return 1
    fi

    if [[ -f "$stats_file" ]]; then
        local removed changed
        removed="$(sed -n '1p' "$stats_file" 2>/dev/null || echo 0)"
        changed="$(sed -n '2p' "$stats_file" 2>/dev/null || echo 0)"
        rm -f "$stats_file"
        if [[ "${removed:-0}" -gt 0 ]]; then
            echo -e "${YELLOW}  ! Sanitized ${label}: removed ${removed} non-letter characters across ${changed} sequence lines.${NC}"
        fi
    fi
}

# ── COG ─────────────────────────────────────────────────────────────────────────
setup_cog() {
    echo -e "${YELLOW}Setting up COG database...${NC}"
    cmd mkdir -p "$DB_DIR/cog"
    [[ "$DRY_RUN" -eq 0 ]] && cd "$DB_DIR/cog"

    if [[ "$DRY_RUN" -eq 0 ]] \
        && file_nonempty cog-24.def.tab \
        && file_nonempty cog-24.fun.tab \
        && file_nonempty cddid.tbl.gz \
        && file_nonempty Cog_LE.tar.gz; then
        echo -e "${GREEN}  ✓ Reusing existing COG files (download skipped)${NC}"
        if [[ ! -f COG.dmnd ]]; then
            if [[ ! -f COGorg24.faa ]]; then
                download_optional "$COG_BASE/COGorg24.faa.gz"
                gunzip_if_present "COGorg24.faa.gz"
            fi
            if [[ -f COGorg24.faa ]]; then
                run_indexer "cog" diamond makedb --in COGorg24.faa -d COG
            else
                echo -e "${YELLOW}  ! COGorg24.faa unavailable; COG.dmnd was not rebuilt.${NC}"
            fi
        fi
        echo "  COGclassifier will use these local resources via --download_dir db/cog."
        echo -e "${GREEN}  ✓ COG database ready${NC}"
        return 0
    fi

    # Official COGclassifier workflow resources (COG + CDD)
    download "$COG_BASE/cog-24.def.tab"
    download "$COG_BASE/cog-24.fun.tab"
    download "$CDD_BASE/cddid.tbl.gz"
    download "$CDD_BASE/little_endian/Cog_LE.tar.gz"

    # Useful metadata tables for downstream enrichment (optional)
    download_optional "$COG_BASE/cog-24.cog.csv"
    download_optional "$COG_BASE/cog-24.mapping.tab"
    download_optional "$COG_BASE/cog-24.org.csv"
    download_optional "$COG_BASE/cog-24.pathways.tab"
    download_optional "$COG_BASE/cog-24.tax.csv"

    # Compatibility artifact used by some downstream wrappers.
    if [[ "$DRY_RUN" -eq 1 || ! -f COG.dmnd ]]; then
        download_optional "$COG_BASE/COGorg24.faa.gz"
        gunzip_if_present "COGorg24.faa.gz"
        if [[ "$DRY_RUN" -eq 1 || -f COGorg24.faa ]]; then
            run_indexer "cog" diamond makedb --in COGorg24.faa -d COG
        fi
    fi

    echo "  COGclassifier will use these local resources via --download_dir db/cog."
    echo -e "${GREEN}  ✓ COG database ready${NC}"
}

# ── Pfam ────────────────────────────────────────────────────────────────────────
setup_pfam() {
    echo -e "${YELLOW}Setting up Pfam database...${NC}"
    cmd mkdir -p "$DB_DIR/pfam"
    [[ "$DRY_RUN" -eq 0 ]] && cd "$DB_DIR/pfam"

    if [[ "$DRY_RUN" -eq 0 ]] && file_nonempty Pfam-A.hmm; then
        echo -e "${GREEN}  ✓ Reusing existing Pfam HMM file (download skipped)${NC}"
        if ! run_indexer "pfam" hmmpress Pfam-A.hmm; then
            echo -e "${YELLOW}  ! Existing Pfam files failed during hmmpress; re-downloading.${NC}"
            rm -f Pfam-A.hmm Pfam-A.hmm.gz Pfam-A.hmm.h3f Pfam-A.hmm.h3i Pfam-A.hmm.h3m Pfam-A.hmm.h3p
            download "$PFAM_BASE/Pfam-A.hmm.gz"
            gunzip_if_present "Pfam-A.hmm.gz"
            run_indexer "pfam" hmmpress Pfam-A.hmm
        fi
    else
        download "$PFAM_BASE/Pfam-A.hmm.gz"
        gunzip_if_present "Pfam-A.hmm.gz"
        run_indexer "pfam" hmmpress Pfam-A.hmm
    fi

    if [[ "$DRY_RUN" -eq 1 ]] || ! file_nonempty Pfam-A.hmm.dat; then
        download_optional "$PFAM_BASE/Pfam-A.hmm.dat.gz"
        gunzip_if_present "Pfam-A.hmm.dat.gz"
    fi

    # Local metadata placeholder matching existing structure
    if [[ "$DRY_RUN" -eq 1 ]]; then
        cmd cp Pfam-A.hmm.dat complete_pfam_definitions.txt
    elif [[ -f Pfam-A.hmm.dat && ! -f complete_pfam_definitions.txt ]]; then
        cmd cp Pfam-A.hmm.dat complete_pfam_definitions.txt
    fi
    echo -e "${GREEN}  ✓ Pfam database ready${NC}"
}

# ── TIGRfam ─────────────────────────────────────────────────────────────────────
setup_tigrfam() {
    echo -e "${YELLOW}Setting up TIGRfam database...${NC}"
    cmd mkdir -p "$DB_DIR/tigrfam"
    [[ "$DRY_RUN" -eq 0 ]] && cd "$DB_DIR/tigrfam"

    if [[ "$DRY_RUN" -eq 0 ]] && file_nonempty TIGRFAMs_15.0_HMM.LIB; then
        echo -e "${GREEN}  ✓ Reusing existing TIGRFAM HMM library (download skipped)${NC}"
        if ! run_indexer "tigrfam" hmmpress TIGRFAMs_15.0_HMM.LIB; then
            echo -e "${YELLOW}  ! Existing TIGRFAM files failed during hmmpress; re-downloading.${NC}"
            rm -f TIGRFAMs_15.0_HMM.LIB TIGRFAMs_15.0_HMM.LIB.gz TIGRFAMs_15.0_HMM.LIB.h3f TIGRFAMs_15.0_HMM.LIB.h3i TIGRFAMs_15.0_HMM.LIB.h3m TIGRFAMs_15.0_HMM.LIB.h3p
            download "$TIGRFAM_BASE/TIGRFAMs_15.0_HMM.LIB.gz"
            gunzip_if_present "TIGRFAMs_15.0_HMM.LIB.gz"
            run_indexer "tigrfam" hmmpress TIGRFAMs_15.0_HMM.LIB
        fi
    else
        download "$TIGRFAM_BASE/TIGRFAMs_15.0_HMM.LIB.gz"
        gunzip_if_present "TIGRFAMs_15.0_HMM.LIB.gz"
        run_indexer "tigrfam" hmmpress TIGRFAMs_15.0_HMM.LIB
    fi

    if [[ "$DRY_RUN" -eq 1 || ! -f TIGRFAMS_GO_LINK ]]; then
        download_optional "$TIGRFAM_BASE/TIGRFAMS_GO_LINK"
    fi
    if [[ "$DRY_RUN" -eq 1 || ! -f TIGRFAMS_ROLE_LINK ]]; then
        download_optional "$TIGRFAM_BASE/TIGRFAMS_ROLE_LINK"
    fi
    if [[ "$DRY_RUN" -eq 1 || ! -f README.md ]]; then
        download_optional "$PGAP_BASE/README.md"
    fi
    echo -e "${GREEN}  ✓ TIGRfam database ready${NC}"
}

# ── dbCAN ───────────────────────────────────────────────────────────────────────
setup_dbcan() {
    echo -e "${YELLOW}Setting up dbCAN database...${NC}"
    cmd mkdir -p "$DB_DIR/dbcan"
    cmd mkdir -p "$DB_DIR/dbcan/tools"
    [[ "$DRY_RUN" -eq 0 ]] && cd "$DB_DIR/dbcan"

    if [[ "$DRY_RUN" -eq 0 ]]; then
        local existing_hmm_db=""
        if dbcan_assets_ready; then
            existing_hmm_db="$(dbcan_pick_hmm_db || true)"
            dbcan_try_press_if_missing "$existing_hmm_db"
            if [[ ! -f dbCAN.txt && -f dbCAN-HMMdb-V14.txt ]]; then
                cp dbCAN-HMMdb-V14.txt dbCAN.txt
            fi
            echo -e "${GREEN}  ✓ Reusing existing dbCAN/CAZy assets from $DB_DIR/dbcan (download skipped)${NC}"
            echo -e "${GREEN}  ✓ dbCAN database ready${NC}"
            return 0
        fi
    fi

    local cpus
    cpus="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)"

    if ! dbcan_build_official "$cpus"; then
        dbcan_build_fallback || {
            echo -e "${RED}  dbCAN setup failed in both official and fallback modes.${NC}"
            return 1
        }
    fi

    if [[ "$DRY_RUN" -eq 0 && ! -f dbCAN.txt && ! -f dbCAN-HMMdb-V12.txt ]]; then
        echo -e "${RED}  dbCAN setup completed but required HMM database file is missing.${NC}"
        return 1
    fi

    echo -e "${GREEN}  ✓ dbCAN database built with official dbCAN runner (includes dbCAN3 datasets and CGC/PUL-compatible assets).${NC}"
    echo -e "${GREEN}  ✓ dbCAN database ready${NC}"
}

# ── KEGG (KofamScan) ────────────────────────────────────────────────────────────
setup_kegg() {
    echo -e "${YELLOW}Setting up KEGG/KofamScan databases...${NC}"
    cmd mkdir -p "$DB_DIR/kegg"
    [[ "$DRY_RUN" -eq 0 ]] && cd "$DB_DIR/kegg"

    if [[ "$DRY_RUN" -eq 0 && -d profiles ]] && file_nonempty ko_list; then
        echo -e "${GREEN}  ✓ Reusing existing KOfam profiles and ko_list (download skipped)${NC}"
    else
        download_from_mirrors "profiles.tar.gz" \
            "$KOFAM_BASE/profiles.tar.gz" \
            "$KOFAM_FALLBACK_BASE/profiles.tar.gz"
        download_from_mirrors "ko_list.gz" \
            "$KOFAM_BASE/ko_list.gz" \
            "$KOFAM_FALLBACK_BASE/ko_list.gz"
        cmd tar -xzf profiles.tar.gz
        gunzip_if_present ko_list.gz
    fi
    echo -e "${GREEN}  ✓ KEGG database ready${NC}"
}

# ── EggNOG ──────────────────────────────────────────────────────────────────────
setup_eggnog() {
    echo -e "${YELLOW}Setting up EggNOG database...${NC}"
    cmd mkdir -p "$DB_DIR/eggnog"
    [[ "$DRY_RUN" -eq 0 ]] && cd "$DB_DIR/eggnog"

    echo "  Method: direct emapperdb download (stable mirror path)"
    echo "  Source: $EGGNOG_EMAPPER_BASE"
    echo "  Target: $DB_DIR/eggnog"
    echo "  Scope:  eggnog.db + eggnog_proteins.dmnd + taxonomy bundle"
    echo ""

    # Remove broken zero-byte artifacts from interrupted/failed downloads.
    if [[ "$DRY_RUN" -eq 0 ]]; then
        local maybe_partial
        for maybe_partial in eggnog.db.gz eggnog_proteins.dmnd.gz eggnog.taxa.tar.gz; do
            if [[ -f "$maybe_partial" && ! -s "$maybe_partial" ]]; then
                rm -f "$maybe_partial"
            fi
        done
    fi

    if [[ "$DRY_RUN" -eq 1 ]]; then
        download "$EGGNOG_EMAPPER_BASE/eggnog.db.gz" "eggnog.db.gz"
        download "$EGGNOG_EMAPPER_BASE/eggnog_proteins.dmnd.gz" "eggnog_proteins.dmnd.gz"
        download "$EGGNOG_EMAPPER_BASE/eggnog.taxa.tar.gz" "eggnog.taxa.tar.gz"
        cmd gunzip -kf eggnog.db.gz
        cmd gunzip -kf eggnog_proteins.dmnd.gz
        cmd tar -xzf eggnog.taxa.tar.gz
        echo -e "${GREEN}  ✓ EggNOG database ready${NC}"
        return 0
    fi

    if file_nonempty eggnog.db && file_nonempty eggnog_proteins.dmnd; then
        if require_min_size_bytes "eggnog.db" 5000000000 "eggnog.db" \
            && require_min_size_bytes "eggnog_proteins.dmnd" 3000000000 "eggnog_proteins.dmnd"; then
            echo -e "${GREEN}  ✓ Reusing existing EggNOG DB assets (download skipped)${NC}"
            echo -e "${GREEN}  ✓ EggNOG database ready${NC}"
            return 0
        fi
        echo -e "${YELLOW}  ! Existing EggNOG assets look incomplete; re-downloading.${NC}"
        rm -f eggnog.db eggnog_proteins.dmnd
    fi

    download "$EGGNOG_EMAPPER_BASE/eggnog.db.gz" "eggnog.db.gz"
    download "$EGGNOG_EMAPPER_BASE/eggnog_proteins.dmnd.gz" "eggnog_proteins.dmnd.gz"
    download "$EGGNOG_EMAPPER_BASE/eggnog.taxa.tar.gz" "eggnog.taxa.tar.gz"

    # Validate large downloads before extraction so partial files are caught early.
    require_min_size_bytes "eggnog.db.gz" 1000000000 "eggnog.db.gz" || return 1
    require_min_size_bytes "eggnog_proteins.dmnd.gz" 1000000000 "eggnog_proteins.dmnd.gz" || return 1
    require_min_size_bytes "eggnog.taxa.tar.gz" 10000000 "eggnog.taxa.tar.gz" || return 1
    validate_gzip_file "eggnog.db.gz" "eggnog.db.gz" || return 1
    validate_gzip_file "eggnog_proteins.dmnd.gz" "eggnog_proteins.dmnd.gz" || return 1
    validate_gzip_file "eggnog.taxa.tar.gz" "eggnog.taxa.tar.gz" || return 1

    gunzip_if_present eggnog.db.gz
    gunzip_if_present eggnog_proteins.dmnd.gz
    if [[ -f eggnog.taxa.tar.gz ]]; then
        cmd tar -xzf eggnog.taxa.tar.gz
    fi

    require_min_size_bytes "eggnog.db" 5000000000 "eggnog.db" || return 1
    require_min_size_bytes "eggnog_proteins.dmnd" 3000000000 "eggnog_proteins.dmnd" || return 1

    if [[ ! -f eggnog.db ]]; then
        echo -e "${RED}  ✗ EggNOG setup failed: missing eggnog.db.${NC}"
        return 1
    fi
    if [[ ! -f eggnog_proteins.dmnd ]]; then
        echo -e "${RED}  ✗ EggNOG setup failed: missing eggnog_proteins.dmnd.${NC}"
        return 1
    fi

    echo -e "${GREEN}  ✓ EggNOG database ready${NC}"
    return 0
}

# ── MEROPS ──────────────────────────────────────────────────────────────────────
setup_merops() {
    echo -e "${YELLOW}Setting up MEROPS database...${NC}"
    cmd mkdir -p "$DB_DIR/merops"
    [[ "$DRY_RUN" -eq 0 ]] && cd "$DB_DIR/merops"

    echo "  Source: $MEROPS_BASE"
    echo "  Method note: MEROPS papers describe BLAST/HMMER-based homologue discovery (Rawlings et al., 2014; 2018)."

    if [[ "$DRY_RUN" -eq 0 ]] && file_nonempty pepunit.lib && file_nonempty meropsscan.lib; then
        echo -e "${GREEN}  ✓ Reusing existing MEROPS source files (download skipped)${NC}"
    else
        download "$MEROPS_BASE/pepunit.lib"
        download "$MEROPS_BASE/meropsscan.lib"
    fi

    if [[ "$DRY_RUN" -eq 1 ]]; then
        cmd cp pepunit.lib pepunit_raw.lib
        cmd awk 'BEGIN{removed=0; changed=0} /^>/{print; next} {orig=$0; gsub(/[^A-Za-z]/, "", $0); if ($0 != orig) changed++; removed += (length(orig)-length($0)); if (length($0)>0) print} END{if (removed>0) printf("Sanitized MEROPS pepunit.lib: removed %d non-letter characters across %d sequence lines\\n", removed, changed) > "/dev/stderr"}' pepunit_raw.lib ">" pepunit.lib
        cmd cp pepunit.lib pepunit.txt
        cmd diamond makedb --in pepunit.lib -d merops
        cmd hmmpress meropsscan.lib
        echo -e "${GREEN}  ✓ MEROPS database ready${NC}"
    elif [[ -f pepunit.lib && -f meropsscan.lib ]]; then
        local merops_hmmpress_ok=1
        cmd cp pepunit.lib pepunit_raw.lib
        sanitize_fasta_sequences pepunit_raw.lib pepunit.lib "MEROPS pepunit.lib" || return 1
        cmd cp pepunit.lib pepunit.txt
        if ! run_indexer "merops" diamond makedb --in pepunit.lib -d merops; then
            echo -e "${YELLOW}  ! Existing MEROPS FASTA failed during DIAMOND indexing; re-downloading source files.${NC}"
            rm -f pepunit.lib meropsscan.lib merops.dmnd
            download "$MEROPS_BASE/pepunit.lib"
            download "$MEROPS_BASE/meropsscan.lib"
            cmd cp pepunit.lib pepunit_raw.lib
            sanitize_fasta_sequences pepunit_raw.lib pepunit.lib "MEROPS pepunit.lib" || return 1
            run_indexer "merops" diamond makedb --in pepunit.lib -d merops
        fi
        if ! run_indexer "merops" hmmpress meropsscan.lib; then
            merops_hmmpress_ok=0
            echo -e "${YELLOW}  ! MEROPS HMMER index creation failed; continuing with DIAMOND index only.${NC}"
        fi
        if [[ ! -f merops.dmnd ]]; then
            echo -e "${RED}  MEROPS DIAMOND index not created (merops.dmnd).${NC}"
            return 1
        fi
        if [[ "$merops_hmmpress_ok" -eq 1 && ! -f meropsscan.lib.h3m ]]; then
            echo -e "${RED}  MEROPS HMMER index not created (meropsscan.lib.h3m).${NC}"
            return 1
        fi
        if [[ "$merops_hmmpress_ok" -eq 0 ]]; then
            echo -e "${YELLOW}  ! meropsscan.lib was downloaded but not pressed (.h3* files missing).${NC}"
        fi
        echo -e "${GREEN}  ✓ MEROPS database ready${NC}"
    else
        echo -e "${RED}  MEROPS required files were not found after download (pepunit.lib, meropsscan.lib).${NC}"
        return 1
    fi
}

# ── TCDB ────────────────────────────────────────────────────────────────────────
setup_tcdb() {
    echo -e "${YELLOW}Setting up TCDB database...${NC}"
    cmd mkdir -p "$DB_DIR/tcdb"
    [[ "$DRY_RUN" -eq 0 ]] && cd "$DB_DIR/tcdb"

    echo "  Source: https://www.tcdb.org"
    echo "  Method note: Official local CLI path uses NCBI BLAST+ against TCDB FASTA (Saier et al., 2006; Saier et al., 2021)."

    # Primary link from provided document
    if [[ "$DRY_RUN" -eq 1 || ! -f tcdb_download.html ]]; then
        download_optional "$TCDB_URL" "tcdb_download.html"
    fi

    # Direct FASTA location used by existing workflows
    if [[ "$DRY_RUN" -eq 0 ]] && is_fasta_like tcdb.fasta; then
        echo -e "${GREEN}  ✓ Reusing existing TCDB FASTA (download skipped)${NC}"
    else
        download "https://www.tcdb.org/public/tcdb" "tcdb.fasta"
    fi
    if [[ "$DRY_RUN" -eq 1 ]]; then
        cmd makeblastdb -in tcdb.fasta -dbtype prot -out tcdb_blast
        echo -e "${GREEN}  ✓ TCDB database ready${NC}"
        return
    fi

    if [[ ! -s tcdb.fasta ]]; then
        echo -e "${RED}  TCDB FASTA could not be downloaded automatically.${NC}"
        echo "  Please download TCDB protein FASTA from https://www.tcdb.org/download.php"
        return 1
    fi

    if head -n 1 tcdb.fasta 2>/dev/null | grep -qi '<!doctype\|<html'; then
        echo -e "${RED}  TCDB download returned HTML instead of FASTA. Check access policy/network.${NC}"
        return 1
    fi

    # Official local CLI-style annotation path: NCBI BLAST+ against TCDB FASTA
    if blast_index_ready tcdb_blast; then
        echo -e "${GREEN}  ✓ Reusing existing TCDB BLAST index (formatting skipped)${NC}"
    elif ! run_indexer "tcdb" makeblastdb -in tcdb.fasta -dbtype prot -out tcdb_blast; then
        echo -e "${YELLOW}  ! Existing TCDB FASTA failed during BLAST formatting; re-downloading.${NC}"
        rm -f tcdb.fasta tcdb_blast.p* tcdb_blast.n* tcdb_blast.*
        download "https://www.tcdb.org/public/tcdb" "tcdb.fasta"
        run_indexer "tcdb" makeblastdb -in tcdb.fasta -dbtype prot -out tcdb_blast
    fi

    # Compatibility artifact used by DIAMOND-based TCDB wrappers.
    if diamond_index_ready tcdb; then
        echo -e "${GREEN}  ✓ Reusing existing TCDB DIAMOND index (formatting skipped)${NC}"
    else
        run_indexer "tcdb" diamond makedb --in tcdb.fasta -d tcdb
    fi
    echo -e "${GREEN}  ✓ TCDB database ready${NC}"
}

# ── UniProt ─────────────────────────────────────────────────────────────────────
setup_uniprot() {
    echo -e "${YELLOW}Setting up UniProt database...${NC}"
    cmd mkdir -p "$DB_DIR/uniprot"
    [[ "$DRY_RUN" -eq 0 ]] && cd "$DB_DIR/uniprot"

    if [[ "$DRY_RUN" -eq 0 ]] && is_fasta_like uniprot_sprot.fasta; then
        echo -e "${GREEN}  ✓ Reusing existing UniProt FASTA (download skipped)${NC}"
    else
        download "$UNIPROT_SPROT_URL"
        gunzip_if_present uniprot_sprot.fasta.gz
    fi

    if diamond_index_ready uniprot_sprot; then
        echo -e "${GREEN}  ✓ Reusing existing UniProt DIAMOND index (formatting skipped)${NC}"
    elif ! run_indexer "uniprot" diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot; then
        echo -e "${YELLOW}  ! Existing UniProt FASTA failed during DIAMOND formatting; re-downloading.${NC}"
        rm -f uniprot_sprot.fasta uniprot_sprot.fasta.gz uniprot_sprot.dmnd
        download "$UNIPROT_SPROT_URL"
        gunzip_if_present uniprot_sprot.fasta.gz
        run_indexer "uniprot" diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot
    fi
    echo "  Using official Swiss-Prot resources for UPIMAPI-backed annotation."
    echo -e "${GREEN}  ✓ UniProt Swiss-Prot database ready${NC}"
}

# ── InterPro ────────────────────────────────────────────────────────────────────
setup_interpro() {
    echo -e "${YELLOW}Setting up InterPro resources (InterProScan + mapping files)...${NC}"
    cmd mkdir -p "$DB_DIR/interpro"
    [[ "$DRY_RUN" -eq 0 ]] && cd "$DB_DIR/interpro"

    # Signature mapping used by existing workflow outputs
    if [[ "$DRY_RUN" -eq 0 ]] && { file_nonempty signature_to_ipr.tsv || file_nonempty signature_to_ipr.dat; }; then
        echo -e "${GREEN}  ✓ Reusing existing InterPro signature mappings (download skipped)${NC}"
    else
        download_optional "$INTERPRO_BASE/current_release/entry.list"
        download_optional "$INTERPRO_BASE/current_release/interpro2go"
        download_optional "$INTERPRO_BASE/current_release/signature_to_ipr.dat.gz"
        if [[ "$DRY_RUN" -eq 1 ]]; then
            cmd gunzip -kf signature_to_ipr.dat.gz
            cmd cp signature_to_ipr.dat signature_to_ipr.tsv
        elif [[ -f signature_to_ipr.dat.gz ]]; then
            gunzip_if_present signature_to_ipr.dat.gz
            cp signature_to_ipr.dat signature_to_ipr.tsv || true
        fi
    fi

    # Official InterProScan core package + checksum verification.
    if [[ "$DRY_RUN" -eq 1 || ! -f "$INTERPROSCAN_ARCHIVE" ]]; then
        download "$INTERPROSCAN_SOFTWARE_BASE/$INTERPROSCAN_ARCHIVE" "$INTERPROSCAN_ARCHIVE"
    fi
    if [[ "$DRY_RUN" -eq 1 || ! -f "$INTERPROSCAN_MD5_FILE" ]]; then
        download "$INTERPROSCAN_SOFTWARE_BASE/$INTERPROSCAN_MD5_FILE" "$INTERPROSCAN_MD5_FILE"
    fi

    verify_md5_checksum_file "$INTERPROSCAN_ARCHIVE" "$INTERPROSCAN_MD5_FILE"

    # Extract exactly as recommended in InterProScan docs (-p preserves permissions).
    if [[ "$DRY_RUN" -eq 1 ]]; then
        cmd tar -pxvzf "$INTERPROSCAN_ARCHIVE"
    else
        if [[ ! -d "interproscan-${INTERPROSCAN_VERSION}" ]]; then
            tar -pxvzf "$INTERPROSCAN_ARCHIVE"
        else
            echo -e "${GREEN}  ✓ Reusing existing interproscan-${INTERPROSCAN_VERSION} directory (extract skipped)${NC}"
        fi
    fi

    if [[ "$DRY_RUN" -eq 0 ]]; then
        require_nonempty_file "interproscan-${INTERPROSCAN_VERSION}/interproscan.sh" "InterProScan launcher" || return 1
    fi

    echo -e "${GREEN}  ✓ InterPro setup complete${NC}"
}

# ── RASTtk ──────────────────────────────────────────────────────────────────────
setup_rasttk() {
    echo -e "${YELLOW}Setting up RASTtk SEED enrichment assets...${NC}"
    cmd mkdir -p "$DB_DIR/rasttk"

    local source_dir="${MARGIE_RASTTK_DB_SOURCE:-}"
    if [[ -z "$source_dir" ]]; then
        for candidate in \
            "$HOME/Desktop/3-13-26/context-based-annotation/db/rasttk" \
            "$SCRIPT_DIR/../3-13-26/context-based-annotation/db/rasttk" \
            "$SCRIPT_DIR/../context-based-annotation/db/rasttk"; do
            if [[ -d "$candidate" ]]; then
                source_dir="$candidate"
                break
            fi
        done
    fi

    if [[ -n "$source_dir" && -d "$source_dir" ]]; then
        echo "  Copying RASTtk enrichment files from: $source_dir"
        if [[ "$DRY_RUN" -eq 1 ]]; then
            cmd cp "$source_dir/subsystem_mapping.tsv" "$DB_DIR/rasttk/"
            cmd cp -R "$source_dir/variant_definitions" "$DB_DIR/rasttk/"
        else
            [[ -f "$source_dir/subsystem_mapping.tsv" ]] && cp "$source_dir/subsystem_mapping.tsv" "$DB_DIR/rasttk/"
            [[ -d "$source_dir/variant_definitions" ]] && rm -rf "$DB_DIR/rasttk/variant_definitions" && cp -R "$source_dir/variant_definitions" "$DB_DIR/rasttk/"
        fi
    fi

    if [[ -f "$DB_DIR/rasttk/subsystem_mapping.tsv" && -d "$DB_DIR/rasttk/variant_definitions" ]]; then
        echo -e "${GREEN}  ✓ RASTtk SEED enrichment assets ready${NC}"
    else
        echo -e "${YELLOW}  ! RASTtk core annotation will work, but full SEED enrichment requires:${NC}"
        echo "    - db/rasttk/subsystem_mapping.tsv"
        echo "    - db/rasttk/variant_definitions/"
        echo "    You can provide them by setting MARGIE_RASTTK_DB_SOURCE=/path/to/db/rasttk and rerunning ./setup_databases.sh --rasttk"
    fi
}

# ── Dispatcher ──────────────────────────────────────────────────────────────────

ALL=(cog pfam tigrfam dbcan kegg eggnog merops tcdb uniprot interpro rasttk)

is_db_ready_minimal() {
    local db="$1"
    case "$db" in
        cog) [[ -f "$DB_DIR/cog/cog-24.def.tab" && -f "$DB_DIR/cog/cddid.tbl.gz" && -f "$DB_DIR/cog/Cog_LE.tar.gz" ]] ;;
        pfam) [[ -f "$DB_DIR/pfam/Pfam-A.hmm" ]] ;;
        tigrfam) [[ -f "$DB_DIR/tigrfam/TIGRFAMs_15.0_HMM.LIB" ]] ;;
        dbcan) [[ -f "$DB_DIR/dbcan/dbCAN.txt" || -f "$DB_DIR/dbcan/dbCAN-HMMdb-V12.txt" || -f "$DB_DIR/dbcan/dbCAN-HMMdb-V14.txt" ]] ;;
        kegg) [[ -f "$DB_DIR/kegg/ko_list" && -d "$DB_DIR/kegg/profiles" ]] ;;
        eggnog) [[ -f "$DB_DIR/eggnog/eggnog.db" && -f "$DB_DIR/eggnog/eggnog_proteins.dmnd" ]] ;;
        merops) [[ -f "$DB_DIR/merops/pepunit.lib" && -f "$DB_DIR/merops/merops.dmnd" ]] ;;
        tcdb) [[ -f "$DB_DIR/tcdb/tcdb.fasta" && ( -f "$DB_DIR/tcdb/tcdb_blast.pin" || -f "$DB_DIR/tcdb/tcdb_blast.phr" ) ]] ;;
        uniprot) [[ -f "$DB_DIR/uniprot/uniprot_sprot.fasta" ]] ;;
        interpro) [[ ( -f "$DB_DIR/interpro/signature_to_ipr.tsv" || -f "$DB_DIR/interpro/signature_to_ipr.dat" ) && -f "$DB_DIR/interpro/interproscan-${INTERPROSCAN_VERSION}/interproscan.sh" ]] ;;
        rasttk) [[ -f "$DB_DIR/rasttk/subsystem_mapping.tsv" && -d "$DB_DIR/rasttk/variant_definitions" ]] ;;
        *) return 1 ;;
    esac
}

all_dbs_ready_minimal() {
    local db
    for db in "${ALL[@]}"; do
        if ! is_db_ready_minimal "$db"; then
            echo "$db"
            return 1
        fi
    done
    return 0
}

if [[ "$SELECTED" == "all" ]]; then
    TARGETS=("${ALL[@]}")
else
    IFS=',' read -ra TARGETS <<< "$SELECTED"
fi

if [[ "$DRY_RUN" -eq 0 && "$SELECTED" == "all" && "${MARGIE_FORCE_SETUP:-0}" != "1" && -f "$STATE_DATABASES_FILE" ]]; then
    if missing_db="$(all_dbs_ready_minimal)"; then
        echo -e "${GREEN}Database setup already recorded as complete; skipping full database setup.${NC}"
        echo "Marker file: $STATE_DATABASES_FILE"
        echo "Set MARGIE_FORCE_SETUP=1 to force re-download and re-index."
        exit 0
    else
        echo -e "${YELLOW}Database state marker exists, but '${missing_db}' is incomplete. Continuing setup.${NC}"
    fi
fi

first_time_guidance

for db in "${TARGETS[@]}"; do
    case "$db" in
        cog)        setup_cog ;;
        pfam)       setup_pfam ;;
        tigrfam)    setup_tigrfam ;;
        dbcan)      setup_dbcan ;;
        kegg)       setup_kegg ;;
        eggnog)     setup_eggnog ;;
        merops)     setup_merops ;;
        tcdb)       setup_tcdb ;;
        uniprot)    setup_uniprot ;;
        interpro)   setup_interpro ;;
        rasttk)     setup_rasttk ;;
        *)          echo -e "${RED}Unknown database: $db${NC}" ;;
    esac
    echo ""
done

echo -e "${GREEN}=== Database setup complete ===${NC}"
if [[ "$DRY_RUN" -eq 0 && "$SELECTED" == "all" ]]; then
    write_database_state
    echo "State recorded: $STATE_DATABASES_FILE"
fi
echo "Next step: run ./setup_containers.sh --skip-db"
