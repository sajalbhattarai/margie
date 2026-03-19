#!/usr/bin/env bash
# setup_containers.sh
# Auto-detects Docker or Apptainer and builds all pipeline containers.
# Usage:
#   ./setup_containers.sh           # auto-detect + download all databases
#   ./setup_containers.sh --docker  # force Docker
#   ./setup_containers.sh --apptainer  # force Apptainer
#   ./setup_containers.sh --skip-db     # build containers only
#   ./setup_containers.sh --operon      # build only the operon container path
#   ./setup_containers.sh --prodigal --rasttk   # build only selected tool containers
#   ./setup_containers.sh --db cog,pfam # build containers + selected databases
#   ./setup_containers.sh --cog --pfam   # build only COG/Pfam container paths and set up those DBs
#   ./setup_containers.sh --all | -all   # all supported databases
#   ./setup_containers.sh --platform auto|linux/amd64|linux/arm64
#   ./setup_containers.sh --dry-run     # preview build/download commands only
#   ./setup_containers.sh --preflight   # check environment and exit
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINERS_DIR="$SCRIPT_DIR/processing/containers"
STATE_DIR="$SCRIPT_DIR/.margie_state"
STATE_CONTAINERS_FILE="$STATE_DIR/containers.ready"
STATE_DATABASES_FILE="$STATE_DIR/databases.ready"

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'

CONTAINERS=(
    "cog:cog-annotation:1.0"
    "kegg:kegg-annotation:1.0"
    "prodigal:prodigal-annotation:1.0"
    "operon:operon-annotation:1.0"
    "tigrfam:tigrfam-annotation:1.0"
    "consolidation:consolidation-annotation:1.0"
    "merops:merops-annotation:1.0"
    "rasttk:rasttk-annotation:1.0"
    "uniprot:uniprot-annotation:1.0"
    "interpro:interpro-annotation:1.0"
    "pfam:pfam-annotation:1.0"
    "tcdb:tcdb-annotation:1.0"
)

OFFICIAL_DBCAN_IMAGE="haidyi/run_dbcan:latest"
OFFICIAL_EGGNOG_IMAGE="quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0"
OFFICIAL_EGGNOG_PLATFORM="linux/amd64"
OFFICIAL_OPERON_IMAGE="ghcr.io/sajalbhattarai/for_oprn:latest"

# Optional prebuilt images hosted in GitHub Container Registry (GHCR).
# These can be overridden via environment variables for forks/custom orgs.
GHCR_OWNER="${MARGIE_GHCR_OWNER:-sajalbhattarai}"
GHCR_REPO="${MARGIE_GHCR_REPO:-margie}"
PREBUILT_PRODIGAL_IMAGE="${MARGIE_GHCR_PRODIGAL_IMAGE:-ghcr.io/${GHCR_OWNER}/${GHCR_REPO}/prodigal-annotation:1.0}"
PREBUILT_RASTTK_IMAGE="${MARGIE_GHCR_RASTTK_IMAGE:-ghcr.io/${GHCR_OWNER}/${GHCR_REPO}/rasttk-annotation:1.0}"
DB_BACKED_TOOLS=(cog pfam tigrfam dbcan kegg eggnog merops tcdb uniprot interpro rasttk)

# ── Runtime detection ───────────────────────────────────────────────────────────
MODE=""
SETUP_DB=1
DB_SELECTION="all"
TOOL_SELECTION="all"
DRY_RUN=0
PREFLIGHT_ONLY=0
DOCKER_PLATFORM="auto"
FORCE_ALL_DBS=0
REQUESTED_DBS=()
REQUESTED_TOOLS=()
BUILD_GUI=0

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

append_requested_tool() {
    local tool_name="$1"
    local existing
    for existing in "${REQUESTED_TOOLS[@]-}"; do
        if [[ "$existing" == "$tool_name" ]]; then
            return 0
        fi
    done
    REQUESTED_TOOLS+=("$tool_name")
}

is_db_backed_tool() {
    local tool_name="$1"
    local existing
    for existing in "${DB_BACKED_TOOLS[@]}"; do
        if [[ "$existing" == "$tool_name" ]]; then
            return 0
        fi
    done
    return 1
}

append_requested_tool_with_db() {
    local tool_name="$1"
    append_requested_tool "$tool_name"
    if is_db_backed_tool "$tool_name"; then
        append_requested_db "$tool_name"
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --docker)
            MODE="docker"
            shift
            ;;
        --apptainer)
            MODE="apptainer"
            shift
            ;;
        --skip-db)
            SETUP_DB=0
            shift
            ;;
        --db)
            DB_SELECTION="${2:-all}"
            if [[ "$DB_SELECTION" == "all" ]]; then
                FORCE_ALL_DBS=1
            else
                append_requested_db_list "$DB_SELECTION"
            fi
            shift 2
            ;;
        --all|-all|all)
            FORCE_ALL_DBS=1
            shift
            ;;
        --cog)
            append_requested_tool_with_db "cog"
            shift
            ;;
        --pfam)
            append_requested_tool_with_db "pfam"
            shift
            ;;
        --tigrfam)
            append_requested_tool_with_db "tigrfam"
            shift
            ;;
        --dbcan)
            append_requested_tool_with_db "dbcan"
            shift
            ;;
        --kegg)
            append_requested_tool_with_db "kegg"
            shift
            ;;
        --eggnog)
            append_requested_tool_with_db "eggnog"
            shift
            ;;
        --merops)
            append_requested_tool_with_db "merops"
            shift
            ;;
        --tcdb)
            append_requested_tool_with_db "tcdb"
            shift
            ;;
        --uniprot)
            append_requested_tool_with_db "uniprot"
            shift
            ;;
        --interpro)
            append_requested_tool_with_db "interpro"
            shift
            ;;
        --rasttk)
            append_requested_tool_with_db "rasttk"
            shift
            ;;
        --prodigal)
            append_requested_tool "prodigal"
            shift
            ;;
        --operon)
            append_requested_tool "operon"
            shift
            ;;
        --consolidation)
            append_requested_tool "consolidation"
            shift
            ;;
        --gui)
            BUILD_GUI=1
            shift
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --platform)
            DOCKER_PLATFORM="${2:-auto}"
            shift 2
            ;;
        --preflight)
            PREFLIGHT_ONLY=1
            shift
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            exit 1
            ;;
    esac
done

if [[ "$FORCE_ALL_DBS" -eq 1 ]]; then
    DB_SELECTION="all"
elif [[ "${#REQUESTED_DBS[@]}" -gt 0 ]]; then
    IFS=',' DB_SELECTION="${REQUESTED_DBS[*]}"
fi

if [[ "${#REQUESTED_TOOLS[@]}" -gt 0 ]]; then
    IFS=',' TOOL_SELECTION="${REQUESTED_TOOLS[*]}"
    if [[ "$FORCE_ALL_DBS" -eq 0 && "${#REQUESTED_DBS[@]}" -eq 0 && "$SETUP_DB" -eq 1 ]]; then
        SETUP_DB=0
    fi
fi

if [[ "$DRY_RUN" -eq 0 && -x "$SCRIPT_DIR/install_commands.sh" ]]; then
    "$SCRIPT_DIR/install_commands.sh" --quiet || true
fi

PASS_COUNT=0
WARN_COUNT=0
FAIL_COUNT=0

pass() {
    echo -e "${GREEN}[PASS]${NC} $1"
    ((PASS_COUNT++)) || true
}

warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
    ((WARN_COUNT++)) || true
}

fail() {
    echo -e "${RED}[FAIL]${NC} $1"
    ((FAIL_COUNT++)) || true
}

run_preflight() {
    echo -e "${GREEN}=== Margie Preflight Check ===${NC}"
    echo "Checking required runtime and setup tools"
    echo ""

    if command -v bash >/dev/null 2>&1; then
        pass "bash found"
    else
        fail "bash not found"
    fi

    if command -v python3 >/dev/null 2>&1; then
        pass "python3 found"
    else
        fail "python3 not found"
    fi

    if command -v wget >/dev/null 2>&1; then
        pass "wget found"
    elif command -v python3 >/dev/null 2>&1; then
        pass "wget not found; Python downloader fallback available"
    else
        fail "Neither wget nor python3 found (one is required for database download)"
    fi

    if command -v docker >/dev/null 2>&1; then
        if docker info >/dev/null 2>&1; then
            pass "Docker daemon available"
        else
            warn "Docker installed but daemon not available"
        fi
    else
        warn "Docker not found"
    fi

    if command -v apptainer >/dev/null 2>&1; then
        pass "Apptainer found"
    elif command -v singularity >/dev/null 2>&1; then
        pass "Singularity found"
    else
        warn "Apptainer/Singularity not found"
    fi

    # Host indexers are optional; setup_databases.sh can use built containers.
    if command -v diamond >/dev/null 2>&1; then
        pass "diamond found"
    else
        warn "diamond not found on host (optional; containerized indexing will be used if images exist)"
    fi

    if command -v hmmpress >/dev/null 2>&1; then
        pass "hmmpress found"
    else
        warn "hmmpress not found on host (optional; containerized indexing will be used if images exist)"
    fi

    AVAIL_GB=$(df -Pk "$SCRIPT_DIR" | awk 'NR==2 {printf "%.1f", $4/1024/1024}')
    if awk "BEGIN {exit !($AVAIL_GB >= 120.0)}"; then
        pass "Disk space check passed (${AVAIL_GB} GB available)"
    else
        warn "Low free space (${AVAIL_GB} GB). Recommended: >= 120 GB"
    fi

    echo ""
    echo "Summary: $PASS_COUNT passed, $WARN_COUNT warnings, $FAIL_COUNT failures"

    if [[ "$FAIL_COUNT" -gt 0 ]]; then
        echo ""
        echo -e "${RED}Preflight failed. Install missing requirements before setup.${NC}"
        return 1
    fi

    echo ""
    echo -e "${GREEN}Preflight complete.${NC}"
    echo "If this is a new machine, next run: ./setup_containers.sh"
    echo "After setup, place genome FASTA files into: $SCRIPT_DIR/input"
    return 0
}

if [[ "$PREFLIGHT_ONLY" -eq 1 ]]; then
    run_preflight
    exit $?
fi

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

# Run preflight before setup to provide actionable warnings.
run_preflight || exit 1

if [[ -z "$MODE" ]]; then
    if command -v docker &>/dev/null && docker info &>/dev/null 2>&1; then
        MODE="docker"
    elif command -v apptainer &>/dev/null; then
        MODE="apptainer"
    elif command -v singularity &>/dev/null; then
        MODE="singularity"
    else
        echo -e "${RED}ERROR: Neither Docker nor Apptainer/Singularity found.${NC}"
        echo "  macOS / Linux workstation: install Docker Desktop"
        echo "  HPC: run 'module load apptainer' or 'module load singularity'"
        exit 1
    fi
fi

echo -e "${GREEN}=== Margie Annotation Pipeline — Container Setup ===${NC}"
echo "Runtime: $MODE"
if [[ "$MODE" == "docker" ]]; then
    echo "Docker platform: $DOCKER_PLATFORM_RESOLVED"
fi
echo "Containers directory: $CONTAINERS_DIR"
echo "Container selection: $TOOL_SELECTION"
if [[ "$DRY_RUN" -eq 1 ]]; then
    echo -e "${YELLOW}DRY-RUN: commands will be printed but not executed${NC}"
fi
if [[ "$SETUP_DB" -eq 1 ]]; then
    echo "Database setup: enabled (selection: $DB_SELECTION)"
else
    echo "Database setup: skipped"
fi
echo ""

BUILT=0; FAILED=0

print_cmd() {
    printf '[DRY-RUN]'
    printf ' %q' "$@"
    printf '\n'
}

ensure_state_dir() {
    mkdir -p "$STATE_DIR"
}

write_state_file() {
    local state_file="$1"
    shift
    ensure_state_dir
    {
        printf 'timestamp=%s\n' "$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
        printf 'runtime=%s\n' "$MODE"
        while [[ "$#" -gt 0 ]]; do
            printf '%s\n' "$1"
            shift
        done
    } > "$state_file"
}

mark_containers_ready() {
    write_state_file "$STATE_CONTAINERS_FILE" \
        "tool_selection=$TOOL_SELECTION" \
        "build_gui=$BUILD_GUI" \
        "docker_platform=$DOCKER_PLATFORM_RESOLVED"
}

ALL_DATABASES=(cog pfam tigrfam dbcan kegg eggnog merops tcdb uniprot interpro rasttk)

is_db_ready_minimal() {
    local db="$1"
    case "$db" in
        cog) [[ -f "$SCRIPT_DIR/db/cog/cog-24.def.tab" && -f "$SCRIPT_DIR/db/cog/cddid.tbl.gz" && -f "$SCRIPT_DIR/db/cog/Cog_LE.tar.gz" ]] ;;
        pfam) [[ -f "$SCRIPT_DIR/db/pfam/Pfam-A.hmm" ]] ;;
        tigrfam) [[ -f "$SCRIPT_DIR/db/tigrfam/TIGRFAMs_15.0_HMM.LIB" ]] ;;
        dbcan) [[ -f "$SCRIPT_DIR/db/dbcan/dbCAN.txt" || -f "$SCRIPT_DIR/db/dbcan/dbCAN-HMMdb-V12.txt" || -f "$SCRIPT_DIR/db/dbcan/dbCAN-HMMdb-V14.txt" ]] ;;
        kegg) [[ -f "$SCRIPT_DIR/db/kegg/ko_list" && -d "$SCRIPT_DIR/db/kegg/profiles" ]] ;;
        eggnog) [[ -f "$SCRIPT_DIR/db/eggnog/eggnog.db" && -f "$SCRIPT_DIR/db/eggnog/eggnog_proteins.dmnd" ]] ;;
        merops) [[ -f "$SCRIPT_DIR/db/merops/pepunit.lib" && -f "$SCRIPT_DIR/db/merops/merops.dmnd" ]] ;;
        tcdb) [[ -f "$SCRIPT_DIR/db/tcdb/tcdb.fasta" && ( -f "$SCRIPT_DIR/db/tcdb/tcdb_blast.pin" || -f "$SCRIPT_DIR/db/tcdb/tcdb_blast.phr" ) ]] ;;
        uniprot) [[ -f "$SCRIPT_DIR/db/uniprot/uniprot_sprot.fasta" ]] ;;
        interpro) [[ -f "$SCRIPT_DIR/db/interpro/signature_to_ipr.tsv" || -f "$SCRIPT_DIR/db/interpro/signature_to_ipr.dat" ]] ;;
        rasttk) [[ -f "$SCRIPT_DIR/db/rasttk/subsystem_mapping.tsv" && -d "$SCRIPT_DIR/db/rasttk/variant_definitions" ]] ;;
        *) return 1 ;;
    esac
}

all_dbs_ready_minimal() {
    local db
    for db in "${ALL_DATABASES[@]}"; do
        if ! is_db_ready_minimal "$db"; then
            echo "$db"
            return 1
        fi
    done
    return 0
}

is_default_full_setup_request() {
    [[ "$TOOL_SELECTION" == "all" && "$SETUP_DB" -eq 1 && "$DB_SELECTION" == "all" && "$BUILD_GUI" -eq 0 ]]
}

is_default_container_only_request() {
    [[ "$TOOL_SELECTION" == "all" && "$SETUP_DB" -eq 0 && "$BUILD_GUI" -eq 0 ]]
}

tool_selected() {
    local tool_name="$1"
    if [[ "$TOOL_SELECTION" == "all" ]]; then
        return 0
    fi
    [[ ",${TOOL_SELECTION}," == *",${tool_name},"* ]]
}

if [[ "$DRY_RUN" -eq 0 && "${MARGIE_FORCE_SETUP:-0}" != "1" ]]; then
    if is_default_full_setup_request && [[ -f "$STATE_CONTAINERS_FILE" && -f "$STATE_DATABASES_FILE" ]]; then
        if missing_db="$(all_dbs_ready_minimal)"; then
            echo -e "${GREEN}Setup already recorded as complete; skipping container and database setup.${NC}"
            echo "Marker files: $STATE_CONTAINERS_FILE and $STATE_DATABASES_FILE"
            echo "Set MARGIE_FORCE_SETUP=1 to force a full rebuild."
            exit 0
        else
            echo -e "${YELLOW}State markers found, but database '${missing_db}' is incomplete. Continuing setup.${NC}"
        fi
    fi

    if is_default_container_only_request && [[ -f "$STATE_CONTAINERS_FILE" ]]; then
        echo -e "${GREEN}Container setup already recorded as complete; skipping container rebuild.${NC}"
        echo "Marker file: $STATE_CONTAINERS_FILE"
        echo "Set MARGIE_FORCE_SETUP=1 to force a rebuild."
        exit 0
    fi
fi

# ── Docker build ────────────────────────────────────────────────────────────────
build_docker() {
    local dir="$1" tag="$2"
    echo -e "${YELLOW}[Docker] Building $tag ...${NC}"
    if [[ "$DRY_RUN" -eq 1 ]]; then
        print_cmd docker build --platform "$DOCKER_PLATFORM_RESOLVED" -t "$tag" "$CONTAINERS_DIR/$dir"
    elif docker image inspect "$tag" >/dev/null 2>&1; then
        echo -e "${YELLOW}  ! Docker image already exists: $tag (skipping)${NC}"
    else
        docker build --platform "$DOCKER_PLATFORM_RESOLVED" -t "$tag" "$CONTAINERS_DIR/$dir"
    fi
    echo -e "${GREEN}  ✓ $tag${NC}"
    ((BUILT++)) || true
}

# ── Apptainer/Singularity build ─────────────────────────────────────────────────
build_apptainer() {
    local dir="$1" tag="$2" cmd="$3"
    local def_file="$CONTAINERS_DIR/$dir/$dir.def"
    local sif_file="$CONTAINERS_DIR/$dir/${dir}.sif"
    if [[ "$DRY_RUN" -eq 0 && ! -f "$def_file" ]]; then
        echo -e "${RED}  ✗ $dir.def not found — skipping${NC}"
        return
    fi
    echo -e "${YELLOW}[$cmd] Building $dir.sif ...${NC}"
    if [[ "$DRY_RUN" -eq 1 ]]; then
        print_cmd "$cmd" build --force "$sif_file" "$def_file"
        echo -e "${GREEN}  ✓ $sif_file${NC}"
        ((BUILT++)) || true
    elif [[ -f "$sif_file" ]]; then
        echo -e "${YELLOW}  ! SIF already exists: $sif_file (skipping)${NC}"
        echo -e "${GREEN}  ✓ $sif_file${NC}"
        ((BUILT++)) || true
    else
        # cd into the container dir so %files relative paths (e.g. scripts/) resolve correctly
        if (cd "$CONTAINERS_DIR/$dir" && "$cmd" build --force "$sif_file" "$def_file"); then
            echo -e "${GREEN}  ✓ $sif_file${NC}"
            ((BUILT++)) || true
        else
            echo -e "${RED}  ✗ Build failed: $dir${NC}"
            return 1
        fi
    fi
}

pull_official_dbcan() {
    local sif_file="$CONTAINERS_DIR/dbcan/dbcan.sif"

    if [[ "$MODE" == "docker" ]]; then
        echo -e "${YELLOW}[Docker] Pulling official dbCAN image $OFFICIAL_DBCAN_IMAGE ...${NC}"
        if [[ "$DRY_RUN" -eq 1 ]]; then
            print_cmd docker pull --platform "$DOCKER_PLATFORM_RESOLVED" "$OFFICIAL_DBCAN_IMAGE"
            echo -e "${GREEN}  ✓ $OFFICIAL_DBCAN_IMAGE${NC}"
            ((BUILT++)) || true
            return 0
        fi

        if docker image inspect "$OFFICIAL_DBCAN_IMAGE" >/dev/null 2>&1; then
            echo -e "${YELLOW}  ! Docker image already exists: $OFFICIAL_DBCAN_IMAGE (skipping)${NC}"
            ((BUILT++)) || true
            return 0
        fi

        if docker pull --platform "$DOCKER_PLATFORM_RESOLVED" "$OFFICIAL_DBCAN_IMAGE"; then
            echo -e "${GREEN}  ✓ $OFFICIAL_DBCAN_IMAGE${NC}"
            ((BUILT++)) || true
            return 0
        fi

        echo -e "${RED}  ✗ Failed to pull $OFFICIAL_DBCAN_IMAGE${NC}"
        return 1
    fi

    if [[ "$MODE" == "apptainer" || "$MODE" == "singularity" ]]; then
        echo -e "${YELLOW}[$MODE] Pulling official dbCAN image into SIF ...${NC}"
        if [[ "$DRY_RUN" -eq 1 ]]; then
            print_cmd "$MODE" pull "$sif_file" "docker://$OFFICIAL_DBCAN_IMAGE"
            echo -e "${GREEN}  ✓ $sif_file${NC}"
            ((BUILT++)) || true
            return 0
        fi

        if "$MODE" pull "$sif_file" "docker://$OFFICIAL_DBCAN_IMAGE"; then
            echo -e "${GREEN}  ✓ $sif_file${NC}"
            ((BUILT++)) || true
            return 0
        fi

        if [[ -f "$sif_file" ]]; then
            echo -e "${YELLOW}  ! Pull failed but $sif_file already exists — using cached SIF.${NC}"
            ((BUILT++)) || true
            return 0
        fi
        echo -e "${RED}  ✗ Failed to pull official dbCAN image into $sif_file${NC}"
        return 1
    fi

    echo -e "${RED}  ✗ Unsupported runtime for official dbCAN image: $MODE${NC}"
    ((FAILED++)) || true
    return 1
}

pull_official_eggnog() {
    local sif_file="$CONTAINERS_DIR/eggnog/eggnog.sif"

    if [[ "$MODE" == "docker" ]]; then
        echo -e "${YELLOW}[Docker] Pulling official eggNOG image $OFFICIAL_EGGNOG_IMAGE ...${NC}"
        if [[ "$DRY_RUN" -eq 1 ]]; then
            print_cmd docker pull --platform "$OFFICIAL_EGGNOG_PLATFORM" "$OFFICIAL_EGGNOG_IMAGE"
            echo -e "${GREEN}  ✓ $OFFICIAL_EGGNOG_IMAGE${NC}"
            ((BUILT++)) || true
            return 0
        fi

        if docker image inspect "$OFFICIAL_EGGNOG_IMAGE" >/dev/null 2>&1; then
            echo -e "${YELLOW}  ! Docker image already exists: $OFFICIAL_EGGNOG_IMAGE (skipping)${NC}"
            ((BUILT++)) || true
            return 0
        fi

        if docker pull --platform "$OFFICIAL_EGGNOG_PLATFORM" "$OFFICIAL_EGGNOG_IMAGE"; then
            echo -e "${GREEN}  ✓ $OFFICIAL_EGGNOG_IMAGE${NC}"
            ((BUILT++)) || true
            return 0
        fi

        echo -e "${RED}  ✗ Failed to pull $OFFICIAL_EGGNOG_IMAGE${NC}"
        return 1
    fi

    if [[ "$MODE" == "apptainer" || "$MODE" == "singularity" ]]; then
        echo -e "${YELLOW}[$MODE] Pulling official eggNOG image into SIF ...${NC}"
        if [[ "$DRY_RUN" -eq 1 ]]; then
            print_cmd "$MODE" pull "$sif_file" "docker://$OFFICIAL_EGGNOG_IMAGE"
            echo -e "${GREEN}  ✓ $sif_file${NC}"
            ((BUILT++)) || true
            return 0
        fi

        if "$MODE" pull "$sif_file" "docker://$OFFICIAL_EGGNOG_IMAGE"; then
            echo -e "${GREEN}  ✓ $sif_file${NC}"
            ((BUILT++)) || true
            return 0
        fi

        if [[ -f "$sif_file" ]]; then
            echo -e "${YELLOW}  ! Pull failed but $sif_file already exists — using cached SIF.${NC}"
            ((BUILT++)) || true
            return 0
        fi
        echo -e "${RED}  ✗ Failed to pull official eggNOG image into $sif_file${NC}"
        return 1
    fi

    echo -e "${RED}  ✗ Unsupported runtime for official eggNOG image: $MODE${NC}"
    ((FAILED++)) || true
    return 1
}

have_local_operon_artifact() {
    local sif_file="$CONTAINERS_DIR/operon/operon.sif"

    if [[ "$MODE" == "docker" ]]; then
        docker image inspect "operon-annotation:1.0" >/dev/null 2>&1
        return $?
    fi

    if [[ "$MODE" == "apptainer" || "$MODE" == "singularity" ]]; then
        [[ -f "$sif_file" ]]
        return $?
    fi

    return 1
}

pull_official_operon() {
    local sif_file="$CONTAINERS_DIR/operon/operon.sif"

    if [[ "$MODE" == "docker" ]]; then
        echo -e "${YELLOW}[Docker] Pulling official operon image $OFFICIAL_OPERON_IMAGE ...${NC}"
        if [[ "$DRY_RUN" -eq 1 ]]; then
            print_cmd docker pull --platform "$DOCKER_PLATFORM_RESOLVED" "$OFFICIAL_OPERON_IMAGE"
            echo -e "${GREEN}  ✓ $OFFICIAL_OPERON_IMAGE${NC}"
            ((BUILT++)) || true
            return 0
        fi

        if docker image inspect "$OFFICIAL_OPERON_IMAGE" >/dev/null 2>&1; then
            echo -e "${YELLOW}  ! Docker image already exists: $OFFICIAL_OPERON_IMAGE (skipping)${NC}"
            ((BUILT++)) || true
            return 0
        fi

        if docker pull --platform "$DOCKER_PLATFORM_RESOLVED" "$OFFICIAL_OPERON_IMAGE"; then
            echo -e "${GREEN}  ✓ $OFFICIAL_OPERON_IMAGE${NC}"
            ((BUILT++)) || true
            return 0
        fi

        if have_local_operon_artifact; then
            echo -e "${YELLOW}  ! Failed to pull $OFFICIAL_OPERON_IMAGE; using locally built operon container instead.${NC}"
            return 0
        fi

        echo -e "${RED}  ✗ Failed to pull $OFFICIAL_OPERON_IMAGE${NC}"
        return 1
    fi

    if [[ "$MODE" == "apptainer" || "$MODE" == "singularity" ]]; then
        echo -e "${YELLOW}[$MODE] Pulling official operon image into SIF ...${NC}"
        if [[ "$DRY_RUN" -eq 1 ]]; then
            print_cmd "$MODE" pull "$sif_file" "docker://$OFFICIAL_OPERON_IMAGE"
            echo -e "${GREEN}  ✓ $sif_file${NC}"
            ((BUILT++)) || true
            return 0
        fi

        if "$MODE" pull "$sif_file" "docker://$OFFICIAL_OPERON_IMAGE"; then
            echo -e "${GREEN}  ✓ $sif_file${NC}"
            ((BUILT++)) || true
            return 0
        fi

        if have_local_operon_artifact; then
            echo -e "${YELLOW}  ! Failed to pull official operon image; using locally built $sif_file instead.${NC}"
            return 0
        fi

        echo -e "${RED}  ✗ Failed to pull official operon image into $sif_file${NC}"
        return 1
    fi

    echo -e "${RED}  ✗ Unsupported runtime for official operon image: $MODE${NC}"
    ((FAILED++)) || true
    return 1
}

pull_prebuilt_local_tag() {
    local local_tag="$1"
    local remote_image="$2"

    if [[ "$MODE" != "docker" ]]; then
        return 1
    fi

    echo -e "${YELLOW}[Docker] Attempting prebuilt image pull for $local_tag from $remote_image ...${NC}"
    if [[ "$DRY_RUN" -eq 1 ]]; then
        print_cmd docker pull --platform "$DOCKER_PLATFORM_RESOLVED" "$remote_image"
        print_cmd docker tag "$remote_image" "$local_tag"
        echo -e "${GREEN}  ✓ $local_tag (from prebuilt image)${NC}"
        ((BUILT++)) || true
        return 0
    fi

    if docker pull --platform "$DOCKER_PLATFORM_RESOLVED" "$remote_image"; then
        docker tag "$remote_image" "$local_tag"
        echo -e "${GREEN}  ✓ $local_tag (from prebuilt image)${NC}"
        ((BUILT++)) || true
        return 0
    fi

    echo -e "${YELLOW}  ! Prebuilt image unavailable for $local_tag; falling back to local build.${NC}"
    return 1
}

# ── Main build loop ──────────────────────────────────────────────────────────────
for entry in "${CONTAINERS[@]}"; do
    dir="${entry%%:*}"
    tag="${entry#*:}"

    if ! tool_selected "$dir"; then
        continue
    fi

    if [[ "$dir" == "prodigal" ]]; then
        if pull_prebuilt_local_tag "$tag" "$PREBUILT_PRODIGAL_IMAGE"; then
            continue
        fi
    fi

    if [[ "$dir" == "rasttk" ]]; then
        if pull_prebuilt_local_tag "$tag" "$PREBUILT_RASTTK_IMAGE"; then
            continue
        fi
    fi

    if [[ "$MODE" == "docker" ]]; then
        if ! build_docker "$dir" "$tag"; then
            if [[ "$dir" == "operon" ]]; then
                echo -e "${YELLOW}  ⚠ Optional container failed: $tag (continuing)${NC}"
            else
                echo -e "${RED}  ✗ Failed: $tag${NC}"
                ((FAILED++)) || true
            fi
        fi
    else
        if ! build_apptainer "$dir" "$tag" "$MODE"; then
            if [[ "$dir" == "operon" ]]; then
                echo -e "${YELLOW}  ⚠ Optional container failed: $dir (continuing)${NC}"
            else
                echo -e "${RED}  ✗ Failed: $dir${NC}"
                ((FAILED++)) || true
            fi
        fi
    fi
done

if tool_selected "dbcan"; then
    pull_official_dbcan || { echo -e "${RED}  ✗ Failed: official dbCAN image${NC}"; ((FAILED++)) || true; }
fi
if tool_selected "eggnog"; then
    pull_official_eggnog || { echo -e "${RED}  ✗ Failed: official eggNOG image${NC}"; ((FAILED++)) || true; }
fi
if tool_selected "operon"; then
    pull_official_operon || { echo -e "${RED}  ✗ Failed: official operon image${NC}"; ((FAILED++)) || true; }
fi

if [[ "$BUILD_GUI" -eq 1 ]]; then
    echo ""
    echo -e "${YELLOW}Building Margie GUI container...${NC}"
    gui_dir="$CONTAINERS_DIR/gui"
    gui_tag="margie-gui:1.0"
    if [[ "$MODE" == "docker" ]]; then
        if [[ "$DRY_RUN" -eq 1 ]]; then
            print_cmd docker build --platform "$DOCKER_PLATFORM_RESOLVED" -t "$gui_tag" "$gui_dir"
            echo -e "${GREEN}  ✓ $gui_tag${NC}"
            ((BUILT++)) || true
        elif docker image inspect "$gui_tag" >/dev/null 2>&1; then
            echo -e "${YELLOW}  ! GUI Docker image already exists: $gui_tag (skipping)${NC}"
            ((BUILT++)) || true
        elif docker build --platform "$DOCKER_PLATFORM_RESOLVED" -t "$gui_tag" "$gui_dir"; then
            echo -e "${GREEN}  ✓ $gui_tag${NC}"
            ((BUILT++)) || true
        else
            echo -e "${RED}  ✗ Failed to build GUI Docker image${NC}"
            ((FAILED++)) || true
        fi
    elif [[ "$MODE" == "apptainer" || "$MODE" == "singularity" ]]; then
        sif_gui="$gui_dir/gui.sif"
        def_gui="$gui_dir/gui.def"
        if [[ ! -f "$def_gui" ]]; then
            echo -e "${RED}  ✗ GUI definition file not found: $def_gui${NC}"
            ((FAILED++)) || true
        elif [[ "$DRY_RUN" -eq 1 ]]; then
            print_cmd "$MODE" build --force "$sif_gui" "$def_gui"
            echo -e "${GREEN}  ✓ $sif_gui${NC}"
            ((BUILT++)) || true
        elif [[ -f "$sif_gui" ]]; then
            echo -e "${YELLOW}  ! GUI SIF already exists: $sif_gui (skipping)${NC}"
            ((BUILT++)) || true
        elif (cd "$gui_dir" && "$MODE" build --force "$sif_gui" "$def_gui"); then
            echo -e "${GREEN}  ✓ $sif_gui${NC}"
            ((BUILT++)) || true
        else
            echo -e "${RED}  ✗ Failed to build GUI SIF${NC}"
            ((FAILED++)) || true
        fi
    else
        echo -e "${RED}  ✗ Unsupported runtime for GUI container: $MODE${NC}"
        ((FAILED++)) || true
    fi
fi

echo ""
echo -e "${GREEN}=== Build complete ===${NC}"
echo "  Built:  $BUILT"
echo "  Failed: $FAILED"

if [[ "$FAILED" -gt 0 ]]; then
    echo -e "${RED}Some containers failed. Check output above.${NC}"
    exit 1
fi

if [[ "$DRY_RUN" -eq 0 && "$TOOL_SELECTION" == "all" ]]; then
    mark_containers_ready
fi

if [[ "$SETUP_DB" -eq 1 ]]; then
    echo ""
    echo -e "${GREEN}=== Downloading databases ===${NC}"
    if [[ "$DB_SELECTION" == "all" ]]; then
        if [[ "$DRY_RUN" -eq 1 ]]; then
            "$SCRIPT_DIR/setup_databases.sh" --runtime "$MODE" --platform "$DOCKER_PLATFORM_RESOLVED" --dry-run
        else
            "$SCRIPT_DIR/setup_databases.sh" --runtime "$MODE" --platform "$DOCKER_PLATFORM_RESOLVED"
        fi
    else
        if [[ "$DRY_RUN" -eq 1 ]]; then
            "$SCRIPT_DIR/setup_databases.sh" --runtime "$MODE" --platform "$DOCKER_PLATFORM_RESOLVED" --db "$DB_SELECTION" --dry-run
        else
            "$SCRIPT_DIR/setup_databases.sh" --runtime "$MODE" --platform "$DOCKER_PLATFORM_RESOLVED" --db "$DB_SELECTION"
        fi
    fi
fi
