#!/usr/bin/env bash
# Multi-genome synteny workflow using MCScanX + DIAMOND in container.

set -euo pipefail

RASTTK_DIR="/container/output/rasttk"
CONSOLIDATED_DIR="/container/output/consolidation"
OUTDIR="/container/output/synteny"
THREADS="8"
GENOME_LIST=""
PREFIX="mcscanx_all"
ALIGNER="diamond"
MAX_HITS_PER_QUERY="30"
INTERGENOME_ONLY="1"
MIN_IDENTITY="25"
MIN_ALN_LEN="40"
GENERATE_GENE_REPORT="1"
RASTTK_REPORT_ROOT="/container/output/rasttk"
OPERON_REPORT_ROOT="/container/output/operon"

usage() {
  cat <<EOF
Usage: $0 [options]

Containerized synteny workflow for existing RASTtk outputs.

Options:
  --rasttk-dir DIR         RASTtk genome output root (default: /container/output/rasttk)
  --consolidated-dir DIR   Consolidated annotation root (default: /container/output/consolidation)
  --outdir DIR             Output root (default: /container/output/synteny)
  --genomes CSV            Comma-separated genome names to include (default: auto-discover)
  --threads N              Threads for DIAMOND (default: 8)
  --prefix NAME            MCScanX prefix basename (default: mcscanx_all)
  --aligner NAME           Sequence aligner: diamond|blastp|auto (default: diamond)
  --max-hits-per-query N   Keep top N hits per query for MCScanX (default: 30)
  --intergenome-only 0|1   Keep only cross-genome hits before MCScanX (default: 1)
  --min-identity N         Minimum percent identity before MCScanX (default: 25)
  --min-aln-len N          Minimum aa alignment length before MCScanX (default: 40)
  --generate-gene-report 0|1  Auto-generate gene JSON/TXT report (default: 1)
  --operon-root DIR        Operon root for gene report (default: /container/output/operon)
  -h, --help               Show this help

Expected per-genome files:
  <rasttk-dir>/<genome>/gene_calls/<genome>.faa
  <rasttk-dir>/<genome>/gene_calls/<genome>.gff

Outputs:
  <outdir>/native/<prefix>.gff
  <outdir>/native/<prefix>.blast
  <outdir>/native/<prefix>.collinearity
  <outdir>/synteny_gene_pairs.tsv
  <outdir>/processed/synteny_gene_pairs_detailed.tsv
  <outdir>/processed/synteny_annotation_comparison.tsv (if consolidated files are available)
  <outdir>/processed/synteny_block_sizes.png           (diagram)
  <outdir>/processed/human_report/                    (gene-centric reports, when enabled)
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --rasttk-dir) RASTTK_DIR="$2"; shift 2 ;;
    --consolidated-dir) CONSOLIDATED_DIR="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --genomes) GENOME_LIST="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    --aligner) ALIGNER="$2"; shift 2 ;;
    --max-hits-per-query) MAX_HITS_PER_QUERY="$2"; shift 2 ;;
    --intergenome-only) INTERGENOME_ONLY="$2"; shift 2 ;;
    --min-identity) MIN_IDENTITY="$2"; shift 2 ;;
    --min-aln-len) MIN_ALN_LEN="$2"; shift 2 ;;
    --generate-gene-report) GENERATE_GENE_REPORT="$2"; shift 2 ;;
    --operon-root) OPERON_REPORT_ROOT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

command -v diamond >/dev/null 2>&1 || true
command -v makeblastdb >/dev/null 2>&1 || true
command -v blastp >/dev/null 2>&1 || true
command -v MCScanX >/dev/null 2>&1 || { echo "ERROR: MCScanX not found in PATH"; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found in PATH"; exit 1; }

if [[ "$ALIGNER" != "diamond" && "$ALIGNER" != "blastp" && "$ALIGNER" != "auto" ]]; then
  echo "ERROR: --aligner must be one of: diamond, blastp, auto"
  exit 1
fi

if [[ "$ALIGNER" == "auto" ]]; then
  if command -v diamond >/dev/null 2>&1; then
    ALIGNER="diamond"
  elif command -v blastp >/dev/null 2>&1 && command -v makeblastdb >/dev/null 2>&1; then
    ALIGNER="blastp"
  else
    echo "ERROR: Neither DIAMOND nor BLAST+ tools are available"
    exit 1
  fi
fi

if [[ "$ALIGNER" == "diamond" ]]; then
  command -v diamond >/dev/null 2>&1 || { echo "ERROR: diamond not found in PATH"; exit 1; }
else
  command -v makeblastdb >/dev/null 2>&1 || { echo "ERROR: makeblastdb not found in PATH"; exit 1; }
  command -v blastp >/dev/null 2>&1 || { echo "ERROR: blastp not found in PATH"; exit 1; }
fi

mkdir -p "$OUTDIR"
NATIVE_DIR="$OUTDIR/native"
PROCESSED_DIR="$OUTDIR/processed"
REPORT_DIR="$PROCESSED_DIR/human_report"
mkdir -p "$NATIVE_DIR" "$PROCESSED_DIR" "$REPORT_DIR"

COMBINED_FAA="$NATIVE_DIR/${PREFIX}.faa"
MC_GFF="$NATIVE_DIR/${PREFIX}.gff"
RAW_BLAST="$NATIVE_DIR/${PREFIX}.blast.raw.tsv"
MC_BLAST="$NATIVE_DIR/${PREFIX}.blast"
PAIRS_TSV="$OUTDIR/synteny_gene_pairs.tsv"
PAIRS_DETAILED_TSV="$PROCESSED_DIR/synteny_gene_pairs_detailed.tsv"
GENE_META_TSV="$NATIVE_DIR/${PREFIX}.gene_loci.tsv"
ANNOTATION_JOIN_TSV="$PROCESSED_DIR/synteny_annotation_comparison.tsv"
LOG_FILE="$NATIVE_DIR/${PREFIX}.log"

: > "$LOG_FILE"
: > "$COMBINED_FAA"
: > "$MC_GFF"
echo -e "gene_id\tgenome\tcontig\tchr\tstart\tend\tstrand\tmidpoint" > "$GENE_META_TSV"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

log "Starting containerized MCScanX synteny workflow"
log "RASTtk dir: $RASTTK_DIR"
log "Consolidated dir: $CONSOLIDATED_DIR"
log "Output dir: $OUTDIR"
log "Aligner: $ALIGNER"
log "Core synteny = position/order/neighborhood from MCScanX anchors and DIAMOND/BLAST homology"
log "No hard requirements on same strand, operon membership, or functional linkage for block detection"
log "Hit filtering: topN=$MAX_HITS_PER_QUERY intergenome_only=$INTERGENOME_ONLY min_identity=$MIN_IDENTITY min_aln_len=$MIN_ALN_LEN"

if [[ ! -d "$RASTTK_DIR" ]]; then
  log "ERROR: RASTtk directory not found: $RASTTK_DIR"
  exit 1
fi

# Resolve genome list.
declare -a GENOMES
if [[ -n "$GENOME_LIST" ]]; then
  IFS=',' read -r -a GENOMES <<< "$GENOME_LIST"
else
  while IFS= read -r d; do
    GENOMES+=("$(basename "$d")")
  done < <(find "$RASTTK_DIR" -mindepth 1 -maxdepth 1 -type d | sort)
fi

if [[ ${#GENOMES[@]} -lt 2 ]]; then
  log "ERROR: Need at least 2 genomes for synteny analysis"
  exit 1
fi

log "Genomes selected (${#GENOMES[@]}): ${GENOMES[*]}"

VALID_GENOMES=0
GENOME_CHR_MAP=""
CHR_IDX=1
for genome in "${GENOMES[@]}"; do
  faa="$RASTTK_DIR/$genome/gene_calls/$genome.faa"
  gff="$RASTTK_DIR/$genome/gene_calls/$genome.gff"

  if [[ ! -f "$faa" ]]; then
    log "WARNING: Missing FAA, skipping genome: $genome"
    continue
  fi
  if [[ ! -f "$gff" ]]; then
    log "WARNING: Missing GFF, skipping genome: $genome"
    continue
  fi

  VALID_GENOMES=$((VALID_GENOMES + 1))

  cat "$faa" >> "$COMBINED_FAA"

  # MCScanX official GFF format is: chromosome_id, gene_id, midpoint
  awk -F'\t' -v OFS='\t' -v genome="$genome" -v meta="$GENE_META_TSV" '
    /^#/ { next }
    NF < 9 { next }
    $3 != "CDS" { next }
    {
      id = ""
      n = split($9, attrs, ";")
      for (i = 1; i <= n; i++) {
        if (attrs[i] ~ /^ID=/) {
          sub(/^ID=/, "", attrs[i])
          id = attrs[i]
          break
        }
      }
      if (id != "") {
        # Prodigal GFF IDs can be like "1_42" while FASTA headers are "contig_42".
        # Normalize IDs so BLAST qseqid/sseqid and MCScanX gene IDs match.
        if (id ~ /^[0-9]+_[0-9]+$/) {
          split(id, p, "_")
          id = $1 "_" p[2]
        }

        chr = genome "__" $1
        gsub(/[^A-Za-z0-9_.-]/, "_", chr)
        mid = int(($4 + $5) / 2)

        print chr, id, mid
        print id, genome, $1, chr, $4, $5, $7, mid >> meta
      }
    }
  ' "$gff" >> "$MC_GFF"
  
  CHR_IDX=$((CHR_IDX + 1))
done

if [[ $VALID_GENOMES -lt 2 ]]; then
  log "ERROR: Fewer than 2 genomes had both FAA and GFF"
  exit 1
fi

if [[ ! -s "$COMBINED_FAA" || ! -s "$MC_GFF" ]]; then
  log "ERROR: Required MCScanX inputs were not generated"
  exit 1
fi

# Sort GFF by chromosome and then by midpoint
sort -t$'\t' -k1,1 -k3,3n "$MC_GFF" > "${MC_GFF}.sorted"
mv "${MC_GFF}.sorted" "$MC_GFF"

if [[ "$ALIGNER" == "diamond" ]]; then
  log "Building DIAMOND database"
  diamond makedb --in "$COMBINED_FAA" --db "$NATIVE_DIR/${PREFIX}" >> "$LOG_FILE" 2>&1

  log "Running DIAMOND all-vs-all"
  diamond blastp \
    --query "$COMBINED_FAA" \
    --db "$NATIVE_DIR/${PREFIX}.dmnd" \
    --out "$RAW_BLAST" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    --threads "$THREADS" \
    --max-target-seqs 50 \
    --evalue 1e-5 >> "$LOG_FILE" 2>&1
else
  log "Building BLASTP database"
  makeblastdb -in "$COMBINED_FAA" -dbtype prot -out "$NATIVE_DIR/${PREFIX}" >> "$LOG_FILE" 2>&1

  log "Running BLASTP all-vs-all"
  blastp \
    -query "$COMBINED_FAA" \
    -db "$NATIVE_DIR/${PREFIX}" \
    -out "$RAW_BLAST" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -num_threads "$THREADS" \
    -max_target_seqs 50 \
    -evalue 1e-5 >> "$LOG_FILE" 2>&1
fi

if [[ ! -s "$RAW_BLAST" ]]; then
  log "ERROR: BLASTP output is empty"
  exit 1
fi

# Filter BLAST/DIAMOND hits for scalability before MCScanX.
log "Filtering alignments for scalable MCScanX input"
# Use streaming text tools to reduce peak memory on very large all-vs-all outputs.
FILTER_TMP1="$NATIVE_DIR/${PREFIX}.blast.prefilter.tsv"
FILTER_TMP2="$NATIVE_DIR/${PREFIX}.blast.sorted.tsv"

if [[ "$INTERGENOME_ONLY" == "1" ]]; then
  awk -F'\t' -v OFS='\t' -v minid="$MIN_IDENTITY" -v minlen="$MIN_ALN_LEN" '
    NR==FNR {
      if (FNR == 1) next
      gene2genome[$1] = $2
      next
    }
    {
      if ($1 == $2) next
      if (($3 + 0) < (minid + 0)) next
      if (($4 + 0) < (minlen + 0)) next
      gq = gene2genome[$1]
      gs = gene2genome[$2]
      if (gq == "" || gs == "") next
      if (gq == gs) next
      print $0
    }
  ' "$GENE_META_TSV" "$RAW_BLAST" > "$FILTER_TMP1"
else
  awk -F'\t' -v OFS='\t' -v minid="$MIN_IDENTITY" -v minlen="$MIN_ALN_LEN" '
    {
      if ($1 == $2) next
      if (($3 + 0) < (minid + 0)) next
      if (($4 + 0) < (minlen + 0)) next
      print $0
    }
  ' "$RAW_BLAST" > "$FILTER_TMP1"
fi

sort -T "$NATIVE_DIR" -S 1G -t $'\t' -k1,1 -k12,12gr -k11,11g "$FILTER_TMP1" > "$FILTER_TMP2"
awk -F'\t' -v OFS='\t' -v topn="$MAX_HITS_PER_QUERY" '
  {
    c[$1]++
    if (c[$1] <= topn) print $0
  }
' "$FILTER_TMP2" > "$MC_BLAST"

rm -f "$FILTER_TMP1" "$FILTER_TMP2"

if [[ ! -s "$MC_BLAST" ]]; then
  log "ERROR: Filtered BLAST file is empty; adjust filtering thresholds"
  exit 1
fi

log "Running MCScanX"
# Run with lenient parameters: min 2 genes per block (-s 2), try to find inter-genome blocks (-b 2)
MCScanX -s 2 -m 50 -w 10 "$NATIVE_DIR/$PREFIX" >> "$LOG_FILE" 2>&1 || {
  log "ERROR: MCScanX failed"
  exit 1
}

COLLINEARITY="$NATIVE_DIR/${PREFIX}.collinearity"
if [[ ! -f "$COLLINEARITY" ]]; then
  log "ERROR: Missing MCScanX collinearity output: $COLLINEARITY"
  exit 1
fi

log "Extracting gene pairs from collinearity output"
python3 /container/scripts/synteny_pairs_from_collinearity.py \
  --collinearity "$COLLINEARITY" \
  --output "$PAIRS_TSV" >> "$LOG_FILE" 2>&1

log "Annotating synteny pairs with genome loci"
python3 /container/scripts/annotate_synteny_pairs_with_loci.py \
  --pairs "$PAIRS_TSV" \
  --loci "$GENE_META_TSV" \
  --output "$PAIRS_DETAILED_TSV" >> "$LOG_FILE" 2>&1 || {
    log "WARNING: Pair loci annotation failed; base pair table still available"
  }

if [[ ! -s "$PAIRS_TSV" ]]; then
  log "WARNING: No syntenic pairs extracted"
fi

if [[ -d "$CONSOLIDATED_DIR" ]]; then
  log "Joining synteny pairs with consolidated annotations"
  python3 /container/work/processing/containers/synteny/scripts/join_synteny_with_annotations.py \
    --pairs "$PAIRS_TSV" \
    --consolidated-dir "$CONSOLIDATED_DIR" \
    --output "$ANNOTATION_JOIN_TSV" >> "$LOG_FILE" 2>&1 || {
      log "WARNING: Annotation join step failed; synteny pairs still available"
    }
else
  log "WARNING: Consolidated directory missing, skipping annotation join"
fi

if [[ "$GENERATE_GENE_REPORT" == "1" ]]; then
  if [[ -f "/container/workspace/processing/scripts/synteny/generate_synteny_gene_report.py" ]]; then
    log "Generating gene-centric JSON/TXT synteny report"
    python3 /container/workspace/processing/scripts/synteny/generate_synteny_gene_report.py \
      --pairs-detailed "$PAIRS_DETAILED_TSV" \
      --rasttk-root "$RASTTK_REPORT_ROOT" \
      --operon-root "$OPERON_REPORT_ROOT" \
      --output-json "$REPORT_DIR/gene_synteny_operon_report.json" \
      --output-txt "$REPORT_DIR/gene_synteny_operon_report.txt" \
      --max-links-per-gene 100 \
      --max-genes-in-txt 300 >> "$LOG_FILE" 2>&1 || {
        log "WARNING: gene-centric JSON/TXT report generation failed"
      }
  else
    log "WARNING: generate_synteny_gene_report.py not mounted; skipping auto gene report"
  fi
fi

log "Rendering synteny diagrams"
python3 /container/scripts/render_synteny_plots.py \
  --pairs "$PAIRS_TSV" \
  --output-dir "$PROCESSED_DIR" >> "$LOG_FILE" 2>&1 || {
    log "WARNING: Plot rendering failed"
  }

log "Workflow complete"
log "Main outputs:"
log "  $MC_GFF"
log "  $MC_BLAST"
log "  $COLLINEARITY"
log "  $PAIRS_TSV"
[[ -f "$PAIRS_DETAILED_TSV" ]] && log "  $PAIRS_DETAILED_TSV"
[[ -f "$ANNOTATION_JOIN_TSV" ]] && log "  $ANNOTATION_JOIN_TSV"
[[ -f "$PROCESSED_DIR/synteny_block_sizes.png" ]] && log "  $PROCESSED_DIR/synteny_block_sizes.png"
