#!/usr/bin/env bash
set -euo pipefail

WORK_DIR="${1:-/container/output/synteny_concept_rasttk}"
RAST_DIR="${2:-/container/output/rasttk}"
THREADS="${3:-8}"
ALIGNER_MODE="${4:-auto}"
GENOMES_CSV="${5:-}"

mkdir -p "$WORK_DIR"
: > "$WORK_DIR/all_proteins.faa"
: > "$WORK_DIR/all.simple.gff"
: > "$WORK_DIR/run.log"

echo "[info] WORK_DIR=$WORK_DIR" | tee -a "$WORK_DIR/run.log"
echo "[info] RAST_DIR=$RAST_DIR" | tee -a "$WORK_DIR/run.log"
echo "[info] aligner_mode=$ALIGNER_MODE" | tee -a "$WORK_DIR/run.log"

if [[ -n "$GENOMES_CSV" ]]; then
  IFS=',' read -r -a genomes <<< "$GENOMES_CSV"
else
  mapfile -t genomes < <(find "$RAST_DIR" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | grep -i "^Bacteroides_" | sort)
fi
echo "[info] selected_genomes=${#genomes[@]}" | tee -a "$WORK_DIR/run.log"

for base in "${genomes[@]}"; do
  gff="$RAST_DIR/$base/gene_calls/$base.gff"
  faa="$RAST_DIR/$base/gene_calls/$base.faa"

  if [[ -f "$gff" && -f "$faa" ]]; then
    echo "[info] processing=$base" >> "$WORK_DIR/run.log"

    # MCScanX official GFF columns: chr_id, gene_id, midpoint
    # Keep gene_id exactly equal to FAA header ID (ID=... in RASTtk GFF).
    awk -F"\t" -v OFS="\t" -v genome="$base" '
      /^#/ { next }
      $3 == "CDS" {
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
          chr = $1
          gsub(/[^A-Za-z0-9_.-]/, "_", chr)
          gsub(/[^A-Za-z0-9_.-]/, "_", genome)
          mid = int(($4 + $5) / 2)
          print genome "__" chr, id, mid
        }
      }
    ' "$gff" >> "$WORK_DIR/all.simple.gff"

    cat "$faa" >> "$WORK_DIR/all_proteins.faa"
  fi
done

cd "$WORK_DIR"
if [[ "$ALIGNER_MODE" == "diamond" || ( "$ALIGNER_MODE" == "auto" && -x "$(command -v diamond 2>/dev/null || true)" ) ]]; then
  command -v diamond >/dev/null 2>&1 || { echo "[error] aligner=diamond requested but diamond not found" | tee -a "$WORK_DIR/run.log"; exit 127; }
  echo "[info] aligner=diamond" | tee -a "$WORK_DIR/run.log"
  diamond makedb --in all_proteins.faa -d all.dmnd >> "$WORK_DIR/run.log" 2>&1
  diamond blastp \
    --db all.dmnd \
    --query all_proteins.faa \
    --out all.blast \
    --outfmt 6 \
    --evalue 1e-5 \
    --max-target-seqs 5 \
    --very-sensitive \
    --threads "$THREADS" >> "$WORK_DIR/run.log" 2>&1
else
  command -v makeblastdb >/dev/null 2>&1 || { echo "[error] makeblastdb not found" | tee -a "$WORK_DIR/run.log"; exit 127; }
  command -v blastp >/dev/null 2>&1 || { echo "[error] blastp not found" | tee -a "$WORK_DIR/run.log"; exit 127; }
  echo "[info] aligner=blastp" | tee -a "$WORK_DIR/run.log"
  makeblastdb -in all_proteins.faa -dbtype prot -out all >> "$WORK_DIR/run.log" 2>&1
  blastp \
    -query all_proteins.faa \
    -db all \
    -out all.blast \
    -outfmt 6 \
    -evalue 1e-5 \
    -max_target_seqs 5 \
    -num_threads "$THREADS" >> "$WORK_DIR/run.log" 2>&1
fi

cp all.simple.gff all.gff
MCScanX -s 2 -m 50 -w 10 all >> "$WORK_DIR/run.log" 2>&1 || true

echo "=== SUMMARY ===" | tee -a "$WORK_DIR/run.log"
wc -l all.gff all.blast all.collinearity 2>/dev/null | tee -a "$WORK_DIR/run.log" || true
grep -E "matches imported|pairwise comparisons|alignments generated|discarded" "$WORK_DIR/run.log" | tee -a "$WORK_DIR/run.log" || true
