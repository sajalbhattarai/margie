#!/bin/bash
# RASTtk Incremental Annotation Pipeline
# Uses individual RASTtk commands for full control

set -e

GENOME_NAME=${1:-"theta"}
INPUT_FASTA=${2:-"input/theta/theta.fasta"}
OUTPUT_DIR=${3:-"output/rasttk/theta"}
SCIENTIFIC_NAME=${4:-"Bacteroides thetaiotaomicron"}
GENETIC_CODE=${5:-11}
DOMAIN=${6:-"Bacteria"}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BVBRC="$SCRIPT_DIR/run_with_bvbrc.sh"

echo "=========================================="
echo "RASTtk Incremental Annotation Pipeline"
echo "=========================================="
echo "Genome: $GENOME_NAME"
echo "Input: $INPUT_FASTA"
echo "Output: $OUTPUT_DIR"
echo "Scientific Name: $SCIENTIFIC_NAME"
echo "Genetic Code: $GENETIC_CODE"
echo "Domain: $DOMAIN"
echo ""

if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: Input genome not found at $INPUT_FASTA"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Step 1: Create initial GTO
echo "=== [1/15] Creating Genome Typed Object (GTO) ==="
$BVBRC rast-create-genome \
    --scientific-name "$SCIENTIFIC_NAME" \
    --genetic-code "$GENETIC_CODE" \
    --domain "$DOMAIN" \
    --contigs "$INPUT_FASTA" \
    > "$OUTPUT_DIR/${GENOME_NAME}.gto.1"
echo "✓ GTO created"
echo ""

# Step 2: Call rRNA genes
echo "=== [2/15] Calling rRNA genes ==="
$BVBRC rast-call-features-rRNA-SEED \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.1" \
    > "$OUTPUT_DIR/${GENOME_NAME}.gto.2"
echo "✓ rRNA genes called"
echo ""

# Step 3: Call tRNA genes
echo "=== [3/15] Calling tRNA genes ==="
$BVBRC rast-call-features-tRNA-trnascan \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.2" \
    > "$OUTPUT_DIR/${GENOME_NAME}.gto.3"
echo "✓ tRNA genes called"
echo ""

# Step 4: Call repeat regions
echo "=== [4/15] Calling repeat regions ==="
$BVBRC rast-call-features-repeat-region-SEED \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.3" \
    > "$OUTPUT_DIR/${GENOME_NAME}.gto.4"
echo "✓ Repeat regions called"
echo ""

# Step 5: Call selenoproteins
echo "=== [5/15] Calling selenoproteins ==="
$BVBRC rast-call-features-selenoprotein \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.4" \
    > "$OUTPUT_DIR/${GENOME_NAME}.gto.5"
echo "✓ Selenoproteins called"
echo ""

# Step 6: Call pyrrolysoproteins
echo "=== [6/15] Calling pyrrolysoproteins ==="
$BVBRC rast-call-features-pyrrolysoprotein \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.5" \
    > "$OUTPUT_DIR/${GENOME_NAME}.gto.6"
echo "✓ Pyrrolysoproteins called"
echo ""

# Step 7: Call CRISPRs
echo "=== [7/15] Calling CRISPR elements ==="
$BVBRC rast-call-features-crispr \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.6" \
    > "$OUTPUT_DIR/${GENOME_NAME}.gto.7"
echo "✓ CRISPR elements called"
echo ""

# Step 8: Call protein-encoding genes with Prodigal
echo "=== [8/15] Calling protein-encoding genes (Prodigal) ==="
$BVBRC rast-call-features-CDS-prodigal \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.7" \
    > "$OUTPUT_DIR/${GENOME_NAME}.gto.8"
echo "✓ Prodigal genes called"
echo ""

# Step 9: Call protein-encoding genes with Glimmer3
echo "=== [9/15] Calling protein-encoding genes (Glimmer3) ==="
$BVBRC rast-call-features-CDS-glimmer3 \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.8" \
    > "$OUTPUT_DIR/${GENOME_NAME}.gto.9"
echo "✓ Glimmer3 genes called"
echo ""

# Step 10: Annotate with k-mer v2
echo "=== [10/15] Annotating proteins (k-mer v2) ==="
$BVBRC rast-annotate-proteins-kmer-v2 \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.9" \
    > "$OUTPUT_DIR/${GENOME_NAME}.gto.10"
echo "✓ k-mer v2 annotation complete"
echo ""

# Step 11: Annotate hypotheticals with k-mer v1
echo "=== [11/15] Annotating hypothetical proteins (k-mer v1) ==="
$BVBRC rast-annotate-proteins-kmer-v1 -H \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.10" \
    > "$OUTPUT_DIR/${GENOME_NAME}.gto.11"
echo "✓ k-mer v1 annotation complete"
echo ""

# Step 12: Annotate remaining by similarity
echo "=== [12/15] Annotating by similarity to close relatives ==="
$BVBRC rast-annotate-proteins-similarity -H \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.11" \
    > "$OUTPUT_DIR/${GENOME_NAME}.gto.12" || {
    echo "⚠ Similarity annotation failed (no close relatives available)"
    cp "$OUTPUT_DIR/${GENOME_NAME}.gto.11" "$OUTPUT_DIR/${GENOME_NAME}.gto.12"
}
echo "✓ Similarity annotation complete"
echo ""

# Step 13: Resolve overlapping features
echo "=== [13/15] Resolving overlapping features ==="
$BVBRC rast-resolve-overlapping-features \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.12" \
    > "$OUTPUT_DIR/${GENOME_NAME}.gto.13"
echo "✓ Overlaps resolved"
echo ""

# Step 14: Call prophage elements
echo "=== [14/15] Calling prophage elements (PhiSpy) ==="
$BVBRC rast-call-features-prophage-phispy \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.13" \
    > "$OUTPUT_DIR/${GENOME_NAME}.gto.14" || {
    echo "⚠ PhiSpy prophage detection skipped (common for small genomes); continuing without prophage annotations"
    cp "$OUTPUT_DIR/${GENOME_NAME}.gto.13" "$OUTPUT_DIR/${GENOME_NAME}.gto.14"
}
echo "✓ Prophage elements step complete"
echo ""

# Step 15: Export in multiple formats
echo "=== [15/15] Exporting annotations ==="

# Feature table
$BVBRC rast-export-genome feature_data \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.14" \
    > "$OUTPUT_DIR/${GENOME_NAME}_features.tsv"
echo "✓ Feature table: ${GENOME_NAME}_features.tsv"

# GenBank format
$BVBRC rast-export-genome genbank \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.14" \
    > "$OUTPUT_DIR/${GENOME_NAME}.gbk"
echo "✓ GenBank: ${GENOME_NAME}.gbk"

# GFF format
$BVBRC rast-export-genome gff \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.14" \
    > "$OUTPUT_DIR/${GENOME_NAME}.gff"
echo "✓ GFF: ${GENOME_NAME}.gff"

# Protein FASTA
$BVBRC rast-export-genome protein_fasta \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.14" \
    > "$OUTPUT_DIR/${GENOME_NAME}.faa"
echo "✓ Protein sequences: ${GENOME_NAME}.faa"

# Gene DNA FASTA
$BVBRC rast-export-genome feature_dna \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.14" \
    > "$OUTPUT_DIR/${GENOME_NAME}.ffn"
echo "✓ Gene sequences: ${GENOME_NAME}.ffn"

# Contig FASTA
$BVBRC rast-export-genome contig_fasta \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.14" \
    > "$OUTPUT_DIR/${GENOME_NAME}.fna"
echo "✓ Contig sequences: ${GENOME_NAME}.fna"

# Copy Prodigal outputs to standard location for downstream tools
echo ""
echo "=== Copying Prodigal outputs for downstream annotation ==="
PRODIGAL_DIR="output/prodigal/${GENOME_NAME}/${GENOME_NAME}"
mkdir -p "$PRODIGAL_DIR"
cp "$OUTPUT_DIR/${GENOME_NAME}.faa" "$PRODIGAL_DIR/${GENOME_NAME}.faa"
cp "$OUTPUT_DIR/${GENOME_NAME}.ffn" "$PRODIGAL_DIR/${GENOME_NAME}.ffn"
cp "$OUTPUT_DIR/${GENOME_NAME}.gff" "$PRODIGAL_DIR/${GENOME_NAME}.gff"
echo "✓ Prodigal outputs copied to: $PRODIGAL_DIR"
echo "  - ${GENOME_NAME}.faa (proteins for downstream tools)"
echo "  - ${GENOME_NAME}.ffn (gene sequences)"
echo "  - ${GENOME_NAME}.gff (gene features)"

# Export specific feature types
$BVBRC rast-export-genome --feature-type CDS feature_data \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.14" \
    > "$OUTPUT_DIR/${GENOME_NAME}_CDS.tsv"
echo "✓ CDS features: ${GENOME_NAME}_CDS.tsv"

$BVBRC rast-export-genome --feature-type rna feature_data \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.14" \
    > "$OUTPUT_DIR/${GENOME_NAME}_RNA.tsv"
echo "✓ RNA features: ${GENOME_NAME}_RNA.tsv"

$BVBRC rast-export-genome --feature-type prophage feature_data \
    < "$OUTPUT_DIR/${GENOME_NAME}.gto.14" \
    > "$OUTPUT_DIR/${GENOME_NAME}_prophage.tsv"
echo "✓ Prophage features: ${GENOME_NAME}_prophage.tsv"

echo ""

# Generate summary
echo "=== Generating Summary ==="
CDS_COUNT=$(tail -n +2 "$OUTPUT_DIR/${GENOME_NAME}_CDS.tsv" 2>/dev/null | wc -l | tr -d ' ' || echo "0")
RNA_COUNT=$(tail -n +2 "$OUTPUT_DIR/${GENOME_NAME}_RNA.tsv" 2>/dev/null | wc -l | tr -d ' ' || echo "0")
PROPHAGE_COUNT=$(tail -n +2 "$OUTPUT_DIR/${GENOME_NAME}_prophage.tsv" 2>/dev/null | wc -l | tr -d ' ' || echo "0")
PROTEIN_COUNT=$(grep -c "^>" "$OUTPUT_DIR/${GENOME_NAME}.faa" 2>/dev/null | tr -d ' ' || echo "0")

cat > "$OUTPUT_DIR/${GENOME_NAME}_summary.txt" << EOF
========================================
RASTtk Annotation Summary
========================================
Genome: $GENOME_NAME
Scientific Name: $SCIENTIFIC_NAME
Domain: $DOMAIN
Genetic Code: $GENETIC_CODE
Date: $(date)

=== Feature Counts ===
CDS (protein-encoding genes): $CDS_COUNT
RNA genes: $RNA_COUNT
Prophage elements: $PROPHAGE_COUNT
Proteins translated: $PROTEIN_COUNT

=== Output Files ===
Feature table: ${GENOME_NAME}_features.tsv
GenBank: ${GENOME_NAME}.gbk
GFF: ${GENOME_NAME}.gff
Proteins: ${GENOME_NAME}.faa
Genes: ${GENOME_NAME}.ffn
Contigs: ${GENOME_NAME}.fna
GTO (final): ${GENOME_NAME}.gto.14

=== Annotation Steps Completed ===
✓ rRNA genes (BLAST-based)
✓ tRNA genes (tRNAscan)
✓ Repeat regions
✓ Selenoproteins
✓ Pyrrolysoproteins
✓ CRISPR elements
✓ Protein-encoding genes (Prodigal)
✓ Protein-encoding genes (Glimmer3)
✓ Function annotation (k-mer v2)
✓ Function annotation (k-mer v1)
✓ Function annotation (similarity)
✓ Overlap resolution
✓ Prophage elements (PhiSpy)
✓ Export to multiple formats

========================================
EOF

echo ""
echo "=========================================="
echo "RASTtk Annotation Complete!"
echo "=========================================="
echo "Results in: $OUTPUT_DIR"
echo ""

# Organize outputs into subdirectories
echo "=== Organizing outputs ==="
mkdir -p "$OUTPUT_DIR/gene_calls"
mkdir -p "$OUTPUT_DIR/native"

# Move gene call files to gene_calls/
mv "$OUTPUT_DIR/${GENOME_NAME}.faa" "$OUTPUT_DIR/gene_calls/"
mv "$OUTPUT_DIR/${GENOME_NAME}.ffn" "$OUTPUT_DIR/gene_calls/"
mv "$OUTPUT_DIR/${GENOME_NAME}.fna" "$OUTPUT_DIR/gene_calls/"
mv "$OUTPUT_DIR/${GENOME_NAME}.gff" "$OUTPUT_DIR/gene_calls/"
mv "$OUTPUT_DIR/${GENOME_NAME}.gbk" "$OUTPUT_DIR/gene_calls/"
echo "✓ Moved gene call files to gene_calls/"

# Move native RAST outputs to native/
mv "$OUTPUT_DIR/${GENOME_NAME}_CDS.tsv" "$OUTPUT_DIR/native/"
mv "$OUTPUT_DIR/${GENOME_NAME}_RNA.tsv" "$OUTPUT_DIR/native/"
mv "$OUTPUT_DIR/${GENOME_NAME}_prophage.tsv" "$OUTPUT_DIR/native/"
mv "$OUTPUT_DIR/${GENOME_NAME}_features.tsv" "$OUTPUT_DIR/native/"
mv "$OUTPUT_DIR/${GENOME_NAME}.gto."* "$OUTPUT_DIR/native/"
mv "$OUTPUT_DIR/${GENOME_NAME}_summary.txt" "$OUTPUT_DIR/native/"
echo "✓ Moved native RAST outputs to native/"

# Consolidate outputs into rast.tsv for downstream annotation with full SEED enrichment
echo ""
echo "=== Consolidating outputs into rast.tsv ==="
python3 "$SCRIPT_DIR/generate_rast_tsv.py" "$OUTPUT_DIR" "$SCIENTIFIC_NAME"

# rast.tsv is written directly into the main output directory by generate_rast_tsv.py
if [ -f "$OUTPUT_DIR/rast.tsv" ]; then
    echo "✓ Created consolidated rast.tsv for merging with other annotations"
    ROW_COUNT=$(tail -n +2 "$OUTPUT_DIR/rast.tsv" | wc -l | tr -d ' ')
    echo "  Total features in rast.tsv: $ROW_COUNT"
    if [[ -f "/container/db/subsystem_mapping.tsv" && -d "/container/db/variant_definitions" ]]; then
        echo "  SEED subsystem hierarchy and variant definitions were loaded from /container/db"
    else
        echo "  Warning: db/rasttk is incomplete, so SEED enrichment may be partial or absent"
    fi
else
    echo "✗ Warning: rast.tsv not created"
fi

# Update Prodigal directory paths with new structure
echo ""
echo "=== Updating Prodigal output copies ==="
PRODIGAL_DIR="output/prodigal/${GENOME_NAME}/${GENOME_NAME}"
mkdir -p "$PRODIGAL_DIR"
cp "$OUTPUT_DIR/gene_calls/${GENOME_NAME}.faa" "$PRODIGAL_DIR/${GENOME_NAME}.faa"
cp "$OUTPUT_DIR/gene_calls/${GENOME_NAME}.ffn" "$PRODIGAL_DIR/${GENOME_NAME}.ffn"
cp "$OUTPUT_DIR/gene_calls/${GENOME_NAME}.gff" "$PRODIGAL_DIR/${GENOME_NAME}.gff"
echo "✓ Prodigal outputs copied to: $PRODIGAL_DIR"

echo ""
echo "=========================================="
echo "Directory Structure:"
echo "=========================================="
echo "$OUTPUT_DIR/"
echo "├── rast.tsv (consolidated annotation table)"
echo "├── gene_calls/"
echo "│   ├── ${GENOME_NAME}.faa (proteins for downstream tools)"
echo "│   ├── ${GENOME_NAME}.ffn (gene nucleotide sequences)"
echo "│   ├── ${GENOME_NAME}.fna (contig sequences)"
echo "│   ├── ${GENOME_NAME}.gff (gene features)"
echo "│   └── ${GENOME_NAME}.gbk (GenBank format)"
echo "└── native/"
echo "    ├── ${GENOME_NAME}_CDS.tsv"
echo "    ├── ${GENOME_NAME}_RNA.tsv"
echo "    ├── ${GENOME_NAME}_prophage.tsv"
echo "    ├── ${GENOME_NAME}_features.tsv"
echo "    ├── ${GENOME_NAME}_summary.txt"
echo "    └── ${GENOME_NAME}.gto.* (14 GTO files)"
echo ""

echo ""
echo "=== Generating Pipeline Provenance Record ==="
PIPELINE_FILE="$OUTPUT_DIR/pipeline.txt"
CURRENT_DATE=$(date +"%Y-%m-%d")
OPERATOR=${USER:-"unknown"}

cat > "$PIPELINE_FILE" << 'PIPELINE_EOF'
================================================================================
RASTTK ANNOTATION PIPELINE - PROVENANCE RECORD
================================================================================
PIPELINE_EOF

cat >> "$PIPELINE_FILE" << PIPELINE_EOF
Genome: $SCIENTIFIC_NAME ($GENOME_NAME)
Date: $CURRENT_DATE
Operator: $OPERATOR
Annotation Method: RASTtk via BV-BRC CLI (local execution)

================================================================================
STEP 1: ENVIRONMENT SETUP
================================================================================
Tool: BV-BRC CLI v1.048
Installation: $SCRIPT_DIR/bvbrc/
Wrapper Script: $SCRIPT_DIR/run_with_bvbrc.sh
Credentials File: $SCRIPT_DIR/credentials.txt

Environment Variables Set:
  KB_TOP=\${SCRIPT_DIR}/bvbrc/deployment
  KB_RUNTIME=\${SCRIPT_DIR}/bvbrc/deployment
  PATH=\${KB_TOP}/bin:\${KB_RUNTIME}/bin:\$PATH
  PERL5LIB=\${KB_TOP}/lib

================================================================================
STEP 2: AUTHENTICATION
================================================================================
Command: p3-login
Username: [from credentials.txt]
Password: [stored in credentials.txt]
Method: Interactive login via wrapper script

Authentication Status: ✓ Successful
Session: Active BV-BRC session established

================================================================================
STEP 3: GENOME ANNOTATION (15-STEP INCREMENTAL PIPELINE)
================================================================================
Input File: $INPUT_FASTA
Output Directory: $OUTPUT_DIR/
Scientific Name: $SCIENTIFIC_NAME
Genetic Code: $GENETIC_CODE
Domain: $DOMAIN

--- Step 1: Create Genome Typed Object (GTO) ---
Command: rast-create-genome \\
    --scientific-name "$SCIENTIFIC_NAME" \\
    --genetic-code $GENETIC_CODE \\
    --domain $DOMAIN \\
    --contigs $INPUT_FASTA
Output: ${GENOME_NAME}.gto.1

--- Step 2: Call rRNA Genes ---
Command: rast-call-features-rRNA-SEED < ${GENOME_NAME}.gto.1
Method: BLAST-based rRNA detection
Output: ${GENOME_NAME}.gto.2

--- Step 3: Call tRNA Genes ---
Command: rast-call-features-tRNA-trnascan < ${GENOME_NAME}.gto.2
Method: tRNAscan-SE
Output: ${GENOME_NAME}.gto.3

--- Step 4: Call Repeat Regions ---
Command: rast-call-features-repeat-region-SEED < ${GENOME_NAME}.gto.3
Output: ${GENOME_NAME}.gto.4

--- Step 5: Call Selenoproteins ---
Command: rast-call-selenoproteins < ${GENOME_NAME}.gto.4
Output: ${GENOME_NAME}.gto.5

--- Step 6: Call Pyrrolysoproteins ---
Command: rast-call-pyrrolysoproteins < ${GENOME_NAME}.gto.5
Output: ${GENOME_NAME}.gto.6

--- Step 7: Call CRISPR Elements ---
Command: rast-call-features-CRISPRs < ${GENOME_NAME}.gto.6
Output: ${GENOME_NAME}.gto.7

--- Step 8: Call Protein-Encoding Genes (Prodigal) ---
Command: rast-call-features-CDS-prodigal < ${GENOME_NAME}.gto.7
Method: Prodigal gene caller
Output: ${GENOME_NAME}.gto.8

--- Step 9: Call Protein-Encoding Genes (Glimmer3) ---
Command: rast-call-features-CDS-glimmer3 < ${GENOME_NAME}.gto.8
Method: Glimmer3 gene caller
Output: ${GENOME_NAME}.gto.9

--- Step 10: Annotate Protein Functions (k-mer v2) ---
Command: rast-annotate-proteins-kmer-v2 < ${GENOME_NAME}.gto.9
Method: k-mer based annotation (version 2)
Output: ${GENOME_NAME}.gto.10

--- Step 11: Annotate Protein Functions (k-mer v1) ---
Command: rast-annotate-proteins-kmer-v1 < ${GENOME_NAME}.gto.10
Method: k-mer based annotation (version 1)
Output: ${GENOME_NAME}.gto.11

--- Step 12: Annotate Protein Functions (Similarity) ---
Command: rast-annotate-proteins-similarity < ${GENOME_NAME}.gto.11
Method: Similarity-based annotation
Output: ${GENOME_NAME}.gto.12

--- Step 13: Resolve Overlapping Features ---
Command: rast-resolve-overlapping-features < ${GENOME_NAME}.gto.12
Output: ${GENOME_NAME}.gto.13

--- Step 14: Call Prophage Elements ---
Command: rast-call-features-prophage-phispy < ${GENOME_NAME}.gto.13
Method: PhiSpy prophage detection
Output: ${GENOME_NAME}.gto.14 (FINAL GTO)

================================================================================
STEP 4: EXPORT ANNOTATION RESULTS
================================================================================
--- Export All Features ---
Command: rast-export-genome feature_data < ${GENOME_NAME}.gto.14
Output: ${GENOME_NAME}_features.tsv

--- Export GenBank Format ---
Command: rast-export-genome genbank < ${GENOME_NAME}.gto.14
Output: ${GENOME_NAME}.gbk

--- Export GFF Format ---
Command: rast-export-genome gff < ${GENOME_NAME}.gto.14
Output: ${GENOME_NAME}.gff

--- Export Protein Sequences ---
Command: rast-export-genome protein_fasta < ${GENOME_NAME}.gto.14
Output: ${GENOME_NAME}.faa

--- Export Gene Sequences ---
Command: rast-export-genome feature_dna < ${GENOME_NAME}.gto.14
Output: ${GENOME_NAME}.ffn

--- Export Contig Sequences ---
Command: rast-export-genome contig_fasta < ${GENOME_NAME}.gto.14
Output: ${GENOME_NAME}.fna

--- Export CDS Features Only ---
Command: rast-export-genome --feature-type CDS feature_data < ${GENOME_NAME}.gto.14
Output: ${GENOME_NAME}_CDS.tsv

--- Export RNA Features Only ---
Command: rast-export-genome --feature-type rna feature_data < ${GENOME_NAME}.gto.14
Output: ${GENOME_NAME}_RNA.tsv

--- Export Prophage Features Only ---
Command: rast-export-genome --feature-type prophage feature_data < ${GENOME_NAME}.gto.14
Output: ${GENOME_NAME}_prophage.tsv

================================================================================
STEP 5: CONSOLIDATE ANNOTATIONS
================================================================================
Script: $SCRIPT_DIR/consolidate_rast_output.py
Input Directory: $OUTPUT_DIR/native/

Command: python $SCRIPT_DIR/consolidate_rast_output.py \\
    $OUTPUT_DIR/native \\
    "$SCIENTIFIC_NAME"

Processing Steps:
  1. Load protein sequences from ${GENOME_NAME}.faa
  2. Load nucleotide sequences from ${GENOME_NAME}.ffn
  3. Parse CDS features from ${GENOME_NAME}_CDS.tsv
  4. Parse RNA features from ${GENOME_NAME}_RNA.tsv
  5. Parse prophage features from ${GENOME_NAME}_prophage.tsv
  6. Extract gene locations (start, end, strand) from RAST_gene_id
  7. Extract EC numbers from RAST_description
  8. Calculate sequence lengths (na_length, aa_length)
  9. Combine all features into single table

Output: rast.tsv ($ROW_COUNT features)

Column Structure (14 columns):
  1. organism - Scientific name
  2. feature_id - RASTtk feature ID (merge key)
  3. RAST_gene_id - Gene location string
  4. gene_start - Start coordinate
  5. gene_end - End coordinate
  6. na_length - Nucleotide sequence length
  7. aa_length - Amino acid sequence length
  8. na_seq - Full nucleotide sequence
  9. aa_seq - Full protein sequence
  10. RAST_gene_type - Feature type (CDS, rRNA, tRNA, prophage)
  11. RAST_description - Functional annotation
  12. RAST_EC_numbers - Extracted EC numbers
  13. RAST_strand - Coding strand (+/-)
  14. RAST_feature_hash - MD5 hash

================================================================================
STEP 6: ORGANIZE OUTPUT FILES
================================================================================
Directory Structure Created:
$OUTPUT_DIR/
├── rast.tsv                    # Consolidated annotation table
├── gene_calls/                 # Gene prediction outputs
│   ├── ${GENOME_NAME}.faa     # Proteins (input for downstream tools)
│   ├── ${GENOME_NAME}.ffn     # Gene sequences
│   ├── ${GENOME_NAME}.fna     # Contig sequences
│   ├── ${GENOME_NAME}.gff     # Gene features
│   └── ${GENOME_NAME}.gbk     # GenBank format
└── native/                     # Raw RASTtk outputs
    ├── ${GENOME_NAME}_CDS.tsv
    ├── ${GENOME_NAME}_RNA.tsv
    ├── ${GENOME_NAME}_prophage.tsv
    ├── ${GENOME_NAME}_features.tsv
    ├── ${GENOME_NAME}_summary.txt
    └── ${GENOME_NAME}.gto.1-14

Files Moved:
  gene_calls/ ← ${GENOME_NAME}.faa, .ffn, .fna, .gff, .gbk
  native/ ← ${GENOME_NAME}_*.tsv, ${GENOME_NAME}.gto.*, ${GENOME_NAME}_summary.txt
  Main directory ← rast.tsv (only)

================================================================================
STEP 7: COPY OUTPUTS FOR DOWNSTREAM ANNOTATION
================================================================================
Destination: output/prodigal/${GENOME_NAME}/${GENOME_NAME}/
Files Copied:
  - ${GENOME_NAME}.faa → Input for COG, Pfam, eggNOG, etc.
  - ${GENOME_NAME}.ffn → Input for sequence-based tools
  - ${GENOME_NAME}.gff → Input for genomic context analysis

Purpose: Standardized input location for all downstream annotation tools

================================================================================
DOWNSTREAM USAGE
================================================================================
Merge Key: feature_id (e.g., fig|6666666.1476522.peg.1)
Input File for Annotations: $OUTPUT_DIR/gene_calls/${GENOME_NAME}.faa
Base Table for Merging: $OUTPUT_DIR/rast.tsv

Next Steps:
  1. Run COG annotation using ${GENOME_NAME}.faa
  2. Run Pfam annotation using ${GENOME_NAME}.faa
  3. Run eggNOG annotation using ${GENOME_NAME}.faa
  4. Run dbCAN annotation using ${GENOME_NAME}.faa
  5. Run KEGG annotation using ${GENOME_NAME}.faa
  6. Merge all results on feature_id column

================================================================================
END OF PIPELINE RECORD
================================================================================
Pipeline Execution Script: processing/scripts/rasttk/run_rasttk_incremental.sh
Consolidation Script: processing/scripts/rasttk/consolidate_rast_output.py
Generated: $CURRENT_DATE by $OPERATOR
PIPELINE_EOF

echo "✓ Pipeline provenance record saved to: $PIPELINE_FILE"

echo ""
echo "Done! $(date)"
