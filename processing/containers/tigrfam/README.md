# TIGRfam Annotation Container

## Overview

TIGRfam protein family annotation using HMMER hmmscan. Identifies TIGRfam families in protein sequences with confidence scoring.

## Container Details

- **Image:** tigrfam-annotation:1.0
- **Size:** 164 MB
- **Platform:** linux/amd64
- **Tools:** HMMER 3.4, Python 3
- **Database:** External mount required

## Architecture

```
Lean container (164 MB)
├── HMMER 3.4 (hmmscan)
├── Python 3 scripts
└── Mount points:
    ├── /container/input  (protein FASTA)
    ├── /container/output (results)
    └── /container/db     (TIGRfam database)
```

## Output Format

### tigrfam.tsv (26 columns)

One row per domain hit with confidence scores:

1. organism - Genome name
2. feature_id - Protein ID
3. protein_length - Full protein length (aa)
4. domain_number - Domain position in protein
5. TIGRFAM_family - Family ID (e.g., TIGR00001)
6. TIGRFAM_description - Family function
7. TIGRFAM_isotype - Type (equivalog, subfamily, etc.)
8. TIGRFAM_tlen - HMM model length
9-14. TIGRFAM_*_from/to - Domain boundaries (HMM, alignment, envelope)
15. TIGRFAM_full_evalue - Full sequence E-value
16. TIGRFAM_full_score - Full sequence bit score
17. TIGRFAM_full_bias - Full sequence bias
18. TIGRFAM_domain_cevalue - Conditional E-value
19. TIGRFAM_domain_ievalue - Independent E-value
20. TIGRFAM_domain_score - Domain bit score
21. TIGRFAM_domain_bias - Domain bias
22. TIGRFAM_acc - Accession number
23. TIGRFAM_evalue_score - E-value component (0-1)
24. TIGRFAM_coverage_score - Coverage component (0-1)
25. TIGRFAM_alignment_score - Bit score component (0-1)
26. TIGRFAM_confidence - Overall confidence (0-1)

### Confidence Scoring

```
confidence = 0.6 × evalue_score + 0.3 × coverage_score + 0.1 × alignment_score

Components:
- evalue_score: -log10(E-value) / 50, capped at 1.0
- coverage_score: (hmm_to - hmm_from + 1) / tlen
- alignment_score: bit_score / 100, capped at 1.0
```

## Database Setup

### Required Files

Place in `db/tigrfam/`:

1. **TIGRFAMs_15.0_HMM.LIB** - Main HMM database (~45 MB)
2. **tigrfam_definitions.txt** - Family descriptions (tab-delimited)
3. **TIGRFAMs_15.0_HMM.LIB.h3[fimp]** - HMMER indices (auto-generated)

### Download Instructions

```bash
# Navigate to database directory
cd db/tigrfam

# Download TIGRfam database
wget ftp://ftp.jcvi.org/data/TIGRFAMs/TIGRFAMs_15.0_HMM.LIB.gz
gunzip TIGRFAMs_15.0_HMM.LIB.gz

# Press database for HMMER
hmmpress TIGRFAMs_15.0_HMM.LIB

# Extract definitions
grep -E "^NAME|^DESC|^ISOTYPE" TIGRFAMs_15.0_HMM.LIB | \
  paste - - - | \
  sed 's/NAME  //; s/  DESC  /\t/; s/  ISOTYPE  /\t/' \
  > tigrfam_definitions.txt
```

See [db/tigrfam/README.md](../../../db/tigrfam/README.md) for details.

## Usage

### Via Pipeline

```bash
# Annotate with TIGRfam
./annotate-genome --select test --tools tigrfam -t 4

# Combine with other tools
./annotate-genome --select test --tools rasttk,kegg,pfam,tigrfam -t 8
```

### Direct Container Usage

```bash
# Step 1: Run HMMER scan
docker run --rm --platform linux/amd64 \
  -v "$(pwd)/proteins.faa:/container/input/proteins.faa:ro" \
  -v "$(pwd)/output:/container/output:rw" \
  -v "$(pwd)/db/tigrfam:/container/db:ro" \
  tigrfam-annotation:1.0 \
  /container/scripts/run_hmmscan.sh \
  "genome_name" \
  "/container/input/proteins.faa" \
  "/container/db/TIGRFAMs_15.0_HMM.LIB" \
  "/container/output" \
  4

# Step 2: Consolidate results
docker run --rm --platform linux/amd64 \
  -v "$(pwd)/output:/container/output:rw" \
  -v "$(pwd)/db/tigrfam:/container/db:ro" \
  tigrfam-annotation:1.0 \
  python3 /container/scripts/consolidate_tigrfam.py \
  "genome_name" \
  "/container/output/tigrfam_domains.tsv" \
  "/container/db/tigrfam_definitions.txt" \
  "/container/output/proteins.faa" \
  "/container/output/tigrfam.tsv"
```

## Database Statistics

TIGRfam v15.0:
- ~4,500 protein families
- Focus on prokaryotic proteins
- Equivalogs: Precise functional assignment
- Subfamilies: Broader functional categories
- Conservative, high-confidence annotations

## Comparison with Pfam

| Database | Families | Coverage | Specificity |
|----------|----------|----------|-------------|
| TIGRfam  | ~4,500   | Lower    | Higher      |
| Pfam     | ~20,000  | Higher   | Variable    |

TIGRfam is more conservative and prokaryote-focused, while Pfam has broader coverage across all domains of life.

## Building Container

```bash
cd processing/containers/tigrfam
./build_container.sh
```

## Testing

```bash
# Test with reference genome
./annotate-genome --select test --tools tigrfam -t 4

# Expected output:
# - output/tigrfam/test/tigrfam.tsv (consolidated)
# - output/tigrfam/test/native/tigrfam_domains.tsv (HMMER raw)
# - output/tigrfam/test/native/hmmscan_output.txt (full HMMER output)
```

## Troubleshooting

### Database not found
```
Error: TIGRfam database not found: db/tigrfam/TIGRFAMs_15.0_HMM.LIB
```
**Solution:** Download TIGRfam database (see Database Setup above)

### Database not pressed
```
Warning: Database not pressed (missing .h3i file)
```
**Solution:** Run `hmmpress db/tigrfam/TIGRFAMs_15.0_HMM.LIB`

### Missing definitions
```
Warning: Definitions file not found
```
**Solution:** Create tigrfam_definitions.txt (see Database Setup above)

## Dependencies

- Docker with linux/amd64 platform support
- TIGRfam database (~45 MB + indices)
- RASTtk gene calls (protein sequences)

## Integration

TIGRfam is integrated into the main annotation pipeline:
- Automatically uses existing gene calls
- Can run standalone or with other tools
- Outputs follow standard format with confidence scoring
- Results consolidated into master annotation table
