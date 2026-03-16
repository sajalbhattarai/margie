# KEGG Annotation Container

Container for KEGG Orthology (KO) annotation using KofamScan and HMMER.

## Overview

- **Tool**: KofamScan 1.3.0 (Ruby-based KEGG annotation)
- **Search**: HMMER 3.3.2 (profile HMM search)
- **Database**: KOfam profiles (mounted externally)
- **Size**: 208 MB
- **Platform**: linux/amd64

## Architecture

**Lean Design**: Databases are **NOT** included in the container. They are mounted at runtime from the host's `db/kegg/` directory.

### Directory Structure

```
processing/containers/kegg/
├── Dockerfile              # Container definition
├── build_container.sh      # Build script
├── README.md              # This file
└── scripts/
    ├── generate_config.sh  # Creates config.yml for KofamScan
    ├── run_kofamscan.sh    # Runs KofamScan annotation
    └── consolidate_kegg.py # Consolidates and enriches output
```

### Container Layout

```
/container/
├── scripts/               # Processing scripts
│   ├── generate_config.sh
│   ├── run_kofamscan.sh
│   └── consolidate_kegg.py
├── db/                    # Database mount point (from host db/kegg/)
│   ├── profiles/          # HMM profiles (~24,000 KOs)
│   └── ko_list            # KO definitions and thresholds
├── input/                 # Input mount point
└── output/                # Output mount point
```

## Database Setup

The KEGG database must be set up in `db/kegg/` at the workspace root.

### Required Files

```
db/kegg/
├── profiles/              # Directory with ~24,000 .hmm files
│   ├── K00001.hmm
│   ├── K00002.hmm
│   └── ... (24,000+ files)
└── ko_list                # KO definitions (tab-separated)
```

### Download Database

```bash
# Create database directory
mkdir -p db/kegg
cd db/kegg

# Download KOfam profiles (~3GB compressed)
wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz

# Download KO list
wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz

# Extract
tar -xzf profiles.tar.gz
gunzip ko_list.gz

# Verify
ls profiles/*.hmm | wc -l  # Should show ~24,000 files
wc -l ko_list              # Should show ~24,000 lines
```

## Building the Container

```bash
cd processing/containers/kegg
./build_container.sh
```

This creates `kegg-annotation:1.0` image (~208 MB).

## Usage

### Via Main Pipeline

```bash
# From workspace root
./annotate-genome -i input/test.fasta --tools rasttk,kegg -t 8
```

This will:
1. Run RASTtk to generate gene calls (test.faa)
2. Run KEGG annotation on the proteins
3. Create consolidated output

### Via Direct Wrapper

```bash
# Requires RASTtk output first
processing/scripts/kegg/run_annotation.sh \
    -i output/rasttk/test/gene_calls/test.faa \
    -o output/kegg \
    -n test \
    -t 8
```

### Manual Container Execution

```bash
# From workspace root
docker run --rm \
    --platform linux/amd64 \
    -v "$(pwd)/output/rasttk/test/gene_calls/test.faa:/container/input/proteins.faa:ro" \
    -v "$(pwd)/output/kegg/test:/container/output:rw" \
    -v "$(pwd)/db/kegg:/container/db:ro" \
    kegg-annotation:1.0 \
    /container/scripts/run_kofamscan.sh \
        --input /container/input/proteins.faa \
        --output /container/output/native/kofamscan_output.txt \
        --threads 8
```

## Output Structure

```
output/kegg/
└── genome_name/
    ├── gene_calls/
    │   └── genome_name.faa      # Copy of protein sequences
    ├── native/
    │   └── kofamscan_output.txt # Raw KofamScan output (detail format)
    └── kegg.tsv                 # Consolidated annotation (18 columns)
```

## Output Format: kegg.tsv

Tab-separated file with 18 columns:

| Column | Description |
|--------|-------------|
| organism | Organism/genome name |
| feature_id | Protein feature ID |
| protein_length | Length in amino acids |
| KEGG_KO_count | Number of KO assignments |
| KEGG_KO_IDs | KO identifiers (semicolon-separated) |
| KEGG_KO_definitions | KO functional definitions |
| KEGG_EC_numbers | EC numbers (if annotated) |
| KEGG_scores | HMMER bit scores |
| KEGG_thresholds | KO-specific thresholds |
| KEGG_evalues | E-values |
| KEGG_score_types | Score types (full/domain) |
| KEGG_best_KO | Best KO assignment |
| KEGG_best_score | Highest score |
| KEGG_best_evalue | Best E-value |
| KEGG_evalue_score | Statistical strength (0-1) |
| KEGG_coverage_score | Score/threshold ratio (0-1) |
| KEGG_alignment_score | Not available (0) |
| KEGG_confidence | Overall confidence (0-1) |

### Confidence Scoring

- **evalue_score** (s1): Statistical strength, `-log10(E-value)/50`
- **coverage_score** (s2): Normalized score/threshold ratio
- **alignment_score** (s3): Not available from KofamScan (set to 0)
- **confidence**: `0.6*s1 + 0.3*s2 + 0.1*s3`

KofamScan uses F-measure optimized thresholds, so hits marked with `*` (above threshold) are high-confidence predictions.

## Example Output

```
organism  feature_id      KEGG_KO_count  KEGG_best_KO  KEGG_confidence
test      fig|123.1.peg.1 2              K00001        0.87
test      fig|123.1.peg.2 1              K00234        0.92
test      fig|123.1.peg.3 0                            
```

## Performance

- **Speed**: ~5-10 minutes for 2,000 proteins (8 threads)
- **Memory**: ~4-8 GB RAM
- **Database**: ~3 GB compressed, ~12 GB extracted

## Technical Details

### KofamScan

- Uses HMMER to search protein sequences against KOfam HMM profiles
- Each KO has a pre-computed F-measure optimized threshold
- Hits above threshold marked with `*` in detail output
- Supports parallel execution via GNU parallel

### HMMER 3.3.2

- Profile Hidden Markov Model search
- Optimized for SSE instructions (x86_64)
- Fast and sensitive sequence similarity search

### Database

- **KOfam**: KEGG Orthology HMM profiles
- **profiles/**: ~24,000 HMM files (one per KO)
- **ko_list**: Metadata with thresholds, definitions, EC numbers
- Updated regularly by KEGG team

## Integration

This container integrates with the main pipeline:

1. **RASTtk** generates gene calls (`.faa` proteins)
2. **KEGG** annotates proteins with KO numbers
3. **Consolidation** merges with other functional annotations
4. **Final output** includes KEGG pathways and modules

## Troubleshooting

### Database Not Found

```
ERROR: profiles directory not found in /container/db
```

**Solution**: Ensure database is set up at `db/kegg/` with `profiles/` and `ko_list`.

### No KofamScan Output

```
ERROR: Output file was not created!
```

**Solution**: Check HMMER is working: `docker run --rm kegg-annotation:1.0 hmmsearch -h`

### Low Hit Rate

If many proteins have no KO assignments:
- Check input is protein sequences (not nucleotides)
- Verify protein file has proper FASTA format
- Check proteins are from prokaryotes (KOfam optimized for bacteria/archaea)

## References

- **KofamScan**: https://www.genome.jp/tools/kofamkoala/
- **Official KOfamScan files**: https://www.genome.jp/ftp/tools/kofam_scan/
- **Official KOfamScan INSTALL**: https://www.genome.jp/ftp/tools/kofam_scan/INSTALL
- **KEGG**: https://www.kegg.jp/
- **HMMER**: http://hmmer.org/
- **KOfam Database**: ftp://ftp.genome.jp/pub/db/kofam/

## Official Config Compatibility

This container follows the official KOfamScan `config-template.yml` approach from Genome.jp INSTALL:

- It copies the official template when available.
- It rewrites only runtime paths (`profile`, `ko_list`, `hmmsearch`, `parallel`, `cpu`) to match mounted container paths.
- It does not modify upstream KOfamScan source files in-place.

Database mount compatibility:

- `/container/db/profiles` and `/container/db/ko_list`
- `/container/db/kegg/profiles` and `/container/db/kegg/ko_list`

This means your host database can live at either:

- `root-dir/db`
- `root-dir/db/kegg`

## Version

- Container: 1.0
- KofamScan: latest official release from `https://www.genome.jp/ftp/tools/kofam_scan/` at build time
- HMMER: 3.3.2
- Build date: 2026-02-11
