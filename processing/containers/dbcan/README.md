# dbCAN CAZyme Annotation Container

## Overview

This container provides CAZyme (Carbohydrate-Active enZyme) annotation using the dbCAN database and HMMER for sequence homology search.

**Version**: 1.0  
**Base Image**: Ubuntu 22.04  
**Tools**: HMMER 3.4, Python 3  
**Database**: dbCAN-HMMdb-V12 (external mount)  
**Size**: ~168 MB  

## Architecture

The container follows a lean design pattern:
- **No embedded databases**: Database mounted from workspace `db/dbcan/`
- **Minimal dependencies**: Only HMMER 3.4 and Python 3
- **Standalone scripts**: All annotation logic inside container
- **Standardized output**: 20-column TSV with CAZyme annotations

## Contents

```
dbcan/
├── Dockerfile                      # Container definition
├── build_container.sh              # Build script
├── README.md                       # This file
└── scripts/
    ├── run_hmmscan.sh              # HMMER execution script
    └── consolidate_dbcan.py        # Output consolidation script
```

## Database Requirements

The container expects the following files in the mounted database directory (`/container/db/`):

1. **dbCAN-HMMdb-V12.txt** (109 MB)  
   - Main HMM database with CAZyme family profiles
   - Should be pressed with `hmmpress` for optimal performance
   - Pressed files: `.h3f`, `.h3i`, `.h3m`, `.h3p`

2. **CAZyDB.08062022.fam.subfam.ec.txt** (651 KB)  
   - Family/subfamily to EC number mappings
   - Format: `subfamily \t protein_id \t ec_number`

3. **substrate-mappings.tsv** (72 KB)  
   - Substrate predictions for CAZyme families
   - Format: `Substrate_high \t Substrate_simple \t Family \t Name \t EC`

## Building

```bash
cd processing/containers/dbcan
./build_container.sh
```

Build time: ~2-3 minutes  
Expected size: ~168 MB

## Usage

### Direct Container Usage

```bash
# Run HMMER scan
docker run --rm \
    -v /path/to/proteins.faa:/container/input/proteins.faa:ro \
    -v /path/to/output:/container/output:rw \
    -v /path/to/db/dbcan:/container/db:ro \
    dbcan-annotation:1.0 \
    /bin/bash /container/scripts/run_hmmscan.sh \
        --input /container/input/proteins.faa \
        --output /container/output \
        --database /container/db/dbCAN-HMMdb-V12.txt \
        --threads 4

# Consolidate results
docker run --rm \
    -v /path/to/output:/container/output:rw \
    -v /path/to/db/dbcan:/container/db:ro \
    dbcan-annotation:1.0 \
    python3 /container/scripts/consolidate_dbcan.py \
        organism_name \
        /container/output/native/dbcan_domains.tsv \
        /container/db/CAZyDB.08062022.fam.subfam.ec.txt \
        /container/db/substrate-mappings.tsv \
        /container/output/gene_calls/organism.faa \
        /container/output/dbcan.tsv
```

### Via Pipeline Wrapper

```bash
# Using wrapper script (recommended)
bash processing/scripts/dbcan/run_annotation.sh \
    --input input/proteins.faa \
    --output output/dbcan \
    --name test_genome \
    --threads 4

# Via annotate-genome pipeline (most convenient)
./annotate-genome --select test --tools dbcan -t 4
```

## Input

- **Protein FASTA file** (`.faa`)
  - Multi-FASTA format
  - Protein sequences from gene calling (RASTtk, Prodigal, etc.)
  - Sequence IDs used for feature mapping

Example:
```
>fig|6666666.12345.peg.1
MTKIAVLGISGSGKSTLLNSILGIEKPTSGIILIDGKEIKKLGFEDWKGDKWGFNLEEKN...
>fig|6666666.12345.peg.2
MVKLNDRYTLVLDEHHFTKEYSQLTKEQGLEFAKKYKEMKLEKDVLMYHFVSNGKSKTK...
```

## Output Structure

```
output_dir/
└── genome_name/
    ├── native/
    │   ├── dbcan_domains.tsv        # HMMER domtblout format (18K)
    │   └── hmmscan_output.txt       # Full HMMER output (3.0M)
    ├── gene_calls/
    │   └── genome_name.faa          # Copy of input proteins
    └── dbcan.tsv                    # Consolidated annotation (137K)
```

### Output Format: dbcan.tsv

20-column tab-separated file with CAZyme annotations:

| Column | Description |
|--------|-------------|
| organism | Organism name |
| feature_id | Protein/gene identifier |
| protein_length | Amino acid length |
| DBCAN_family_count | Number of CAZyme domains found |
| DBCAN_families | Family assignments (roman numeral indexed) |
| DBCAN_family_names | Common names for families |
| DBCAN_domain_positions | Alignment positions (start-end) |
| DBCAN_substrates | Predicted substrates |
| DBCAN_EC_numbers | Associated EC numbers |
| DBCAN_domain_evalues | Domain E-values (independent) |
| DBCAN_full_evalues | Full sequence E-values |
| DBCAN_domain_scores | Domain bit scores |
| DBCAN_coverages | HMM coverage values |
| DBCAN_best_family | Best scoring family |
| DBCAN_best_evalue | Best E-value |
| DBCAN_best_score | Best bit score |
| DBCAN_confidence | Composite confidence (0-1) |
| DBCAN_evalue_score | E-value component |
| DBCAN_coverage_score | Coverage component |
| DBCAN_bitscore_score | Bit score component |

**Example rows:**
```
organism  feature_id                    protein_length  DBCAN_family_count  DBCAN_families  ...
test      fig|6666666.1481409.peg.1129  341             1                   (i) GH3         ...
test      fig|6666666.1481409.peg.1184  653             2                   (i) GH153; (ii) CE4  ...
```

## Confidence Scoring

Multi-component confidence calculation (0-1 scale):

1. **E-value confidence** (50% weight):
   - < 1e-10: 1.0 (very high)
   - < 1e-5: 0.8 (high)
   - < 1e-3: 0.5 (medium)
   - ≥ 1e-3: 0.2 (low)

2. **Coverage confidence** (30% weight):
   - HMM profile coverage (0-1)

3. **Score confidence** (20% weight):
   - Normalized bit score (score / 100, capped at 1.0)

**Final confidence** = 0.5×E-value + 0.3×coverage + 0.2×score

## CAZyme Families

The dbCAN database contains profiles for:
- **GH**: Glycoside Hydrolases (breakdown glycosidic bonds)
- **GT**: GlycosylTransferases (form glycosidic bonds)
- **PL**: Polysaccharide Lyases (eliminate glycosidic bonds)
- **CE**: Carbohydrate Esterases (remove ester linkages)
- **AA**: Auxiliary Activities (redox enzymes)
- **CBM**: Carbohydrate-Binding Modules (binding domains)

Example: `GH3` = Glycoside Hydrolase family 3 (β-glucosidases, β-xylosidases)

## Performance

**Test genome** (2,196 proteins, 921 KB):
- HMMER scan: ~38 seconds (4 threads)
- Consolidation: ~1 second
- **Total runtime**: ~41 seconds

**Resource usage**:
- Memory: ~500 MB
- CPU: Scales linearly with threads
- Disk: ~4 MB per genome (output)

**Database size**:
- HMM profiles: 109 MB (uncompressed)
- Pressed indices: 118 MB (.h3f/.h3i/.h3m/.h3p)
- Metadata: 723 KB

## Statistics (Test Genome)

```
Total proteins: 2,196
Proteins with CAZyme domains: 39 (1.8%)
Total CAZyme domain hits: 91
Unique CAZyme families found: 62
Unique substrates predicted: 52
Average domains per CAZyme protein: 2.33
```

## Integration

The container integrates with the annotate-genome pipeline:

```bash
# Check available genomes
./annotate-genome --list

# Annotate specific genome
./annotate-genome --select test --tools dbcan -t 4

# Multiple genomes
./annotate-genome --select "test,ref_theta" --tools dbcan -t 8

# All genomes
./annotate-genome --select all --tools dbcan -t 16
```

## Citations

**dbCAN2**:  
Zhang H, Yohe T, Huang L, et al. (2018). dbCAN2: a meta server for automated carbohydrate-active enzyme annotation. *Nucleic Acids Research*, 46(W1):W95-W101. doi:10.1093/nar/gky418

**CAZy Database**:  
Cantarel BL, Coutinho PM, Rancurel C, et al. (2009). The Carbohydrate-Active EnZymes database (CAZy): an expert resource for Glycogenomics. *Nucleic Acids Research*, 37:D233-D238. doi:10.1093/nar/gkn663

**HMMER**:  
Eddy SR (2011). Accelerated profile HMM searches. *PLoS Computational Biology*, 7(10):e1002195. doi:10.1371/journal.pcbi.1002195

## Troubleshooting

**Problem**: HMMER scan is slow  
**Solution**: Press the HMM database with `hmmpress db/dbcan/dbCAN-HMMdb-V12.txt`

**Problem**: No substrates predicted  
**Solution**: Ensure `substrate-mappings.tsv` exists in database directory

**Problem**: Missing EC numbers  
**Solution**: Ensure `CAZyDB.08062022.fam.subfam.ec.txt` exists in database directory

**Problem**: Container build fails  
**Solution**: Check Docker version, ensure internet access for HMMER download

## Notes

- The dbCAN database is updated periodically. Current version is V12 (2023).
- Some CAZyme families may not have substrate predictions.
- Multiple domains in a single protein are common (modular enzymes).
- E-value cutoff default is 1e-15 (very stringent).
- Coverage cutoff default is 0.35 (35% of HMM must align).

## Version History

- **1.0** (2026-02-11): Initial release with HMMER 3.4, dbCAN-HMMdb-V12
