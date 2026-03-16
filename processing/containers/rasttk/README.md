# BV-BRC RASTtk Container

Independent container for gene calling and functional annotation using BV-BRC RASTtk toolkit.

## Container Design

**Philosophy:** Self-contained, lightweight container with only BV-BRC tools.
- **Size:** ~500MB (no databases included)
- **Platform:** linux/amd64 (compatible with macOS ARM via Rosetta)
- **Base:** Ubuntu 22.04

## What's Included

### BV-BRC Tools
- **rast-create-genome** - Create Genome Typed Object (GTO)
- **rast-call-features-*** - All RASTtk feature calling tools
  - rRNA genes (SEED)
  - tRNA genes (tRNAscan-SE)
  - Repeat regions
  - CDS (Prodigal, Glimmer3)
  - Functional annotation (kmer, similarity)
- **rast-export-genome** - Export to GFF, GenBank, FASTA formats

### Processing Scripts
Located in `/container/scripts/`:
- `run_rasttk_incremental.sh` - Main RASTtk pipeline (15 steps)
- `run_with_bvbrc.sh` - BV-BRC command wrapper
- `generate_rast_tsv.py` - SEED subsystem enrichment
- `consolidate_rast_output.py` - Output formatter
- `finalize_structure.py` - Directory organizer
- `process_genome.sh` - Single genome processor
- `run_batch_genomes.sh` - Batch processor
- `pipeline_status.sh` - Status checker

## What's External (Mounted)

### Databases
**Location:** Mounted from `db/rasttk/` → `/container/db/`
- `subsystem_mapping.tsv` - SEED subsystem hierarchy (1,055 subsystems)
- `variant_definitions/` - Complete SEED variant definitions (~450MB)

If `variant_definitions/` is missing, the pipeline still completes and writes placeholder
SEED variant fields in `rast.tsv` (status: `pending`) so downstream steps remain non-blocking.

**Why external?**
- Keeps container small
- Allows database updates without rebuild
- Shared across multiple containers if needed

### Input/Output
- **Input:** Mounted at `/container/input/` (read-only)
- **Output:** Mounted at `/container/output/` (read-write)

## Build

```bash
cd processing/containers/rasttk
./build_container.sh
```

**Result:** `rasttk-annotation:1.0` image

## Usage

### Via Wrapper Script (Recommended)
```bash
cd processing/scripts/rasttk
./run_annotation.sh \
    --input genome.fasta \
    --output ../../output/rasttk \
    --name my_genome \
    --scientific "Escherichia coli" \
    --threads 4
```

### Direct Container Run
```bash
docker run --rm \
    -v $(pwd)/input:/container/input:ro \
    -v $(pwd)/output:/container/output:rw \
    -v $(pwd)/db/rasttk:/container/db:ro \
    rasttk-annotation:1.0 \
    /container/scripts/run_rasttk_incremental.sh \
        genome_name \
        /container/input/genome.fasta \
        /container/output \
        "Organism name" \
        11 \
        Bacteria
```

### Using run_container.sh
```bash
./run_container.sh \
    --input /path/to/input \
    --output /path/to/output \
    --db /path/to/db \
    bash -c "command to run"
```

## Output Structure

```
output/rasttk/<genome_name>/
├── gene_call/
│   ├── <genome>.faa       # Protein sequences
│   ├── <genome>.gff       # GFF3 format
│   └── <genome>.gbk       # GenBank format
├── processed/
│   └── <genome>_rast.tsv  # Consolidated annotation with SEED enrichment
└── intermediate/
    └── <genome>.gto.*     # Genome Typed Objects (15 steps)
```

## Annotation Output Format

**File:** `<genome>_rast.tsv`

**Columns:**
1. `feature_id` - Feature identifier
2. `feature_type` - CDS, rRNA, tRNA, etc.
3. `contig_id` - Contig/scaffold ID
4. `start` - Start position
5. `end` - End position
6. `strand` - Plus (+) or minus (-)
7. `length_bp` - Feature length
8. `length_aa` - Protein length (CDS only)
9. `function` - Functional annotation
10. `ec_numbers` - EC numbers (`;` separated)
11. `subsystems` - SEED subsystems (`;` separated)
12. `subsystem_classes` - Hierarchy (`;` separated)
13. `roles` - SEED roles (`;` separated)
14. `aliases` - Alternative IDs
15. `sequence` - Nucleotide/amino acid sequence

## Thread Usage

RASTtk doesn't parallelize internally (single-threaded BV-BRC tools), but setting `--threads` allows:
- Genome-level parallelization (when processing multiple genomes)
- Resource allocation hints for downstream tools

## Dependencies

None! Container is fully self-contained with:
- Ubuntu 22.04 base
- Python 3.10
- BV-BRC CLI 1.040 (.deb package)
- All required Perl modules
- System libraries

## Size Optimization

**Container only:** ~500MB
- BV-BRC CLI: ~300MB
- Python + dependencies: ~150MB
- System packages: ~50MB

**External databases:** ~450MB (not in container)
- Mounted at runtime
- No container bloat

## Platform Notes

**macOS ARM64 (M1/M2/M3):**
- Uses `--platform linux/amd64`
- Runs via Rosetta 2 emulation
- ~10-20% performance overhead
- Fully functional

**Linux x86_64:**
- Native performance
- No emulation needed

**Windows WSL2:**
- Requires Docker Desktop
- Use WSL2 paths for volumes

## Troubleshooting

### Container not found
```bash
cd processing/containers/rasttk
./build_container.sh
```

### Permission errors
```bash
chmod +x processing/containers/rasttk/*.sh
chmod +x processing/scripts/rasttk/*.sh
```

### Database not found
```bash
# Ensure db/rasttk/ contains:
ls -la db/rasttk/
# Should show:
# - subsystem_mapping.tsv
# - variant_definitions/
```

### BV-BRC tools not found
```bash
# Test container:
docker run --rm rasttk-annotation:1.0 which rast-create-genome
# Should output: /usr/bin/rast-create-genome
```

## Integration

### Called by Main Pipeline
The main `annotate-genome` command calls:
```bash
processing/scripts/rasttk/run_annotation.sh
```

This wrapper:
1. Validates inputs
2. Builds container if needed
3. Runs RASTtk annotation
4. Generates consolidated TSV
5. Returns to main pipeline

### Standalone Usage
Can be used independently:
```bash
cd processing/scripts/rasttk
./run_annotation.sh -i genome.fasta -o output/
```

## Version

- **Container:** 1.0
- **BV-BRC CLI:** 1.040
- **Updated:** February 11, 2026

## Next Tools

After RASTtk, the pipeline continues with:
1. **COG** - Clusters of Orthologous Groups  
2. **Pfam** - Protein families
3. **KEGG** - Metabolic pathways
4. **dbCAN** - Carbohydrate-active enzymes
5. ...and more

Each tool has its own independent container following the same design principles.
