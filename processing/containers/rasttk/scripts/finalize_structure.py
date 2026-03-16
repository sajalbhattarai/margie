
import os

def generate_documentation(output_dir, genome_name, scientific_name):
    """
    Generates FILE_GUIDE.md and README_rast_tsv.md in the output directory,
    replacing 'theta' with the genome_name and customizing content.
    """
    
    # 1. FILE_GUIDE.md
    file_guide_content = f"""# RASTtk Output Files Explained

## Quick Reference Table

| File | Format | Content | Use For |
|------|--------|---------|---------|
| **{genome_name}_features.tsv** | TSV | ALL features (CDS, RNA, repeats, prophage, etc.) | Complete feature list |
| **{genome_name}_CDS.tsv** | TSV | Only protein-coding genes | Protein analysis |
| **{genome_name}_RNA.tsv** | TSV | Only RNA genes (rRNA, tRNA) | RNA analysis |
| **{genome_name}_prophage.tsv** | TSV | Only prophage elements | Phage analysis |
| **{genome_name}_summary.txt** | Text | Human-readable summary with counts | Quick overview |
| **{genome_name}.faa** | FASTA | Protein sequences (amino acids) | **Input for annotation tools** |
| **{genome_name}.ffn** | FASTA | Gene sequences (nucleotides) | Gene expression analysis |
| **{genome_name}.fna** | FASTA | Contig sequences | Genome assembly |
| **{genome_name}.gff** | GFF3 | Gene features (standard format) | Genome browsers, tools |
| **{genome_name}.gbk** | GenBank | Complete annotation | NCBI submission, visualization |
| **{genome_name}.gto.14** | JSON | Genome Typed Object (final) | Further RASTtk processing |

---

## Detailed Explanations

### 1. TSV Files (Tab-Separated Values)

#### **{genome_name}_features.tsv** - Complete Feature List
- **Contains**: ALL features annotated
  - Protein-coding genes (CDS)
  - rRNA genes
  - tRNA genes
  - Repeat regions
  - CRISPR elements
  - Prophage regions
  - Selenoproteins, pyrrolysoproteins
- **Columns**: Feature ID, Location, Type, Function, Aliases
- **Use**: Complete inventory of all genomic features

#### **{genome_name}_CDS.tsv** - Protein-Coding Genes Only
- **Contains**: Only CDS (protein-coding) features
- **Use**: When you only need protein-coding genes
- **Same columns as features.tsv** but filtered

#### **{genome_name}_RNA.tsv** - RNA Genes Only
- **Contains**: Only RNA features (rRNA, tRNA)
- **Use**: RNA gene analysis, ribosome studies
- **Includes**: 16S, 23S rRNA; tRNA genes

#### **{genome_name}_prophage.tsv** - Prophage Elements Only
- **Contains**: Only prophage regions
- **Use**: Phage analysis, horizontal gene transfer studies
- **Identified by**: PhiSpy algorithm

---

### 2. FASTA Sequence Files

#### **{genome_name}.faa** - Protein Sequences ⭐ MOST IMPORTANT
- **Format**: FASTA amino acid sequences
- **Contains**: Translations of all protein-coding genes
- **Use**: 
  - ✓ **INPUT FOR ALL DOWNSTREAM ANNOTATION TOOLS**
  - COG, eggNOG, Pfam, dbCAN, KEGG, etc.

#### **{genome_name}.ffn** - Gene Nucleotide Sequences
- **Format**: FASTA DNA sequences
- **Contains**: Nucleotide sequences of all genes
- **Use**: 
  - Gene expression analysis
  - Codon usage studies
  - PCR primer design

#### **{genome_name}.fna** - Contig Sequences
- **Format**: FASTA DNA sequences
- **Contains**: Complete genome contigs (same as input)
- **Use**: 
  - Reference genome
  - Alignment reference
  - Assembly validation

---

### 3. Standard Annotation Formats

#### **{genome_name}.gff** - Gene Feature Format
- **Format**: GFF3 (standard bioinformatics format)
- **Contains**: Gene locations, types, attributes
- **Use**:
  - Genome browsers (IGV, Artemis, JBrowse)
  - Bioinformatics pipelines
  - Cross-tool compatibility

#### **{genome_name}.gbk** - GenBank Format
- **Format**: GenBank flat file
- **Contains**: Complete annotation with all metadata
- **Use**:
  - NCBI GenBank submission
  - Visualization tools (Artemis, Geneious)
  - Comprehensive archival format

---

### 4. Processing Files

#### **{genome_name}.gto.1** through **{genome_name}.gto.14** - Genome Typed Objects
- **Format**: JSON (KBase-compatible)
- **Contains**: Incremental annotation data
- **Keep**: Only gto.14 (final) is needed; others are intermediate

#### **{genome_name}_summary.txt** - Human-Readable Summary
- **Format**: Plain text
- **Contains**: 
  - Feature counts
  - File list
  - Annotation steps completed
- **Use**: Quick overview, documentation

---

## Which Files to Use When?

### For Downstream Annotation Tools (COG, Pfam, etc.)
✓ Use: **{genome_name}.faa** (already copied to `output/prodigal/{genome_name}/{genome_name}/{genome_name}.faa`)

### For Genome Visualization
✓ Use: **{genome_name}.gff** or **{genome_name}.gbk**

### For Data Analysis/Statistics
✓ Use: **{genome_name}_features.tsv** (all features)
✓ Use: **{genome_name}_CDS.tsv** (proteins only)

### For Publication/Submission
✓ Use: **{genome_name}.gbk** (GenBank format)
✓ Use: **{genome_name}.gff** (standard format)

### For Further RASTtk Processing
✓ Use: **{genome_name}.gto.14** (final GTO)
"""
    
    with open(os.path.join(output_dir, "FILE_GUIDE.md"), "w") as f:
        f.write(file_guide_content)

    # 2. README_rast_tsv.md
    readme_content = f"""# RAST.TSV Output File

## Overview
The `rast.tsv` file consolidates all RASTtk annotation outputs into a single tab-separated file for easy merging with downstream annotation results.

## Column Structure

| Column | Description | Example | Purpose |
|--------|-------------|---------|---------|
| `organism` | Scientific name | {scientific_name} | Organism identifier |
| `feature_id` | RASTtk feature ID | fig|... | **Key column for merging with other annotations** |
| `RAST_gene_id` | Gene location string | ... | Contig and position info |
| `gene_start` | Start position | ... | Gene start coordinate |
| `gene_end` | End position | ... | Gene end coordinate |
| `na_length` | Nucleotide length | ... | Length of nucleotide sequence (bp) |
| `aa_length` | Amino acid length | ... | Length of protein sequence (aa) |
| `na_seq` | Nucleotide sequence | ... | Full gene DNA sequence |
| `aa_seq` | Amino acid sequence | ... | Full protein sequence |
| `RAST_gene_type` | Feature type | CDS, rRNA, tRNA, prophage | Type of genomic feature |
| `RAST_description` | Function annotation | ... | RAST-assigned function |
| `RAST_EC_numbers` | EC numbers | ... | Extracted enzyme classification |
| `RAST_strand` | DNA strand | +, - | Coding strand orientation |
| `RAST_feature_hash` | Feature hash | ... | MD5 hash for tracking |

## Key Points

### Merging with Other Annotations
- **Use the `feature_id` column** to merge with downstream annotation results (COG, Pfam, eggNOG, etc.)
- The `feature_id` matches exactly with the sequence IDs in `{genome_name}.faa` file

### RAST Prefix Convention
- Columns prefixed with `RAST_` contain RAST-specific annotations
- Non-prefixed columns (`organism`, `feature_id`, `gene_start`, `gene_end`, `na_length`, `aa_length`, `na_seq`, `aa_seq`) are common identifiers for merging

### EC Numbers
- Automatically extracted from RAST descriptions
- Multiple EC numbers separated by semicolons
- Empty if no EC number found

### Gene Location Parsing
- `RAST_gene_id` format: `CONTIG_START[+/-]LENGTH`
- `+` indicates forward strand, `-` indicates reverse strand
- `gene_start` and `gene_end` are absolute coordinates on the contig

## Generation
This file is automatically generated by:
```bash
python processing/scripts/rasttk/consolidate_rast_output.py output/rasttk/{genome_name} "{scientific_name}"
```

Or as part of the full pipeline:
```bash
./processing/scripts/rasttk/run_rasttk_incremental.sh
```

## Usage in Downstream Analysis
```python
import pandas as pd

# Load RAST annotations
rast = pd.read_csv('output/rasttk/{genome_name}/rast.tsv', sep='\\t')

# Load other annotations (e.g., COG)
cog = pd.read_csv('output/cog/{genome_name}/cog_results.tsv', sep='\\t')

# Merge on feature_id
merged = rast.merge(cog, on='feature_id', how='left')
```
"""
    with open(os.path.join(output_dir, "README_rast_tsv.md"), "w") as f:
        f.write(readme_content)

    print(f"Generated documentation files for {genome_name}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 4:
        print("Usage: python finalize_structure.py <output_dir> <genome_name> <scientific_name>")
    else:
        generate_documentation(sys.argv[1], sys.argv[2], sys.argv[3])
