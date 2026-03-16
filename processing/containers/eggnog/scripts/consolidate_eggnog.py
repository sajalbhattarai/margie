#!/usr/bin/env python3
"""
Consolidate eggNOG-mapper output into comprehensive TSV format.
Parses official eggnog-mapper annotations and formats for the pipeline.
"""

import sys
import csv
from pathlib import Path
from datetime import datetime


def log(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}]  {message}", flush=True)


def load_protein_lengths(fasta_file):
    """Calculate protein lengths from input FASTA file."""
    proteins = {}
    current_id = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    proteins[current_id] = len(''.join(current_seq))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            proteins[current_id] = len(''.join(current_seq))
    
    return proteins


def parse_eggnog_annotations(annotations_file):
    """
    Parse eggnog-mapper .emapper.annotations file.
    Returns dict: protein_id -> annotation_dict
    """
    annotations = {}
    
    with open(annotations_file, 'r') as f:
        for line in f:
            # Skip header and comment lines
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 21:
                continue
            
            protein_id = parts[0]
            
            # Parse eggnog-mapper columns
            # Standard eggnog-mapper output format (v2.1+):
            # 0: query
            # 1: seed_ortholog
            # 2: evalue
            # 3: score
            # 4: eggNOG_OGs
            # 5: max_annot_lvl
            # 6: COG_category
            # 7: Description
            # 8: Preferred_name
            # 9: GOs
            # 10: EC
            # 11: KEGG_ko
            # 12: KEGG_Pathway
            # 13: KEGG_Module
            # 14: KEGG_Reaction
            # 15: KEGG_rclass
            # 16: BRITE
            # 17: KEGG_TC
            # 18: CAZy
            # 19: BiGG_Reaction
            # 20: PFAMs
            
            annotations[protein_id] = {
                'seed_ortholog': parts[1] if len(parts) > 1 else '',
                'evalue': parts[2] if len(parts) > 2 else '',
                'score': parts[3] if len(parts) > 3 else '',
                'eggNOG_OGs': parts[4] if len(parts) > 4 else '',
                'max_annot_lvl': parts[5] if len(parts) > 5 else '',
                'COG_category': parts[6] if len(parts) > 6 else '',
                'description': parts[7] if len(parts) > 7 else '',
                'preferred_name': parts[8] if len(parts) > 8 else '',
                'GOs': parts[9] if len(parts) > 9 else '',
                'EC': parts[10] if len(parts) > 10 else '',
                'KEGG_ko': parts[11] if len(parts) > 11 else '',
                'KEGG_pathway': parts[12] if len(parts) > 12 else '',
                'KEGG_module': parts[13] if len(parts) > 13 else '',
                'KEGG_reaction': parts[14] if len(parts) > 14 else '',
                'KEGG_rclass': parts[15] if len(parts) > 15 else '',
                'BRITE': parts[16] if len(parts) > 16 else '',
                'KEGG_TC': parts[17] if len(parts) > 17 else '',
                'CAZy': parts[18] if len(parts) > 18 else '',
                'BiGG_reaction': parts[19] if len(parts) > 19 else '',
                'PFAMs': parts[20] if len(parts) > 20 else '',
            }
    
    return annotations


def get_cog_category_description(category_code):
    """Convert COG category code to full description."""
    category_map = {
        'J': 'Translation, ribosomal structure and biogenesis',
        'A': 'RNA processing and modification',
        'K': 'Transcription',
        'L': 'Replication, recombination and repair',
        'B': 'Chromatin structure and dynamics',
        'D': 'Cell cycle control, cell division, chromosome partitioning',
        'Y': 'Nuclear structure',
        'V': 'Defense mechanisms',
        'T': 'Signal transduction mechanisms',
        'M': 'Cell wall/membrane/envelope biogenesis',
        'N': 'Cell motility',
        'Z': 'Cytoskeleton',
        'W': 'Extracellular structures',
        'U': 'Intracellular trafficking, secretion, and vesicular transport',
        'O': 'Posttranslational modification, protein turnover, chaperones',
        'X': 'Mobilome: prophages, transposons',
        'C': 'Energy production and conversion',
        'G': 'Carbohydrate transport and metabolism',
        'E': 'Amino acid transport and metabolism',
        'F': 'Nucleotide transport and metabolism',
        'H': 'Coenzyme transport and metabolism',
        'I': 'Lipid transport and metabolism',
        'P': 'Inorganic ion transport and metabolism',
        'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
        'R': 'General function prediction only',
        'S': 'Function unknown',
    }
    
    # Handle multiple categories
    if len(category_code) == 0:
        return ''
    elif len(category_code) == 1:
        return category_map.get(category_code, '')
    else:
        # Multiple categories
        descriptions = [category_map.get(c, c) for c in category_code]
        return '; '.join(descriptions)


def main():
    if len(sys.argv) < 5:
        log("Usage: consolidate_eggnog.py <annotations_file> <input_faa> <organism_name> <output_tsv>")
        sys.exit(1)
    
    annotations_file = sys.argv[1]
    input_faa = sys.argv[2]
    organism_name = sys.argv[3]
    output_tsv = sys.argv[4]
    
    log("="*70)
    log("eggNOG Annotation Consolidation")
    log("="*70)
    log("")
    log(f"Organism:            {organism_name}")
    log(f"Annotations file:    {annotations_file}")
    log(f"Input FASTA:         {input_faa}")
    log(f"Output TSV:          {output_tsv}")
    log("")
    
    # Check if files exist
    if not Path(annotations_file).exists():
        log(f"ERROR: Annotations file not found: {annotations_file}")
        sys.exit(1)
    
    if not Path(input_faa).exists():
        log(f"ERROR: Input FASTA not found: {input_faa}")
        sys.exit(1)
    
    # Load data
    log("Loading reference data...")
    log("")
    
    log("  Loading protein lengths...")
    protein_lengths = load_protein_lengths(input_faa)
    log(f"    Loaded {len(protein_lengths)} proteins")
    
    log("  Parsing eggNOG annotations...")
    annotations = parse_eggnog_annotations(annotations_file)
    log(f"    Parsed {len(annotations)} annotations")
    log("")
    
    # Write consolidated output
    log("Writing consolidated output...")
    log(f"  Output file: {output_tsv}")
    log("")
    
    with open(output_tsv, 'w', newline='') as out_f:
        writer = csv.writer(out_f, delimiter='\t')
        
        # Write header
        header = [
            'organism',
            'feature_id',
            'protein_length',
            'eggNOG_seed_ortholog',
            'eggNOG_evalue',
            'eggNOG_score',
            'eggNOG_OGs',
            'eggNOG_max_annot_lvl',
            'eggNOG_COG_category',
            'eggNOG_COG_category_description',
            'eggNOG_description',
            'eggNOG_preferred_name',
            'eggNOG_GOs',
            'eggNOG_EC',
            'eggNOG_KEGG_ko',
            'eggNOG_KEGG_pathway',
            'eggNOG_KEGG_module',
            'eggNOG_KEGG_reaction',
            'eggNOG_KEGG_rclass',
            'eggNOG_BRITE',
            'eggNOG_KEGG_TC',
            'eggNOG_CAZy',
            'eggNOG_BiGG_reaction',
            'eggNOG_PFAMs',
        ]
        writer.writerow(header)
        
        # Write data rows
        rows_written = 0
        for protein_id in sorted(protein_lengths.keys()):
            prot_len = protein_lengths[protein_id]
            
            if protein_id in annotations:
                annot = annotations[protein_id]
                
                # Get COG category description
                cog_desc = get_cog_category_description(annot['COG_category'])
                
                row = [
                    organism_name,
                    protein_id,
                    prot_len,
                    annot['seed_ortholog'],
                    annot['evalue'],
                    annot['score'],
                    annot['eggNOG_OGs'],
                    annot['max_annot_lvl'],
                    annot['COG_category'],
                    cog_desc,
                    annot['description'],
                    annot['preferred_name'],
                    annot['GOs'],
                    annot['EC'],
                    annot['KEGG_ko'],
                    annot['KEGG_pathway'],
                    annot['KEGG_module'],
                    annot['KEGG_reaction'],
                    annot['KEGG_rclass'],
                    annot['BRITE'],
                    annot['KEGG_TC'],
                    annot['CAZy'],
                    annot['BiGG_reaction'],
                    annot['PFAMs'],
                ]
                writer.writerow(row)
                rows_written += 1
        
        log(f"  Wrote {rows_written} annotation rows")
        log("")
    
    log("Results summary:")
    log(f"  Total proteins:            {len(protein_lengths)}")
    log(f"  Proteins with annotations: {len(annotations)}")
    if len(protein_lengths) > 0:
        coverage = 100 * len(annotations) / len(protein_lengths)
        log(f"  Coverage:                  {coverage:.1f}%")
    log(f"  Output rows:               {rows_written}")
    log("")
    log("✓ eggNOG consolidation complete!")

if __name__ == "__main__":
    main()
