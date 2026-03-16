#!/usr/bin/env python3
"""
Consolidate Diamond-based EggNOG annotations

Parses Diamond blastp hits against bacteria.dmnd and looks up comprehensive
annotations from eggNOG database files including:
- Ortholog group (OG) assignments
- COG functional categories
- Functional descriptions
- Taxonomy information
"""

import sys
import csv
from datetime import datetime
from collections import defaultdict


def log(message):
    """Log with timestamp."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}]  {message}", flush=True)


def load_protein_lengths(fasta_file):
    """Load protein sequences and compute lengths."""
    protein_lengths = {}
    current_id = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_id:
                    protein_lengths[current_id] = len(''.join(current_seq))
                # Start new sequence
                header = line[1:].split()[0]
                current_id = header
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_id:
            protein_lengths[current_id] = len(''.join(current_seq))
    
    return protein_lengths


def load_eggnog_annotations(annotations_file):
    """
    Load EggNOG ortholog group annotations.
    Returns dict: og_id -> {category, description}
    """
    annotations = {}
    
    log(f"  Loading EggNOG OG annotations from: {annotations_file}")
    
    with open(annotations_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            
            # Format: tax_level, og_id, category, [description]
            tax_level = parts[0]
            og_id = parts[1]
            category = parts[2] if len(parts) > 2 else ''
            description = parts[3] if len(parts) > 3 else ''
            
            # We're using bacterial annotations (tax_level=2)
            if tax_level == '2':
                annotations[og_id] = {
                    'category': category,
                    'description': description
                }
    
    log(f"    Loaded {len(annotations)} ortholog group annotations")
    return annotations


def load_protein_to_og_mapping(members_file):
    """
    Load mapping from protein IDs to ortholog groups.
    Returns dict: protein_id -> og_id
    
    The members file format is:
    tax_level  og_id  count1  count2  comma_separated_protein_ids
    """
    protein_to_og = {}
    
    log(f"  Loading protein to OG mapping from: {members_file}")
    log("    This may take a moment (394MB file)...")
    
    lines_processed = 0
    with open(members_file, 'r') as f:
        for line in f:
            lines_processed += 1
            if lines_processed % 50000 == 0:
                log(f"    Processed {lines_processed:,} OGs, mapped {len(protein_to_og):,} proteins...")
            
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            
            tax_level = parts[0]
            og_id = parts[1]
            members_str = parts[4]  # Comma-separated protein IDs
            
            # Only process bacterial OGs (tax_level=2)
            if tax_level != '2':
                continue
            
            # Parse member proteins
            members = members_str.split(',')
            for member in members:
                member = member.strip()
                if member:
                    protein_to_og[member] = og_id
    
    log(f"    Loaded {len(protein_to_og):,} protein -> OG mappings")
    return protein_to_og


def load_taxonomy_names(taxid_file):
    """
    Load taxonomy ID to name mapping.
    Returns dict: taxid -> organism_name
    """
    taxid_to_name = {}
    
    log(f"  Loading taxonomy names from: {taxid_file}")
    
    try:
        with open(taxid_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    taxid = parts[0]
                    name = parts[1]
                    taxid_to_name[taxid] = name
        
        log(f"    Loaded {len(taxid_to_name):,} taxonomy names")
    except FileNotFoundError:
        log("    Warning: Taxonomy file not found, will use taxids only")
    
    return taxid_to_name


def parse_diamond_hits(diamond_file):
    """
    Parse Diamond output.
    Returns dict: protein_id -> best_hit_dict with all alignment details
    """
    hits = {}
    
    with open(diamond_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        for row in reader:
            protein_id = row['qseqid']
            
            # Take only best hit per protein (Diamond already provides max-target-seqs 1)
            if protein_id not in hits:
                hits[protein_id] = {
                    'og_id': row['sseqid'],  # Full protein ID from database
                    'pident': row.get('pident', ''),
                    'length': row.get('length', ''),
                    'evalue': row.get('evalue', ''),
                    'bitscore': row.get('bitscore', ''),
                    'qcovhsp': row.get('qcovhsp', ''),
                    'scovhsp': row.get('scovhsp', ''),
                    'qstart': row.get('qstart', ''),
                    'qend': row.get('qend', ''),
                    'sstart': row.get('sstart', ''),
                    'send': row.get('send', '')
                }
    
    return hits


def main():
    if len(sys.argv) != 7:
        print("Usage: consolidate_diamond_eggnog.py <diamond_hits.tsv> <input.faa> <annotations.tsv> <members.tsv> <organism> <output.tsv>")
        sys.exit(1)
    
    diamond_file = sys.argv[1]
    fasta_file = sys.argv[2]
    annotations_file = sys.argv[3]
    members_file = sys.argv[4]
    organism_name = sys.argv[5]
    output_tsv = sys.argv[6]
    
    # Try to load taxonomy names (optional)
    taxid_file = '/container/db/taxid_to_name.tsv'
    
    log("="*70)
    log("EggNOG Diamond-based Annotation Consolidation")
    log("")
    log(f"Organism:            {organism_name}")
    log(f"Diamond hits:        {diamond_file}")
    log(f"Input FASTA:         {fasta_file}")
    log(f"Annotations DB:      {annotations_file}")
    log(f"Members DB:          {members_file}")
    log(f"Taxonomy DB:         {taxid_file}")
    log(f"Output:              {output_tsv}")
    log("")
    
    # Load reference data
    log("Loading reference data...")
    log("  Loading protein lengths...")
    protein_lengths = load_protein_lengths(fasta_file)
    log(f"    Loaded {len(protein_lengths)} proteins")
    log("")
    
    log("  Loading EggNOG annotations...")
    og_annotations = load_eggnog_annotations(annotations_file)
    log("")
    
    log("  Loading protein to OG mappings...")
    protein_to_og = load_protein_to_og_mapping(members_file)
    log("")
    
    log("  Loading taxonomy names...")
    taxid_to_name = load_taxonomy_names(taxid_file)
    log("")
    
    log("  Parsing Diamond hits...")
    diamond_hits = parse_diamond_hits(diamond_file)
    log(f"    Found hits for {len(diamond_hits)} proteins")
    log("")
    
    # Write consolidated output
    log("Writing consolidated output...")
    
    # Output format with comprehensive information
    with open(output_tsv, 'w', newline='') as out_f:
        writer = csv.writer(out_f, delimiter='\t')
        
        # Header - following eggNOG-mapper naming conventions
        header = [
            'organism',
            'feature_id',
            'protein_length',
            'EGGNOG_hit_protein',
            'EGGNOG_hit_taxid',
            'EGGNOG_hit_organism',
            'EGGNOG_OG',
            'EGGNOG_COG_category',
            'EGGNOG_COG_description',
            'EGGNOG_description',
            'EGGNOG_pident',
            'EGGNOG_alignment_length',
            'EGGNOG_evalue',
            'EGGNOG_bitscore',
            'EGGNOG_query_coverage',
            'EGGNOG_subject_coverage',
            'EGGNOG_query_start',
            'EGGNOG_query_end',
            'EGGNOG_subject_start',
            'EGGNOG_subject_end'
        ]
        writer.writerow(header)
        
        rows_written = 0
        proteins_with_hits = 0
        proteins_with_og = 0
        proteins_with_annotation = 0
        
        for protein_id in sorted(protein_lengths.keys()):
            prot_len = protein_lengths[protein_id]
            
            if protein_id in diamond_hits:
                hit = diamond_hits[protein_id]
                proteins_with_hits += 1
                
                # Extract hit protein ID
                hit_protein = hit['og_id']  # This is actually the full protein ID from Diamond
                
                # Extract taxid from protein ID (format: taxid.protein)
                hit_taxid = ''
                hit_organism = ''
                if '.' in hit_protein:
                    hit_taxid = hit_protein.split('.')[0]
                    hit_organism = taxid_to_name.get(hit_taxid, f"taxid:{hit_taxid}")
                
                # Look up OG assignment for the hit protein
                og_id = protein_to_og.get(hit_protein, '')
                
                if og_id:
                    proteins_with_og += 1
                
                # Look up OG annotation
                category = ''
                og_description = ''
                cog_desc = ''
                
                if og_id and og_id in og_annotations:
                    proteins_with_annotation += 1
                    annot = og_annotations[og_id]
                    category = annot['category']
                    og_description = annot['description']
                    cog_desc = get_cog_category_description(category)
                
                row = [
                    organism_name,
                    protein_id,
                    prot_len,
                    hit_protein,
                    hit_taxid,
                    hit_organism,
                    og_id,
                    category,
                    cog_desc,
                    og_description,
                    hit['pident'],
                    hit['length'],
                    hit['evalue'],
                    hit['bitscore'],
                    hit['qcovhsp'],
                    hit['scovhsp'],
                    hit['qstart'],
                    hit['qend'],
                    hit['sstart'],
                    hit['send']
                ]
                
                writer.writerow(row)
                rows_written += 1
    
    log(f"  Wrote {rows_written} annotation rows")
    log("")
    
    log("Results summary:")
    log(f"  Total proteins:              {len(protein_lengths)}")
    log(f"  Proteins with Diamond hits:  {proteins_with_hits}")
    log(f"  Proteins with OG assignment: {proteins_with_og}")
    log(f"  Proteins with annotations:   {proteins_with_annotation}")
    if len(protein_lengths) > 0:
        hit_coverage = 100 * proteins_with_hits / len(protein_lengths)
        og_coverage = 100 * proteins_with_og / len(protein_lengths)
        annot_coverage = 100 * proteins_with_annotation / len(protein_lengths)
        log(f"  Diamond hit coverage:        {hit_coverage:.1f}%")
        log(f"  OG assignment coverage:      {og_coverage:.1f}%")
        log(f"  Annotation coverage:         {annot_coverage:.1f}%")
    log(f"  Output rows:                 {rows_written}")
    log("")
    log("✓ EggNOG consolidation complete!")


def get_cog_category_description(category):
    """Get description for COG functional category."""
    cog_categories = {
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
    
    if not category:
        return ''
    
    # Handle multiple categories (e.g., "KL")
    descriptions = []
    for cat in category:
        if cat in cog_categories:
            descriptions.append(cog_categories[cat])
    
    return '; '.join(descriptions) if descriptions else category


if __name__ == '__main__':
    main()
