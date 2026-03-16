#!/usr/bin/env python3
"""
Consolidate COG annotation output into comprehensive TSV format.
Parses Diamond output and enriches with COG metadata including reference organism.
"""

import sys
import csv
import math
from pathlib import Path


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


def load_protein_to_cog_mapping(cog_csv):
    """Load protein ID (WP_*) to COG ID and assembly ID mappings from cog-24.cog.csv."""
    protein_to_cog = {}
    protein_to_assembly = {}
    print("  Loading protein-to-COG mappings (this may take a moment)...")
    
    with open(cog_csv, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) >= 7:
                assembly_id = parts[1]  # GCF_ assembly ID
                protein_id = parts[2]   # WP_ ID
                cog_id = parts[6]       # COG ID
                protein_to_cog[protein_id] = cog_id
                protein_to_assembly[protein_id] = assembly_id
    
    return protein_to_cog, protein_to_assembly


def load_cog_definitions(def_file):
    """Load COG ID to description mappings."""
    cog_defs = {}
    
    with open(def_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                cog_id = parts[0]
                category = parts[1]
                description = parts[2]
                gene_name = parts[3] if len(parts) > 3 else ''
                pathway = parts[4] if len(parts) > 4 else ''
                cog_defs[cog_id] = {
                    'category': category,
                    'description': description,
                    'gene_name': gene_name,
                    'pathway': pathway
                }
    
    return cog_defs


def load_cog_categories(fun_file):
    """Load COG category codes to full descriptions."""
    categories = {}
    
    with open(fun_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                code = parts[0]
                description = parts[2]
                categories[code] = description
    
    return categories


def load_organism_names(org_csv):
    """Load assembly ID to organism details from cog-24.org.csv.
    Returns dict: assembly_id -> {name, taxid, lineage}
    """
    organisms = {}
    
    with open(org_csv, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) >= 4:
                assembly_id = parts[0]  # GCF_* or GCA_*
                org_name = parts[1].replace('_', ' ')  # Replace underscores with spaces
                taxid = parts[2]
                lineage = parts[3]
                organisms[assembly_id] = {
                    'name': org_name,
                    'taxid': taxid,
                    'lineage': lineage
                }
    
    return organisms


def parse_diamond_output(diamond_file):
    """Parse Diamond output file."""
    hits = {}
    
    with open(diamond_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            protein_id = row['qseqid']
            
            # Keep best hit only
            if protein_id not in hits:
                hits[protein_id] = {
                    'sseqid': row['sseqid'],
                    'pident': float(row['pident']),
                    'length': int(row['length']),
                    'qlen': int(row['qlen']),
                    'slen': int(row['slen']),
                    'evalue': float(row['evalue']),
                    'bitscore': float(row['bitscore']),
                    'qcovhsp': float(row['qcovhsp']),
                    'scovhsp': float(row['scovhsp']),
                    'stitle': row.get('stitle', '')
                }
    
    return hits


def main():
    if len(sys.argv) != 9:
        print("Usage: consolidate_cog.py <diamond_tsv> <input_faa> <cog_csv> <def_tsv> <fun_tsv> <organism> <org_csv> <output_tsv>")
        sys.exit(1)
    
    diamond_file = sys.argv[1]
    fasta_file = sys.argv[2]
    cog_csv = sys.argv[3]
    def_file = sys.argv[4]
    fun_file = sys.argv[5]
    organism = sys.argv[6]
    org_csv = sys.argv[7]
    output_file = sys.argv[8]
    
    print(f"Consolidating COG results for {organism}...")
    
    # Load all data
    print("  Loading protein lengths...")
    proteins = load_protein_lengths(fasta_file)
    print(f"    {len(proteins)} proteins")
    
    # Load COG mappings and metadata
    protein_to_cog, protein_to_assembly = load_protein_to_cog_mapping(cog_csv)
    print(f"    {len(protein_to_cog)} protein-COG mappings")
    
    print("  Loading COG definitions...")
    cog_defs = load_cog_definitions(def_file)
    print(f"    {len(cog_defs)} COG definitions")
    
    print("  Loading COG categories...")
    categories = load_cog_categories(fun_file)
    print(f"    {len(categories)} categories")
    
    print("  Loading organism details...")
    organisms = load_organism_names(org_csv)
    print(f"    {len(organisms)} organisms")
    
    print("  Parsing Diamond output...")
    hits = parse_diamond_output(diamond_file)
    print(f"    {len(hits)} COG hits")
    
    # Write consolidated output
    print(f"  Writing output to {output_file}...")
    
    with open(output_file, 'w') as out:
        # Header - Added COG_reference_genome after COG_id
        header = [
            'organism', 'feature_id', 'protein_length',
            'COG_id', 'COG_reference_genome', 'COG_category_code', 'COG_category_description',
            'COG_description', 'COG_gene_name', 'COG_pathway',
            'COG_pident', 'COG_alignment_length',
            'COG_evalue', 'COG_bitscore', 'COG_qcovhsp', 'COG_scovhsp',
            'COG_qcovhsp_threshold', 'COG_scovhsp_threshold', 'COG_threshold_formula'
        ]
        out.write('\t'.join(header) + '\n')
        
        # Process each protein with hits
        for protein_id in sorted(hits.keys()):
            hit = hits[protein_id]
            protein_length = proteins.get(protein_id, 0)
            
            # Map subject protein ID to COG ID and reference assembly
            subject_id = hit['sseqid']
            cog_id = protein_to_cog.get(subject_id, subject_id)
            assembly_id = protein_to_assembly.get(subject_id, '')
            
            # Format reference genome with all organism details
            if assembly_id and assembly_id in organisms:
                org_info = organisms[assembly_id]
                reference_genome = f"{assembly_id}: {org_info['name']}, {org_info['taxid']}, {org_info['lineage']}"
            else:
                reference_genome = assembly_id
            
            # Get COG metadata
            if cog_id in cog_defs:
                meta = cog_defs[cog_id]
                category_code = meta['category']
                description = meta['description']
                gene_name = meta['gene_name']
                pathway = meta['pathway']
                
                # Get full category descriptions
                category_descs = []
                for char in category_code:
                    if char in categories:
                        category_descs.append(categories[char])
                category_description = '; '.join(category_descs) if category_descs else ''
            else:
                category_code = ''
                category_description = ''
                description = ''
                gene_name = ''
                pathway = ''
            
            # Format row - Added reference_genome after cog_id
            row = [
                organism,
                protein_id,
                protein_length,
                cog_id,
                reference_genome,
                category_code,
                category_description,
                description,
                gene_name,
                pathway,
                f"{hit['pident']:.1f}",
                hit['length'],
                f"{hit['evalue']:.2e}" if hit['evalue'] > 0 else "0.0",
                f"{hit['bitscore']:.1f}",
                f"{hit['qcovhsp']:.2f}",
                f"{hit['scovhsp']:.2f}",
                "0",  # COG_qcovhsp_threshold: No filter (COGclassifier approach)
                "0",  # COG_scovhsp_threshold: No filter (COGclassifier approach)
                "E-value<=0.01, best-hit only (COGclassifier approach)"  # Formula
            ]
            
            out.write('\t'.join(str(x) for x in row) + '\n')
    
    print(f"\n  Summary:")
    print(f"    Total proteins: {len(proteins)}")
    print(f"    Proteins with COG hits: {len(hits)} ({100*len(hits)/len(proteins):.1f}%)")
    print(f"\n✓ Consolidation complete!")


if __name__ == '__main__':
    main()
