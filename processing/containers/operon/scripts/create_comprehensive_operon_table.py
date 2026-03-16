#!/usr/bin/env python3
"""
Create comprehensive gene-centric operon annotation table.

This script parses UniOP outputs and creates a detailed table where each row represents
a gene with full operon context, prediction scores, and genomic distances.

Output columns:
- feature_id: Gene identifier
- operon_id: Operon identifier (or NA if singleton)
- operon_size: Number of genes in operon
- operon_position: Position within operon (1-based)
- member_genes: Comma-separated list of all genes in operon
- upstream_gene: Adjacent gene upstream (by genomic position)
- downstream_gene: Adjacent gene downstream (by genomic position)
- upstream_operon_score: UniOP prediction score with upstream gene
- downstream_operon_score: UniOP prediction score with downstream gene
- intergenic_distance_upstream: Distance to upstream gene (bp, negative if overlapping)
- intergenic_distance_downstream: Distance to downstream gene (bp, negative if overlapping)
- contig: Contig/chromosome name
- start: Gene start coordinate
- end: Gene end coordinate
- strand: Strand (+/-)
"""

import sys
import re
from pathlib import Path

def parse_prodigal_faa(faa_file):
    """Parse Prodigal-style FAA to extract gene coordinates and mapping."""
    gene_info = {}
    gene_order = []
    
    with open(faa_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # >NODE_10_length_104094_cov_126.636434_100 # 3 # 644 # -1 # ID=fig|6666666.1481409.peg.1211
                parts = line[1:].strip().split(' # ')
                if len(parts) >= 5:
                    contig_with_num = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    strand_val = int(parts[3])
                    strand = '+' if strand_val == 1 else '-'
                    
                    # Extract gene ID from last part
                    id_part = parts[4]
                    match = re.search(r'ID=(.+?)(?:\s|$)', id_part)
                    if match:
                        gene_id = match.group(1)
                        
                        # Extract contig name (remove _genenum suffix)
                        contig = '_'.join(contig_with_num.split('_')[:-1])
                        
                        gene_info[gene_id] = {
                            'contig': contig,
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'contig_with_num': contig_with_num
                        }
                        gene_order.append(gene_id)
    
    return gene_info, gene_order

def parse_prediction_scores(pred_file, gene_order):
    """Parse UniOP prediction scores between adjacent genes."""
    scores = {}
    
    with open(pred_file, 'r') as f:
        # Skip header lines
        for line in f:
            if line.startswith('Gene A'):
                break
        
        # Read predictions
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                try:
                    gene_a_idx = int(parts[0]) - 1  # Convert to 0-based
                    gene_b_idx = int(parts[1]) - 1
                    score = float(parts[2]) if parts[2].strip() else None
                    
                    if score is not None and gene_a_idx < len(gene_order) and gene_b_idx < len(gene_order):
                        gene_a = gene_order[gene_a_idx]
                        gene_b = gene_order[gene_b_idx]
                        scores[(gene_a, gene_b)] = score
                except (ValueError, IndexError):
                    continue
    
    return scores

def parse_operon_file(operon_file, gene_order):
    """Parse reformatted operon assignments (processed_uniop_*.operon format)."""
    operons = {}
    
    with open(operon_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('Organism'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                operon_id = parts[1]
                gene_ids_str = parts[3]
                
                # Parse gene IDs (comma-separated with spaces)
                gene_ids = [gid.strip() for gid in gene_ids_str.split(',')]
                
                operons[operon_id] = gene_ids
    
    return operons

def calculate_intergenic_distance(gene1_info, gene2_info):
    """Calculate intergenic distance between two genes."""
    # Only calculate if on same contig
    if gene1_info['contig'] != gene2_info['contig']:
        return None
    
    # Get boundaries
    gene1_left = min(gene1_info['start'], gene1_info['end'])
    gene1_right = max(gene1_info['start'], gene1_info['end'])
    gene2_left = min(gene2_info['start'], gene2_info['end'])
    gene2_right = max(gene2_info['start'], gene2_info['end'])
    
    # Calculate distance (negative if overlapping)
    if gene1_right < gene2_left:
        return gene2_left - gene1_right - 1
    elif gene2_right < gene1_left:
        return gene1_left - gene2_right - 1
    else:
        # Overlapping
        overlap = min(gene1_right, gene2_right) - max(gene1_left, gene2_left) + 1
        return -overlap

def find_adjacent_genes(gene_id, gene_info, gene_order):
    """Find upstream and downstream adjacent genes on same contig."""
    current_info = gene_info[gene_id]
    current_contig = current_info['contig']
    current_start = current_info['start']
    current_end = current_info['end']
    current_pos = (current_start + current_end) // 2
    
    # Find genes on same contig
    contig_genes = [(gid, gene_info[gid]) for gid in gene_order if gene_info[gid]['contig'] == current_contig]
    
    # Sort by position
    contig_genes.sort(key=lambda x: (x[1]['start'] + x[1]['end']) // 2)
    
    # Find current gene index
    gene_positions = [((g[1]['start'] + g[1]['end']) // 2, g[0]) for g in contig_genes]
    current_idx = next((i for i, (pos, gid) in enumerate(gene_positions) if gid == gene_id), None)
    
    if current_idx is None:
        return None, None
    
    upstream = gene_positions[current_idx - 1][1] if current_idx > 0 else None
    downstream = gene_positions[current_idx + 1][1] if current_idx < len(gene_positions) - 1 else None
    
    return upstream, downstream

def create_comprehensive_table(genome_name, prodigal_faa, pred_file, processed_operon_file, output_file):
    """Create comprehensive operon annotation table."""
    
    print(f"Creating comprehensive operon table for: {genome_name}")
    print(f"  Input FAA: {prodigal_faa}")
    print(f"  Input predictions: {pred_file}")
    print(f"  Input operons: {processed_operon_file}")
    print()
    
    # Parse input files
    print("Parsing gene coordinates...")
    gene_info, gene_order = parse_prodigal_faa(prodigal_faa)
    print(f"  Found {len(gene_info)} genes")
    
    print("Parsing prediction scores...")
    scores = parse_prediction_scores(pred_file, gene_order)
    print(f"  Found {len(scores)} pairwise predictions")
    
    print("Parsing operon assignments...")
    operons = parse_operon_file(processed_operon_file, gene_order)
    print(f"  Found {len(operons)} operons")
    
    # Create gene-to-operon mapping
    gene_to_operon = {}
    for operon_id, members in operons.items():
        for i, gene_id in enumerate(members, start=1):
            gene_to_operon[gene_id] = {
                'operon_id': operon_id,
                'operon_size': len(members),
                'operon_position': i,
                'member_genes': ','.join(members)
            }
    
    # Build comprehensive table
    print("\nBuilding comprehensive table...")
    rows = []
    
    for gene_id in gene_order:
        info = gene_info[gene_id]
        
        # Operon membership
        if gene_id in gene_to_operon:
            operon_data = gene_to_operon[gene_id]
            operon_id = operon_data['operon_id']
            operon_size = operon_data['operon_size']
            operon_position = operon_data['operon_position']
            member_genes = operon_data['member_genes']
        else:
            operon_id = 'NA'
            operon_size = 1
            operon_position = 1
            member_genes = gene_id
        
        # Find adjacent genes
        upstream_gene, downstream_gene = find_adjacent_genes(gene_id, gene_info, gene_order)
        
        # Get prediction scores
        upstream_score = scores.get((upstream_gene, gene_id), None) if upstream_gene else None
        downstream_score = scores.get((gene_id, downstream_gene), None) if downstream_gene else None
        
        # Calculate intergenic distances
        upstream_dist = None
        downstream_dist = None
        
        if upstream_gene:
            upstream_dist = calculate_intergenic_distance(gene_info[upstream_gene], info)
        
        if downstream_gene:
            downstream_dist = calculate_intergenic_distance(info, gene_info[downstream_gene])
        
        # Format values
        upstream_gene_str = upstream_gene if upstream_gene else 'NA'
        downstream_gene_str = downstream_gene if downstream_gene else 'NA'
        upstream_score_str = f"{upstream_score:.6f}" if upstream_score is not None else 'NA'
        downstream_score_str = f"{downstream_score:.6f}" if downstream_score is not None else 'NA'
        upstream_dist_str = str(upstream_dist) if upstream_dist is not None else 'NA'
        downstream_dist_str = str(downstream_dist) if downstream_dist is not None else 'NA'
        
        row = {
            'feature_id': gene_id,
            'operon_id': operon_id,
            'operon_size': operon_size,
            'operon_position': operon_position,
            'member_genes': member_genes,
            'upstream_gene': upstream_gene_str,
            'downstream_gene': downstream_gene_str,
            'upstream_operon_score': upstream_score_str,
            'downstream_operon_score': downstream_score_str,
            'intergenic_distance_upstream': upstream_dist_str,
            'intergenic_distance_downstream': downstream_dist_str,
            'contig': info['contig'],
            'start': info['start'],
            'end': info['end'],
            'strand': info['strand']
        }
        
        rows.append(row)
    
    # Write output
    print(f"\nWriting output: {output_file}")
    with open(output_file, 'w') as f:
        # Write header
        header = [
            'feature_id', 'operon_id', 'operon_size', 'operon_position', 'member_genes',
            'upstream_gene', 'downstream_gene', 
            'upstream_operon_score', 'downstream_operon_score',
            'intergenic_distance_upstream', 'intergenic_distance_downstream',
            'contig', 'start', 'end', 'strand'
        ]
        f.write('\t'.join(header) + '\n')
        
        # Write rows
        for row in rows:
            values = [str(row[col]) for col in header]
            f.write('\t'.join(values) + '\n')
    
    # Statistics
    genes_in_operons = sum(1 for row in rows if row['operon_id'] != 'NA')
    total_genes = len(rows)
    coverage = 100.0 * genes_in_operons / total_genes if total_genes > 0 else 0
    
    print(f"\n✓ Created comprehensive operon table")
    print(f"  Total genes: {total_genes}")
    print(f"  Genes in operons: {genes_in_operons}")
    print(f"  Coverage: {coverage:.1f}%")
    print(f"  Operons: {len(operons)}")

if __name__ == '__main__':
    if len(sys.argv) != 6:
        print("Usage: create_comprehensive_operon_table.py <genome_name> <prodigal_faa> <pred_file> <processed_operon_file> <output_file>")
        sys.exit(1)
    
    genome_name = sys.argv[1]
    prodigal_faa = sys.argv[2]
    pred_file = sys.argv[3]
    processed_operon_file = sys.argv[4]
    output_file = sys.argv[5]
    
    create_comprehensive_table(genome_name, prodigal_faa, pred_file, processed_operon_file, output_file)
