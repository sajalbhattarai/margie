#!/usr/bin/env python3
"""
Reformat UniOP output to use original gene IDs and add organism name.

Usage:
    python reformat_uniop_output.py <organism_name> <input.faa> <uniop.pred> <uniop.operon> [output_dir]

Example:
    python reformat_uniop_output.py "Bacteroides thetaiotaomicron" theta.faa uniop.pred uniop.operon output/
"""

import sys
import re
from pathlib import Path


def parse_faa_gene_ids(faa_file):
    """
    Extract gene IDs from Prodigal-style FASTA file in sequential order.
    Returns dict: {1: 'fig|6666666.1476522.peg.XXX', 2: 'fig|6666666.1476522.peg.YYY', ...}
    
    Keys are SEQUENTIAL POSITIONS in the FASTA file (1-indexed), as UniOP uses these positions.
    
    Parses headers like:
    >NODE_10_length_197040_cov_20.399664_152 # 1095 # 2867 # -1 # ID=fig|6666666.1476522.peg.2281
    
    Extracts the original gene ID from the ID= field.
    """
    gene_ids = {}
    seq_position = 0  # Sequential position counter (will be 1-indexed)
    
    with open(faa_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                seq_position += 1  # Increment for each gene (1-indexed)
                
                # Parse Prodigal-style header
                # Format: >contig_name_genenum # start # stop # strand # ID=original_gene_id
                parts = line.strip().split(' # ')
                if len(parts) >= 5:
                    # Extract original gene ID from ID= field
                    id_field = parts[4].strip()
                    if id_field.startswith('ID='):
                        original_id = id_field[3:]  # Remove 'ID='
                        gene_ids[seq_position] = original_id
                    else:
                        # Fallback: use contig name
                        contig_gene = parts[0][1:]  # Remove '>'
                        gene_ids[seq_position] = contig_gene
                else:
                    # Fallback: use first word as gene ID
                    gene_id = line.strip().split()[0][1:]
                    gene_ids[seq_position] = gene_id
    
    return gene_ids


def reformat_pred_file(pred_file, gene_ids, organism, output_file):
    """
    Reformat uniop.pred to use actual gene IDs and add organism name.
    """
    with open(pred_file, 'r') as f_in, open(output_file, 'w') as f_out:
        lines = f_in.readlines()
        
        # Write header with organism column
        f_out.write("# Reformatted operon prediction scores\n")
        f_out.write("Organism\tGene_A\tGene_B\tPrediction_Confidence\n")
        
        # Process data lines (skip first 2 header lines)
        for line in lines[2:]:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                gene_a = parts[0].strip()
                gene_b = parts[1].strip()
                prediction = parts[2].strip() if len(parts) > 2 else ''
                
                # Convert to gene IDs
                if gene_a and gene_a.isdigit():
                    gene_a_id = gene_ids.get(int(gene_a), gene_a)
                else:
                    gene_a_id = gene_a
                    
                if gene_b and gene_b.isdigit():
                    gene_b_id = gene_ids.get(int(gene_b), gene_b)
                else:
                    gene_b_id = gene_b
                
                # Add organism name as first column
                f_out.write(f"{organism}\t{gene_a_id}\t{gene_b_id}\t{prediction}\n")
    
    print(f"✓ Created: {output_file}")


def reformat_operon_file(operon_file, gene_ids, organism, output_file):
    """
    Reformat uniop.operon to use actual gene IDs and add organism name.
    """
    with open(operon_file, 'r') as f_in, open(output_file, 'w') as f_out:
        lines = f_in.readlines()
        
        # Write header with organism column
        f_out.write("# Reformatted operon clusters\n")
        f_out.write("Organism\tOperon_ID\tNum_Genes\tGene_IDs\n")
        
        # Process data lines (skip header)
        for line in lines[1:]:
            # Parse the CSV line: "[np.int64(2), np.int64(3), np.int64(4)]",0
            # Split by last comma to separate genes from operon_id
            parts = line.strip().rsplit(',', 1)
            if len(parts) == 2:
                genes_part = parts[0]
                operon_id = parts[1].strip()
                
                # Extract all numbers from np.int64(...) patterns
                gene_nums = re.findall(r'np\.int64\((\d+)\)', genes_part)
                
                # Convert to gene IDs
                gene_id_list = []
                for gene_num in gene_nums:
                    gene_id = gene_ids.get(int(gene_num), gene_num)
                    gene_id_list.append(gene_id)
                
                # Format output
                num_genes = len(gene_id_list)
                gene_ids_str = ', '.join(gene_id_list)
                
                f_out.write(f"{organism}\tOperon_{operon_id}\t{num_genes}\t{gene_ids_str}\n")
    
    print(f"✓ Created: {output_file}")


def main():
    if len(sys.argv) < 5:
        print("Usage: python reformat_uniop_output.py <organism_name> <input.faa> <uniop.pred> <uniop.operon> [output_dir]")
        sys.exit(1)
    
    organism = sys.argv[1]
    faa_file = sys.argv[2]
    pred_file = sys.argv[3]
    operon_file = sys.argv[4]
    output_dir = sys.argv[5] if len(sys.argv) > 5 else '.'
    
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Extract organism short name from full name (e.g., "theta" from "Bacteroides thetaiotaomicron")
    organism_short = organism.split()[-1].lower()
    
    # Output files
    pred_output = f"{output_dir}/processed_uniop_{organism_short}.pred"
    operon_output = f"{output_dir}/processed_uniop_{organism_short}.operon"
    
    print(f"Reformatting UniOP output for: {organism}")
    print(f"  Input FAA: {faa_file}")
    print(f"  Input pred: {pred_file}")
    print(f"  Input operon: {operon_file}")
    print("")
    
    # Parse gene IDs from FASTA
    print("Parsing gene IDs from FASTA...")
    gene_ids = parse_faa_gene_ids(faa_file)
    print(f"  Found {len(gene_ids)} genes")
    print("")
    
    # Reformat pred file
    print("Reformatting prediction file...")
    reformat_pred_file(pred_file, gene_ids, organism, pred_output)
    print("")
    
    # Reformat operon file
    print("Reformatting operon file...")
    reformat_operon_file(operon_file, gene_ids, organism, operon_output)
    print("")
    
    print("✓ Reformatting complete!")


if __name__ == "__main__":
    main()
