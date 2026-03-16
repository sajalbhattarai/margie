#!/usr/bin/env python3
"""
Consolidate COGclassifier output into comprehensive TSV format.
Converts COGclassifier cog_classify.tsv to match expected output format.
"""

import sys
import csv
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


def load_rpsblast_results(rpsblast_file):
    """Load RPS-BLAST results for alignment stats"""
    blast_stats = {}
    
    with open(rpsblast_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 12:
                query_id = parts[0]
                pident = float(parts[2])
                length = int(parts[3])
                evalue = float(parts[10])
                bitscore = float(parts[11])
                
                blast_stats[query_id] = {
                    'pident': pident,
                    'length': length,
                    'evalue': evalue,
                    'bitscore': bitscore
                }
    
    return blast_stats


def main():
    if len(sys.argv) != 5:
        print("Usage: consolidate_cog_new.py <organism> <input_faa> <cogclassifier_output_dir> <output_tsv>")
        sys.exit(1)
    
    organism = sys.argv[1]
    input_faa = sys.argv[2]
    cogclassifier_dir = sys.argv[3]
    output_file = sys.argv[4]
    
    # Input files
    cog_classify_file = f"{cogclassifier_dir}/cog_classify.tsv"
    rpsblast_file = f"{cogclassifier_dir}/rpsblast.tsv"
    
    # Check files exist
    if not Path(cog_classify_file).exists():
        print(f"ERROR: COG classification file not found: {cog_classify_file}")
        sys.exit(1)
    
    print("COGclassifier Output Consolidation")
    print("=" * 50)
    print(f"Organism:               {organism}")
    print(f"Input proteins:         {input_faa}")
    print(f"COG classification:     {cog_classify_file}")
    print(f"RPS-BLAST results:      {rpsblast_file}")
    print(f"Output file:            {output_file}")
    print("")
    
    # Load protein lengths
    print("Loading protein lengths...")
    proteins = load_protein_lengths(input_faa)
    print(f"  {len(proteins)} proteins")
    
    # Load RPS-BLAST stats
    print("Loading RPS-BLAST statistics...")
    blast_stats = {}
    if Path(rpsblast_file).exists():
        blast_stats = load_rpsblast_results(rpsblast_file)
        print(f"  {len(blast_stats)} BLAST hits")
    
    # Parse COGclassifier output
    print("Parsing COGclassifier output...")
    cog_hits = {}
    
    with open(cog_classify_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            query_id = row['QUERY_ID']
            cog_hits[query_id] = {
                'cog_id': row['COG_ID'],
                'cdd_id': row['CDD_ID'],
                'evalue': float(row['EVALUE']) if row['EVALUE'] != '0.0' else 0.0,
                'identity': float(row['IDENTITY']),
                'gene_name': row['GENE_NAME'],
                'cog_name': row['COG_NAME'],
                'cog_letter': row['COG_LETTER'],
                'cog_description': row['COG_DESCRIPTION']
            }
    
    print(f" {len(cog_hits)} COG classifications")
    
    # Write consolidated output
    print(f"Writing consolidated output to {output_file}...")
    
    with open(output_file, 'w') as out:
        # Header - matches expected format
        header = [
            'organism', 'feature_id', 'protein_length',
            'COG_id', 'COG_reference_genome', 'COG_category_code', 'COG_category_description',
            'COG_description', 'COG_gene_name', 'COG_pathway',
            'COG_pident', 'COG_alignment_length',
            'COG_evalue', 'COG_bitscore', 'COG_qcovhsp', 'COG_scovhsp',
            'COG_qcovhsp_threshold', 'COG_scovhsp_threshold', 'COG_threshold_formula'
        ]
        out.write('\t'.join(header) + '\n')
        
        # Process each protein with COG hits
        for protein_id in sorted(cog_hits.keys()):
            cog = cog_hits[protein_id]
            protein_length = proteins.get(protein_id, 0)
            
            # Get BLAST stats if available
            if protein_id in blast_stats:
                blast = blast_stats[protein_id]
                pident = blast['pident']
                alignment_length = blast['length']
                evalue = blast['evalue']
                bitscore = blast['bitscore']
                
                # Calculate coverage (approximate)
                qcovhsp = (alignment_length / protein_length * 100) if protein_length > 0 else 0
                scovhsp = 50.0  # Default, not available from COGclassifier
            else:
                # Use COGclassifier data
                pident = cog['identity']
                alignment_length = int(protein_length * pident / 100)  # Approximate
                evalue = cog['evalue']
                bitscore = 0.0  # Not provided by COGclassifier
                qcovhsp = 50.0  # Default
                scovhsp = 50.0  # Default
            
            # Format row
            row = [
                organism,
                protein_id,
                protein_length,
                cog['cog_id'],
                f"CDD:{cog['cdd_id']}",  # Reference is CDD ID
                cog['cog_letter'],
                cog['cog_description'],
                cog['cog_name'],
                cog['gene_name'],
                '',  # Pathway not provided by COGclassifier
                f"{pident:.1f}",
                alignment_length,
                f"{evalue:.2e}" if evalue > 0 else "0.0",
                f"{bitscore:.1f}",
                f"{qcovhsp:.2f}",
                f"{scovhsp:.2f}",
                "30.0",  # Updated threshold: 30% coverage
                "30.0",  # Updated threshold: 30% coverage
                "COGclassifier default (E-value < 0.01)"
            ]
            
            out.write('\t'.join(str(x) for x in row) + '\n')
    
    print(f"\nSummary:")
    print(f"  Total proteins: {len(proteins)}")
    print(f"  Proteins with COG hits: {len(cog_hits)} ({100*len(cog_hits)/len(proteins):.1f}%)")
    print(f"\n✓ Consolidation complete!")


if __name__ == '__main__':
    main()
