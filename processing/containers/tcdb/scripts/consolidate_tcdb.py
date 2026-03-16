#!/usr/bin/env python3
"""
Consolidate TCDB annotation output into comprehensive TSV format.
Parses Diamond output and extracts TCDB metadata (TC number, family, substrates).
Includes confidence scoring components.
"""

import sys
import csv
import math
import re
from datetime import datetime


def log(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}]  {message}", flush=True)


def load_protein_lengths(fasta_file):
    """
    Calculate protein lengths from input FASTA file.
    Returns dict: protein_id -> length
    """
    proteins = {}
    current_id = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_id:
                    proteins[current_id] = len(''.join(current_seq))
                # Start new sequence
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_id:
            proteins[current_id] = len(''.join(current_seq))
    
    return proteins


def load_tcdb_families(families_file):
    """
    Load TCDB family descriptions.
    Format: TC_number \t family_description
    Returns dict: TC_family -> description
    """
    families = {}
    
    with open(families_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t', 1)
            if len(parts) == 2:
                tc_family = parts[0].strip()
                description = parts[1].strip()
                families[tc_family] = description
    
    return families


def load_tcdb_substrates(substrates_file):
    """
    Load TCDB substrate information.
    Format: TC_id \t CHEBI:id;name|CHEBI:id;name
    Returns dict: TC_id -> substrates_string
    """
    substrates = {}
    
    with open(substrates_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t', 1)
            if len(parts) == 2:
                tc_id = parts[0].strip()
                substrate_info = parts[1].strip()
                substrates[tc_id] = substrate_info
    
    return substrates


def parse_tcdb_stitle(stitle):
    """
    Parse TCDB stitle to extract metadata.
    Format: gnl|TC-DB|accession|TC_number description
    
    Example: gnl|TC-DB|1CN3_F|1.A.83.1.5 Chain F, FRAGMENT OF COAT PROTEIN VP2
    """
    result = {
        'accession': '',
        'tc_id': '',
        'description': ''
    }
    
    if not stitle:
        return result
    
    # Extract accession and TC number
    # Format: gnl|TC-DB|accession|TC_number description
    match = re.match(r'gnl\|TC-DB\|([^|]+)\|([0-9.A-Z]+)\s+(.*)', stitle)
    if match:
        result['accession'] = match.group(1).strip()
        result['tc_id'] = match.group(2).strip()
        result['description'] = match.group(3).strip()
    
    return result


def get_tc_family_from_id(tc_id):
    """
    Extract family from TC ID.
    Example: 1.A.83.1.5 -> 1.A.83 (first three levels)
    """
    parts = tc_id.split('.')
    if len(parts) >= 3:
        return '.'.join(parts[:3])
    return tc_id


def calculate_confidence(evalue, qstart, qend, protein_length, pident):
    """
    Calculate confidence score components.
    
    Returns:
        evalue_score: Statistical significance score [0-1]
        coverage_score: Query coverage score [0-1]
        identity_score: Sequence identity score [0-1]
        confidence: Weighted final score [0-1]
    """
    # E-value score (statistical significance)
    if evalue == 0.0:
        evalue_score = 1.0
    else:
        evalue_score = min(1.0, max(0.0, -math.log10(evalue) / 50.0))
    
    # Coverage score
    query_coverage = (qend - qstart + 1) / protein_length if protein_length > 0 else 0
    coverage_score = min(1.0, max(0.0, (query_coverage - 0.3) / 0.7))
    
    # Identity score
    identity_score = min(1.0, max(0.0, pident / 100.0))
    
    # Final confidence (weighted sum: 60% evalue, 30% coverage, 10% identity)
    confidence = 0.6 * evalue_score + 0.3 * coverage_score + 0.1 * identity_score
    
    return evalue_score, coverage_score, identity_score, confidence


def parse_tcdb_hits_output(tcdb_hits_file):
    """
    Parse Diamond output file.
    Returns dict: protein_id -> hit_info
    """
    hits = {}
    
    with open(tcdb_hits_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            protein_id = row['qseqid']
            
            # Keep best hit only (max-target-seqs 1, but just in case)
            if protein_id not in hits:
                hits[protein_id] = {
                    'sseqid': row['sseqid'],
                    'pident': float(row['pident']),
                    'length': int(row['length']),
                    'qstart': int(row['qstart']),
                    'qend': int(row['qend']),
                    'evalue': float(row['evalue']),
                    'bitscore': float(row['bitscore']),
                    'qlen': int(row['qlen']),
                    'slen': int(row['slen']),
                    'qcovhsp': float(row['qcovhsp']),
                    'scovhsp': float(row['scovhsp']),
                    'stitle': row['stitle']
                }
    
    return hits


def main():
    if len(sys.argv) != 7:
        log("Usage: consolidate_tcdb.py <tcdb_hits_tsv> <input_faa> <families_tsv> <substrates_tsv> <organism> <output_tsv>")
        sys.exit(1)
    
    tcdb_hits_file = sys.argv[1]
    fasta_file = sys.argv[2]
    families_file = sys.argv[3]
    substrates_file = sys.argv[4]
    organism = sys.argv[5]
    output_file = sys.argv[6]
    
    log("="*70)
    log("TCDB Annotation Consolidation")
    log("="*70)
    log("")
    log(f"Organism:        {organism}")
    log(f"Input files:")
    log(f"  TCDB hits TSV: {tcdb_hits_file}")
    log(f"  Protein FASTA: {fasta_file}")
    log(f"  Families:      {families_file}")
    log(f"  Substrates:    {substrates_file}")
    log(f"Output file:     {output_file}")
    log("")
    
    # Load all data
    log("Loading reference data...")
    log("")
    
    log("  Loading protein lengths...")
    proteins = load_protein_lengths(fasta_file)
    log(f"    Loaded {len(proteins)} proteins")
    
    log("  Loading TCDB families...")
    families = load_tcdb_families(families_file)
    log(f"    Loaded {len(families)} families")
    
    log("  Loading TCDB substrates...")
    substrates = load_tcdb_substrates(substrates_file)
    log(f"    Loaded {len(substrates)} substrate annotations")
    log("")
    
    log("Parsing TCDB hits output...")
    hits = parse_tcdb_hits_output(tcdb_hits_file)
    log(f"  Found {len(hits)} TCDB hits")
    log("")
    
    # Write consolidated output
    log(f"Writing consolidated output...")
    log(f"  Output file: {output_file}")
    log("")
    
    with open(output_file, 'w') as out:
        # Header
        header = [
            'organism', 'feature_id', 'protein_length',
            'TCDB_id', 'TCDB_accession', 'TCDB_description', 'TCDB_family',
            'TCDB_family_description', 'TCDB_substrates',
            'TCDB_pident', 'TCDB_alignment_length', 'TCDB_qstart', 'TCDB_qend',
            'TCDB_evalue', 'TCDB_bitscore', 'TCDB_qcovhsp', 'TCDB_scovhsp',
            'TCDB_qcovhsp_threshold', 'TCDB_scovhsp_threshold', 'TCDB_threshold_formula'
        ]
        out.write('\t'.join(header) + '\n')
        
        # Process each protein (only those with hits)
        for protein_id in sorted(hits.keys()):
            hit = hits[protein_id]
            protein_length = proteins.get(protein_id, 0)
            
            # Parse stitle
            stitle_info = parse_tcdb_stitle(hit['stitle'])
            tc_id = stitle_info['tc_id']
            accession = stitle_info['accession']
            description = stitle_info['description']
            
            # Get family and family description
            tc_family = get_tc_family_from_id(tc_id)
            family_description = families.get(tc_family, '')
            
            # Get substrates
            substrate_info = substrates.get(tc_id, '')
            
            # Format row
            row = [
                organism,
                protein_id,
                protein_length,
                tc_id,
                accession,
                description,
                tc_family,
                family_description,
                substrate_info,
                f"{hit['pident']:.1f}",
                hit['length'],
                hit['qstart'],
                hit['qend'],
                f"{hit['evalue']:.2e}" if hit['evalue'] > 0 else "0.0",
                f"{hit['bitscore']:.1f}",
                f"{hit['qcovhsp']:.2f}",
                f"{hit['scovhsp']:.2f}",
                "50.0",  # TCDB_qcovhsp_threshold: fixed 50% filter
                "50.0",  # TCDB_scovhsp_threshold: fixed 50% filter
                "qcovhsp>=50% AND scovhsp>=50%"  # Formula from TCDB search runner
            ]
            
            out.write('\t'.join(str(x) for x in row) + '\n')
    
    log("Output file structure:")
    log("  Columns: 20 (organism, feature_id, protein_length, TCDB annotations)")
    log("")
    log("Results summary:")
    log(f"  Total proteins:             {len(proteins)}")
    log(f"  Proteins with TCDB hits:    {len(hits)}")
    if len(proteins) > 0:
        coverage = 100 * len(hits) / len(proteins)
        log(f"  Coverage:                   {coverage:.1f}%")
    log("")
    log("✓ TCDB consolidation complete!")


if __name__ == '__main__':
    main()
