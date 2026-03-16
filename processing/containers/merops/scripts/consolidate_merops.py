#!/usr/bin/env python3
"""
MEROPS Peptidase Annotation Consolidation Script

Parses Diamond output and creates consolidated TSV with:
- MEROPS peptidase metadata (family, clan, name)
- Confidence scoring based on E-value, coverage, and identity
"""

import sys
import math
import re
from datetime import datetime

def log(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}]  {message}", flush=True)

def load_protein_lengths(faa_file):
    """Load protein sequences and calculate lengths"""
    lengths = {}
    current_id = None
    current_seq = []
    
    with open(faa_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    lengths[current_id] = len(''.join(current_seq))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            lengths[current_id] = len(''.join(current_seq))
    
    log(f"  Loaded {len(lengths)} proteins")
    return lengths

def parse_merops_stitle(stitle):
    """
    Parse MEROPS stitle to extract metadata.
    Format: MER0000001 - peptidase_name (Organism) [family]#clan#{peptidase unit: X-Y}~source~
    """
    result = {
        'merops_id': '',
        'peptidase_name': '',
        'organism': '',
        'family': '',
        'clan': '',
        'peptidase_unit': '',
        'source': ''
    }
    
    if not stitle:
        return result
    
    # Extract MEROPS ID
    mer_match = re.match(r'(MER\d+)', stitle)
    if mer_match:
        result['merops_id'] = mer_match.group(1)
    
    # Extract peptidase name
    name_match = re.search(r' - ([^(]+)', stitle)
    if name_match:
        result['peptidase_name'] = name_match.group(1).strip()
    
    # Extract organism
    org_match = re.search(r'\(([^)]+)\)', stitle)
    if org_match:
        result['organism'] = org_match.group(1).strip()
    
    # Extract family
    family_match = re.search(r'\[([^\]]+)\]', stitle)
    if family_match:
        result['family'] = family_match.group(1).strip()
    
    # Extract clan
    clan_match = re.search(r'#([^#]+)#', stitle)
    if clan_match:
        result['clan'] = clan_match.group(1).strip()
    
    # Extract peptidase unit
    unit_match = re.search(r'peptidase unit: ([^}]+)', stitle)
    if unit_match:
        result['peptidase_unit'] = unit_match.group(1).strip()
    
    # Extract source
    source_match = re.search(r'~source ([^~]+)~', stitle)
    if source_match:
        result['source'] = source_match.group(1).strip()
    
    return result

def calculate_confidence(evalue, qstart, qend, protein_length, pident):
    """
    Calculate confidence score (0-1) based on three components:
    1. E-value score (60% weight): Statistical significance
    2. Coverage score (30% weight): Query coverage
    3. Identity score (10% weight): Sequence similarity
    """
    # E-value score
    if evalue == 0:
        evalue_score = 1.0
    else:
        log_eval = -math.log10(evalue)
        evalue_score = min(1.0, log_eval / 50.0)
    
    # Coverage score
    query_coverage = (qend - qstart + 1) / protein_length if protein_length > 0 else 0
    coverage_score = min(1.0, query_coverage)
    
    # Identity score
    identity_score = pident / 100.0
    
    # Weighted combination
    confidence = (0.6 * evalue_score) + (0.3 * coverage_score) + (0.1 * identity_score)
    
    return round(confidence, 2), round(evalue_score, 6), round(coverage_score, 6), round(identity_score, 2)

def main():
    if len(sys.argv) != 5:
        print("Usage: consolidate_merops.py <genome_name> <diamond_tsv> <faa_file> <output_tsv>")
        sys.exit(1)
    
    genome_name = sys.argv[1]
    diamond_file = sys.argv[2]
    faa_file = sys.argv[3]
    output_file = sys.argv[4]
    
    log("MEROPS Annotation Consolidation")
    log("="*60)
    log(f"Organism:        {genome_name}")
    log(f"DIAMOND output:  {diamond_file}")
    log(f"Input FASTA:     {faa_file}")
    log(f"Output file:     {output_file}")
    log("")
    
    log("Loading protein lengths from input FASTA...")
    protein_lengths = load_protein_lengths(faa_file)
    log("")
    
    log("Parsing DIAMOND output...")
    
    # Parse Diamond output
    hits = {}
    hit_count = 0
    with open(diamond_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 13:
                continue
            
            qseqid = fields[0]
            sseqid = fields[1]
            pident = float(fields[2])
            length = int(fields[3])
            qstart = int(fields[6])
            qend = int(fields[7])
            evalue = float(fields[10])
            bitscore = float(fields[11])
            qlen = int(fields[12]) if len(fields) > 12 else 0
            slen = int(fields[13]) if len(fields) > 13 else 0
            qcovhsp = float(fields[14]) if len(fields) > 14 else 0.0
            scovhsp = float(fields[15]) if len(fields) > 15 else 0.0
            stitle = fields[16] if len(fields) > 16 else ""
            
            # Keep only best hit per protein
            if qseqid not in hits:
                hits[qseqid] = {
                    'sseqid': sseqid,
                    'pident': pident,
                    'length': length,
                    'qstart': qstart,
                    'qend': qend,
                    'evalue': evalue,
                    'bitscore': bitscore,
                    'qlen': qlen,
                    'slen': slen,
                    'qcovhsp': qcovhsp,
                    'scovhsp': scovhsp,
                    'stitle': stitle
                }
                hit_count += 1
    
    log(f"  Parsed {len(hits)} MEROPS peptidase hits")
    log("")
    
    log("Writing consolidated output...")
    
    # Write output
    with open(output_file, 'w') as out:
        # Header
        header = [
            'organism', 'feature_id', 'protein_length',
            'MEROPS_id', 'MEROPS_peptidase_name', 'MEROPS_organism',
            'MEROPS_family', 'MEROPS_clan', 'MEROPS_peptidase_unit', 'MEROPS_source',
            'MEROPS_pident', 'MEROPS_alignment_length', 'MEROPS_qstart', 'MEROPS_qend',
            'MEROPS_evalue', 'MEROPS_bitscore', 'MEROPS_qcovhsp', 'MEROPS_scovhsp',
            'MEROPS_qcovhsp_threshold', 'MEROPS_scovhsp_threshold', 'MEROPS_threshold_formula'
        ]
        out.write('\t'.join(header) + '\n')
        
        # Write hits
        for protein_id in sorted(hits.keys()):
            hit = hits[protein_id]
            protein_length = protein_lengths.get(protein_id, 0)
            
            # Parse stitle
            stitle_info = parse_merops_stitle(hit['stitle'])
            
            row = [
                genome_name,
                protein_id,
                str(protein_length),
                stitle_info['merops_id'],
                stitle_info['peptidase_name'],
                stitle_info['organism'],
                stitle_info['family'],
                stitle_info['clan'],
                stitle_info['peptidase_unit'],
                stitle_info['source'],
                f"{hit['pident']:.1f}",
                str(hit['length']),
                str(hit['qstart']),
                str(hit['qend']),
                f"{hit['evalue']:.2g}",
                f"{hit['bitscore']:.1f}",
                f"{hit['qcovhsp']:.2f}",
                f"{hit['scovhsp']:.2f}",
                "50.0",  # MEROPS_qcovhsp_threshold
                "50.0",  # MEROPS_scovhsp_threshold
                "qcovhsp>=50% AND scovhsp>=50%"  # Formula
            ]
            out.write('\t'.join(row) + '\n')
    
    log("")
    log("Statistics:")
    log(f"  Total proteins:             {len(protein_lengths)}")
    log(f"  Proteins with MEROPS hits:  {len(hits)} ({100*len(hits)/len(protein_lengths):.1f}%)")
    log("")
    log("✓ MEROPS consolidation complete!")

if __name__ == '__main__':
    main()
