#!/usr/bin/env python3
"""
Consolidate UniProt DIAMOND output with GO terms and metadata.
Maps proteins to Swiss-Prot curated annotations with Gene Ontology terms.
"""

import sys
import csv
import math
import re
from datetime import datetime

def log(message):
    """Log with timestamp."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}]  {message}", flush=True)

def load_go_mapping(go_mapping_file, needed_accessions):
    """
    Load UniProt to GO term mapping ONLY for accessions we have hits for.
    This avoids loading the entire 6GB file into memory.
    Format: uniprot_id \t GO:term1; GO:term2; ...
    Returns dict: uniprot_id -> [list of GO terms]
    """
    go_map = {}
    log(f"Loading UniProt to GO mapping for {len(needed_accessions)} accessions...")
    log(f"  (Scanning large file, this may take a few minutes)")
    
    accession_set = set(needed_accessions)
    loaded_count = 0
    
    with open(go_mapping_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                uniprot_id = parts[0]
                if uniprot_id in accession_set:
                    go_terms = parts[1].split('; ') if parts[1] else []
                    go_map[uniprot_id] = go_terms
                    loaded_count += 1
                    # Early exit if we found all needed accessions
                    if loaded_count == len(needed_accessions):
                        break
    
    log(f"  Loaded GO mappings for {loaded_count} proteins")
    return go_map

def load_protein_lengths(fasta_file):
    """
    Load protein IDs and sequence lengths from input FASTA.
    Returns dict: protein_id -> length
    """
    proteins = {}
    log(f"Loading protein lengths from input FASTA...")
    
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
    
    log(f"  Loaded {len(proteins)} proteins")
    return proteins

def parse_uniprot_stitle(stitle):
    """
    Parse UniProt stitle to extract accession, name, organism, and description.
    Format: sp|P12345|NAME_ORGANISM Description OS=Organism name OX=taxid GN=gene PE=level SV=version
           tr|P12345|NAME_ORGANISM Description OS=Organism name OX=taxid GN=gene PE=level SV=version
    """
    result = {
        'database': '',
        'accession': '',
        'entry_name': '',
        'protein_name': '',
        'organism': '',
        'gene': ''
    }
    
    if not stitle:
        return result
    
    # Extract database type (sp or tr) and identifiers
    match = re.match(r'(sp|tr)\|([A-Z0-9]+)\|([A-Z0-9_]+)', stitle)
    if match:
        result['database'] = 'SwissProt' if match.group(1) == 'sp' else 'TrEMBL'
        result['accession'] = match.group(2)
        result['entry_name'] = match.group(3)
    
    # Extract protein name (before OS=)
    name_match = re.search(r'\|([^\|]+?)\s+OS=', stitle)
    if name_match:
        result['protein_name'] = name_match.group(1).strip()
    
    # Extract organism (OS=...)
    org_match = re.search(r'OS=([^=]+?)(?:\s+OX=|\s+GN=|\s+PE=|\s+SV=|$)', stitle)
    if org_match:
        result['organism'] = org_match.group(1).strip()
    
    # Extract gene name (GN=...)
    gene_match = re.search(r'GN=([^=]+?)(?:\s+PE=|\s+SV=|$)', stitle)
    if gene_match:
        result['gene'] = gene_match.group(1).strip()
    
    return result

def parse_diamond_output(diamond_file):
    """
    Parse DIAMOND output file.
    Returns dict: query_id -> {subject_id, pident, length, evalue, bitscore, qlen, slen, qcovhsp, scovhsp, stitle}
    Only keeps best hit per query.
    Format: qseqid sseqid pident length qlen slen evalue bitscore qcovhsp scovhsp stitle
    """
    hits = {}
    log(f"Parsing DIAMOND output...")
    
    with open(diamond_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 11:
                continue
            
            query_id = parts[0]
            
            # Only keep first (best) hit for each query
            if query_id in hits:
                continue
            
            hits[query_id] = {
                'subject_id': parts[1],
                'pident': float(parts[2]),
                'length': int(parts[3]),
                'qlen': int(parts[4]),
                'slen': int(parts[5]),
                'evalue': float(parts[6]) if parts[6] != '0' else 0.0,
                'bitscore': float(parts[7]),
                'qcovhsp': float(parts[8]),
                'scovhsp': float(parts[9]),
                'stitle': parts[10] if len(parts) > 10 else parts[1]  # Fallback to subject_id if no stitle
            }
    
    log(f"  Parsed {len(hits)} UniProt hits")
    return hits

def calculate_scoring_components(evalue, qstart, qend, protein_length, pident):
    """
    Calculate MARGIE scoring components.
    
    evalue_score: Statistical strength from E-value [0-1]
    coverage_score: Query coverage [0-1]
    identity_score: Sequence similarity [0-1]
    """
    # E-value score (statistical strength)
    if evalue == 0.0:
        evalue_score = 1.0
    else:
        evalue_score = min(1.0, max(0.0, -math.log10(evalue) / 50.0))
    
    # Coverage score
    query_coverage = (qend - qstart + 1) / protein_length if protein_length > 0 else 0
    coverage_score = min(1.0, max(0.0, (query_coverage - 0.3) / 0.7))
    
    # Identity score
    identity_score = min(1.0, max(0.0, pident / 100.0))
    
    return evalue_score, coverage_score, identity_score

def calculate_confidence_score(evalue_score, coverage_score, identity_score):
    """
    Calculate final confidence score as weighted sum.
    Weights: E-value (0.6), Coverage (0.3), Identity (0.1)
    """
    return 0.6 * evalue_score + 0.3 * coverage_score + 0.1 * identity_score

def consolidate_uniprot_annotations(diamond_file, fasta_file, go_mapping_file, 
                                    organism_name, output_file):
    """
    Main consolidation function.
    Creates comprehensive annotation file with GO terms.
    """
    log(f"UniProt Annotation Consolidation")
    log(f"="*60)
    log(f"Organism:        {organism_name}")
    log(f"DIAMOND output:  {diamond_file}")
    log(f"Input FASTA:     {fasta_file}")
    log(f"GO mapping:      {go_mapping_file}")
    log(f"Output file:     {output_file}")
    log("")
    
    # Load protein lengths and hits first
    proteins = load_protein_lengths(fasta_file)
    log("")
    hits = parse_diamond_output(diamond_file)
    log("")
    
    # Extract accessions we need GO terms for
    needed_accessions = set()
    for hit in hits.values():
        stitle_info = parse_uniprot_stitle(hit['stitle'])
        if stitle_info['accession']:
            needed_accessions.add(stitle_info['accession'])
    
    log(f"Need GO terms for {len(needed_accessions)} unique accessions")
    
    log("")
    # Only load GO terms for proteins we have hits for
    go_map = load_go_mapping(go_mapping_file, needed_accessions)
    log("")
    
    # Write consolidated output
    log(f"Writing consolidated results...")
    
    with open(output_file, 'w') as out:
        # Header
        header = [
            'organism', 'feature_id', 'protein_length',
            'UNIPROT_database', 'UNIPROT_accession', 'UNIPROT_entry_name', 'UNIPROT_protein_name',
            'UNIPROT_organism', 'UNIPROT_gene', 'UNIPROT_go_terms',
            'UNIPROT_pident', 'UNIPROT_alignment_length',
            'UNIPROT_evalue', 'UNIPROT_bitscore', 'UNIPROT_qcovhsp', 'UNIPROT_scovhsp',
            'UNIPROT_qcovhsp_threshold', 'UNIPROT_scovhsp_threshold', 'UNIPROT_threshold_formula'
        ]
        out.write('\t'.join(header) + '\n')
        
        # Process all proteins
        hit_count = 0
        for protein_id in sorted(proteins.keys()):
            protein_length = proteins[protein_id]
            
            if protein_id in hits:
                hit = hits[protein_id]
                hit_count += 1
                
                # Parse stitle to extract metadata
                stitle_info = parse_uniprot_stitle(hit['stitle'])
                
                # Get GO terms for this UniProt entry
                accession = stitle_info['accession']
                go_terms = go_map.get(accession, [])
                go_terms_str = '; '.join(go_terms[:20]) if go_terms else ''  # Limit to 20 GO terms
                
                row = [
                    organism_name,
                    protein_id,
                    protein_length,
                    stitle_info['database'],
                    accession,
                    stitle_info['entry_name'],
                    stitle_info['protein_name'],
                    stitle_info['organism'],
                    stitle_info['gene'],
                    go_terms_str,
                    f"{hit['pident']:.1f}",
                    hit['length'],
                    f"{hit['evalue']:.2e}" if hit['evalue'] > 0 else "0.0",
                    f"{hit['bitscore']:.1f}",
                    f"{hit['qcovhsp']:.2f}",
                    f"{hit['scovhsp']:.2f}",
                    "50.0",  # UNIPROT_qcovhsp_threshold
                    "50.0",  # UNIPROT_scovhsp_threshold
                    "qcovhsp>=50% AND scovhsp>=50%"  # Formula
                ]
            else:
                # No hit for this protein
                row = [
                    organism_name,
                    protein_id,
                    protein_length,
                    '', '', '', '', '', '', '',  # Empty UniProt fields
                    '', '', '', '', '', '',  # Empty alignment fields (no qstart/qend)
                    '', '', ''  # Empty threshold fields
                ]
            
            out.write('\t'.join(str(x) for x in row) + '\n')
    
    # Print summary
    log("")
    log(f"Statistics:")
    log(f"  Total proteins:             {len(proteins)}")
    log(f"  Proteins with UniProt hits: {hit_count} ({100*hit_count/len(proteins):.1f}%)")
    log(f"  Unique UniProt entries:     {len(needed_accessions)}")
    log(f"  Proteins with GO terms:     {sum(1 for h in hits.values() if go_map.get(parse_uniprot_stitle(h['stitle'])['accession']))}")
    log("")
    log(f"✓ UniProt consolidation complete!")

if __name__ == '__main__':
    if len(sys.argv) != 6:
        print("Usage: consolidate_uniprot.py <diamond_tsv> <input_faa> <go_mapping> <organism> <output_tsv>")
        sys.exit(1)
    
    consolidate_uniprot_annotations(
        diamond_file=sys.argv[1],
        fasta_file=sys.argv[2],
        go_mapping_file=sys.argv[3],
        organism_name=sys.argv[4],
        output_file=sys.argv[5]
    )
