#!/usr/bin/env python3
"""
Consolidate and enrich KEGG KofamScan annotations
Parses detailed KofamScan output and enriches with KO metadata
Produces standardized TSV format: organism, feature_id, KO assignments, scores, confidence
"""

import argparse
import sys
import re
import math
from collections import defaultdict
from pathlib import Path


def parse_kofamscan_detail(kofamscan_file):
    """
    Parse KofamScan detailed output format.
    
    Format: Each line is a hit
    * gene_name  KO  threshold  score  E-value  KO_definition  (above threshold)
      gene_name  KO  threshold  score  E-value  KO_definition  (below threshold)
    """
    hits = defaultdict(list)
    
    with open(kofamscan_file) as f:
        for line in f:
            line = line.rstrip('\n')
            
            if not line or line.startswith('#'):
                continue
            
            # Check if this is a threshold-passing hit (starts with *)
            above_threshold = line.startswith('*')
            if above_threshold:
                line = line[1:].strip()  # Remove the asterisk
            else:
                line = line.strip()
            
            # Parse the hit line - whitespace separated
            parts = line.split(None, 5)  # Split on whitespace, max 6 parts
            if len(parts) < 5:
                continue
            
            gene_id = parts[0]
            ko_id = parts[1]
            
            # Only include if it looks like a KO ID
            if ko_id.startswith('K'):
                threshold = parts[2].strip() if parts[2].strip() != '-' else '0'
                score = parts[3].strip()
                evalue = parts[4].strip()
                
                hits[gene_id].append({
                    'ko_id': ko_id,
                    'threshold': threshold,
                    'score': score,
                    'evalue': evalue,
                    'above_threshold': above_threshold
                })
    
    return hits


def load_ko_definitions(ko_list_file):
    """
    Load KO definitions from ko_list file.
    
    Format: knum \t threshold \t score_type \t profile_type \t ... \t definition
    """
    ko_info = {}
    
    with open(ko_list_file) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 12:
                ko_id = parts[0]
                definition = parts[11]
                
                # Extract EC numbers from definition
                ec_numbers = []
                ec_matches = re.findall(r'\[EC:([\d\.\-]+)\]', definition)
                if ec_matches:
                    ec_numbers = ec_matches
                
                ko_info[ko_id] = {
                    'definition': definition,
                    'ec_numbers': ';'.join(ec_numbers) if ec_numbers else '',
                    'threshold': parts[1],
                    'score_type': parts[2]
                }
    
    return ko_info


def get_protein_lengths(protein_file):
    """Extract protein lengths from FASTA file."""
    lengths = {}
    current_id = None
    current_seq = []
    
    with open(protein_file) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if current_id:
                    lengths[current_id] = len(''.join(current_seq))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id:
            lengths[current_id] = len(''.join(current_seq))
    
    return lengths


def calculate_confidence_scores(hits_data, ko_info):
    """
    Calculate annotation confidence scores for KofamScan.
    
    KofamScan-specific scoring (uses pre-computed F-measure optimized thresholds):
    - s1 (evalue_score): Statistical strength from E-value
    - s2 (coverage_score): Normalized score/threshold ratio (accounts for model quality)
    - s3 (alignment_score): Not available from KofamScan detail format (set to 0)
    - confidence: 0.6*s1 + 0.3*s2 + 0.1*s3
    
    The score/threshold ratio is meaningful because KEGG thresholds are F-measure optimized
    and account for profile quality, model completeness, and training set characteristics.
    """
    
    for feature_id, hits in hits_data.items():
        for hit in hits:
            evalue = float(hit['evalue']) if hit['evalue'] != '-' else 1.0
            score = float(hit['score'])
            threshold = float(hit['threshold']) if hit['threshold'] not in ['0', '-'] else 100.0
            
            # s1: Statistical Strength from E-value
            try:
                S = -math.log10(max(evalue, 1e-200))
                s1 = min(1.0, max(0.0, S / 50.0))
            except:
                s1 = 0.0
            
            # s2: Coverage/Quality - use normalized score/threshold ratio
            # KofamScan thresholds are F-measure optimized and encode model quality
            try:
                score_ratio = score / threshold if threshold > 0 else 0.0
                s2 = min(1.0, max(0.0, score_ratio))
            except:
                s2 = 0.0
            
            # s3: Alignment reliability - not available from KofamScan detail format  
            s3 = 0.0
            
            # Final confidence (weighted: 0.6, 0.3, 0.1)
            confidence = 0.6*s1 + 0.3*s2 + 0.1*s3
            
            hit['evalue_score'] = round(s1, 6)
            hit['coverage_score'] = round(s2, 6)
            hit['alignment_score'] = round(s3, 6)
            hit['confidence'] = round(confidence, 6)


def consolidate_kegg_annotations(organism, kofamscan_output, ko_list, protein_file, output_file):
    """Main consolidation function."""
    
    print(f"Loading protein lengths from {protein_file}...")
    protein_lengths = get_protein_lengths(protein_file)
    print(f"  Loaded {len(protein_lengths)} proteins")
    
    print(f"Loading KO definitions from {ko_list}...")
    ko_info = load_ko_definitions(ko_list)
    print(f"  Loaded {len(ko_info)} KO definitions")
    
    print(f"Parsing KofamScan output from {kofamscan_output}...")
    hits_data = parse_kofamscan_detail(kofamscan_output)
    print(f"  Found hits for {len(hits_data)} proteins")
    
    print(f"Writing consolidated output to {output_file}...")
    
    # Define output columns
    header = [
        'organism',
        'feature_id',
        'protein_length',
        'KEGG_KO_count',
        'KEGG_KO_IDs',
        'KEGG_KO_definitions',
        'KEGG_EC_numbers',
        'KEGG_scores',
        'KEGG_thresholds',
        'KEGG_threshold_type',
        'KEGG_evalues',
        'KEGG_score_types',
        'KEGG_best_KO',
        'KEGG_best_score',
        'KEGG_best_evalue'
    ]
    
    with open(output_file, 'w') as out:
        out.write('\t'.join(header) + '\n')
        
        # Write all proteins (hits and non-hits)
        for feature_id in sorted(protein_lengths.keys()):
            protein_length = protein_lengths[feature_id]
            all_hits = hits_data.get(feature_id, [])
            
            # Filter to only above-threshold hits (marked with * in KofamScan output)
            hits = [h for h in all_hits if h.get('above_threshold', False)]
            
            if hits:
                # Sort hits by score (descending)
                hits = sorted(hits, key=lambda x: float(x['score']), reverse=True)
                
                ko_count = len(hits)
                ko_ids = ';'.join([h['ko_id'] for h in hits])
                ko_defs = ';'.join([ko_info.get(h['ko_id'], {}).get('definition', '') for h in hits])
                ec_numbers = ';'.join([ko_info.get(h['ko_id'], {}).get('ec_numbers', '') for h in hits])
                scores = ';'.join([h['score'] for h in hits])
                thresholds = ';'.join([h['threshold'] for h in hits])
                evalues = ';'.join([h['evalue'] for h in hits])
                score_types = ';'.join([ko_info.get(h['ko_id'], {}).get('score_type', '') for h in hits])
                
                # Best hit (highest score)
                best_hit = hits[0]
                best_ko = best_hit['ko_id']
                best_score = best_hit['score']
                best_evalue = best_hit['evalue']
                
                row = [
                    organism,
                    feature_id,
                    str(protein_length),
                    str(ko_count),
                    ko_ids,
                    ko_defs,
                    ec_numbers,
                    scores,
                    thresholds,  # KO-specific curated thresholds from ko_list
                    'KO-specific curated',  # KEGG_threshold_type: each KO has curated HMM threshold
                    evalues,
                    score_types,
                    best_ko,
                    best_score,
                    best_evalue
                ]
            else:
                # No hits
                row = [
                    organism,
                    feature_id,
                    str(protein_length),
                    '0',
                    '', '', '', '', '', '', '', '', '', '', ''
                ]
            
            out.write('\t'.join(row) + '\n')
    
    # Statistics
    total_proteins = len(protein_lengths)
    proteins_with_hits = sum(1 for fid, hits in hits_data.items() 
                            if any(h.get('above_threshold', False) for h in hits))
    total_hits = sum(sum(1 for h in hits if h.get('above_threshold', False)) 
                    for hits in hits_data.values())
    
    print(f"\nStatistics:")
    print(f"  Total proteins: {total_proteins}")
    print(f"  Proteins with KO assignments: {proteins_with_hits} ({100*proteins_with_hits/total_proteins:.1f}%)")
    print(f"  Total KO assignments: {total_hits}")
    if proteins_with_hits > 0:
        print(f"  Average KOs per annotated protein: {total_hits/proteins_with_hits:.2f}")
    print(f"\n✓ Consolidation complete!")


def main():
    parser = argparse.ArgumentParser(
        description='Consolidate and enrich KEGG KofamScan annotations'
    )
    parser.add_argument('output_dir',
                       help='Output directory containing kofamscan_output.txt')
    parser.add_argument('genome_name',
                       help='Genome name for output file naming')
    
    args = parser.parse_args()
    
    # Construct paths (container paths)
    output_dir = args.output_dir
    genome_name = args.genome_name
    
    kofamscan_output = f"{output_dir}/native/kofamscan_output.txt"
    ko_list = "/container/db/ko_list"
    protein_file = f"{output_dir}/gene_calls/{genome_name}.faa"
    output_file = f"{output_dir}/kegg.tsv"
    
    # Validate inputs
    if not Path(kofamscan_output).exists():
        print(f"ERROR: KofamScan output not found: {kofamscan_output}")
        sys.exit(1)
    
    if not Path(ko_list).exists():
        print(f"ERROR: KO list not found: {ko_list}")
        sys.exit(1)
    
    if not Path(protein_file).exists():
        print(f"ERROR: Protein file not found: {protein_file}")
        sys.exit(1)
    
    consolidate_kegg_annotations(
        genome_name,
        kofamscan_output,
        ko_list,
        protein_file,
        output_file
    )


if __name__ == '__main__':
    main()
