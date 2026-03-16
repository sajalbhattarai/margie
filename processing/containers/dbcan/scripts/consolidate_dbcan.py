#!/usr/bin/env python3
"""
Consolidate and enrich dbCAN CAZyme annotations
Parses HMMER domain table output and enriches with subfamily, EC, and substrate information
"""

import sys
from collections import defaultdict
from pathlib import Path

def _to_roman(num):
    """Convert number to lowercase roman numerals for indexing."""
    val = [1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1]
    syms = ['m', 'cm', 'd', 'cd', 'c', 'xc', 'l', 'xl', 'x', 'ix', 'v', 'iv', 'i']
    roman_num = ''
    i = 0
    while num > 0:
        for _ in range(num // val[i]):
            roman_num += syms[i]
            num -= val[i]
        i += 1
    return roman_num

def parse_hmmer_domtblout(hmmer_file):
    """
    Parse HMMER domtblout format.
    
    domtblout format (space-separated):
    0: target name  1: accession  2: tlen  3: query name  4: accession  5: qlen
    6: E-value  7: score  8: bias  9: #  10: of  11: c-Evalue  12: i-Evalue
    13: domain score  14: domain bias  15: hmm from  16: hmm to  17: ali from  18: ali to
    19: env from  20: env to  21: acc  22+: description
    """
    hits = defaultdict(list)
    
    with open(hmmer_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) < 22:
                continue
            
            target_name = parts[0]   # CAZyme family (e.g., GH1.hmm)
            query_name = parts[3]    # Protein ID
            full_evalue = float(parts[6]) if parts[6] != '-' else 1.0
            full_score = float(parts[7]) if parts[7] != '-' else 0.0
            domain_ievalue = float(parts[12]) if parts[12] != '-' else 1.0
            domain_score = float(parts[13]) if parts[13] != '-' else 0.0
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            ali_from = int(parts[17])
            ali_to = int(parts[18])
            coverage = float(parts[21])
            
            # Extract family name (remove .hmm extension)
            family = target_name.replace('.hmm', '')
            
            hits[query_name].append({
                'family': family,
                'full_evalue': full_evalue,
                'full_score': full_score,
                'domain_evalue': domain_ievalue,
                'domain_score': domain_score,
                'hmm_from': hmm_from,
                'hmm_to': hmm_to,
                'ali_from': ali_from,
                'ali_to': ali_to,
                'coverage': coverage
            })
    
    return hits

def load_fam_subfam_ec(fam_subfam_ec_file):
    """
    Load subfamily and EC number mappings.
    
    Format: subfamily \t protein_id \t ec_number
    Example: AA1_1 \t AAA33103.1 \t 1.10.3.2
    """
    fam_ec_map = defaultdict(set)
    
    with open(fam_subfam_ec_file) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                subfamily = parts[0]
                ec_number = parts[2]
                
                # Extract family from subfamily (e.g., AA1_1 -> AA1)
                family = subfamily.split('_')[0]
                
                if ec_number and ec_number != '-':
                    fam_ec_map[family].add(ec_number)
    
    # Convert sets to sorted lists
    for family in fam_ec_map:
        fam_ec_map[family] = sorted(list(fam_ec_map[family]))
    
    return fam_ec_map

def load_substrate_mapping(substrate_file):
    """
    Load substrate predictions for CAZyme families.
    
    Format (tab-separated):
    Substrate_high_level \t Substrate_Simple \t Family \t Name \t EC_Number
    """
    substrate_map = defaultdict(lambda: {'substrates': set(), 'names': set(), 'ec_numbers': set()})
    
    try:
        with open(substrate_file) as f:
            header = f.readline()  # Skip header
            
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 5:
                    substrate_high = parts[0]
                    substrate_simple = parts[1]
                    family = parts[2]
                    name = parts[3]
                    ec_number = parts[4]
                    
                    if substrate_simple and substrate_simple != '-':
                        substrate_map[family]['substrates'].add(substrate_simple)
                    
                    if name and name != '-':
                        substrate_map[family]['names'].add(name)
                    
                    if ec_number and ec_number != '-':
                        substrate_map[family]['ec_numbers'].add(ec_number)
        
        # Convert sets to sorted lists
        for family in substrate_map:
            substrate_map[family]['substrates'] = sorted(list(substrate_map[family]['substrates']))
            substrate_map[family]['names'] = sorted(list(substrate_map[family]['names']))
            substrate_map[family]['ec_numbers'] = sorted(list(substrate_map[family]['ec_numbers']))
    except FileNotFoundError:
        print(f"Warning: Substrate mapping file not found: {substrate_file}", file=sys.stderr)
        print("Continuing without substrate predictions...", file=sys.stderr)
    
    return substrate_map

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

def calculate_confidence_scores(hits_data):
    """
    Calculate confidence scores for CAZyme assignments.
    
    Based on:
    - Domain E-value (lower is better)
    - Coverage of HMM profile
    - Domain score
    """
    for feature_id, hits in hits_data.items():
        for hit in hits:
            domain_evalue = hit['domain_evalue']
            coverage = hit['coverage']
            score = hit['domain_score']
            
            # E-value confidence (log-scale, normalized)
            if domain_evalue == 0:
                evalue_conf = 1.0
            elif domain_evalue < 1e-10:
                evalue_conf = 1.0
            elif domain_evalue < 1e-5:
                evalue_conf = 0.8
            elif domain_evalue < 1e-3:
                evalue_conf = 0.5
            else:
                evalue_conf = 0.2
            
            # Coverage confidence (0-1)
            coverage_conf = min(coverage, 1.0)
            
            # Score confidence (normalize by typical scores 20-500+)
            score_conf = min(score / 100.0, 1.0)
            
            # Composite confidence: 50% E-value, 30% coverage, 20% score
            confidence = (0.5 * evalue_conf) + (0.3 * coverage_conf) + (0.2 * score_conf)
            
            hit['confidence'] = round(confidence, 6)
            hit['evalue_conf'] = round(evalue_conf, 6)
            hit['coverage_conf'] = round(coverage_conf, 6)
            hit['score_conf'] = round(score_conf, 6)

def main():
    if len(sys.argv) != 7:
        print("Usage: consolidate_dbcan.py <organism> <hmmer_domains.tsv> <fam_subfam_ec.txt> <substrate_mapping.tsv> <proteins.faa> <output.tsv>")
        sys.exit(1)
    
    organism = sys.argv[1]
    hmmer_output = sys.argv[2]
    fam_subfam_ec = sys.argv[3]
    substrate_mapping = sys.argv[4]
    protein_file = sys.argv[5]
    output_file = sys.argv[6]
    
    print(f"Loading protein lengths from {protein_file}...")
    protein_lengths = get_protein_lengths(protein_file)
    print(f"  Loaded {len(protein_lengths)} proteins")
    
    print(f"Loading subfamily/EC mappings from {fam_subfam_ec}...")
    fam_ec_map = load_fam_subfam_ec(fam_subfam_ec)
    print(f"  Loaded EC mappings for {len(fam_ec_map)} families")
    
    print(f"Loading substrate mappings from {substrate_mapping}...")
    substrate_map = load_substrate_mapping(substrate_mapping)
    print(f"  Loaded substrate data for {len(substrate_map)} families")
    
    print(f"Parsing HMMER output from {hmmer_output}...")
    hits_data = parse_hmmer_domtblout(hmmer_output)
    print(f"  Found hits for {len(hits_data)} proteins")
    
    print(f"Writing consolidated output to {output_file}...")
    
    # Define output columns
    header = [
        'organism', 'feature_id', 'protein_length', 'DBCAN_family_count',
        'DBCAN_families', 'DBCAN_family_names', 'DBCAN_domain_positions',
        'DBCAN_substrates', 'DBCAN_EC_numbers', 'DBCAN_domain_evalues',
        'DBCAN_full_evalues', 'DBCAN_domain_scores', 'DBCAN_coverages',
        'DBCAN_best_family', 'DBCAN_best_evalue', 'DBCAN_best_score',
        'DBCAN_threshold_type', 'DBCAN_threshold_value'
    ]
    
    with open(output_file, 'w') as out:
        out.write('\t'.join(header) + '\n')
        
        # Write all proteins (with and without hits)
        for feature_id in sorted(protein_lengths.keys()):
            protein_length = protein_lengths[feature_id]
            
            if feature_id in hits_data:
                hits = hits_data[feature_id]
                family_count = len(hits)
                
                # Sort by best e-value (lowest first)
                hits.sort(key=lambda h: h['domain_evalue'])
                
                # Build numbered lists with roman numerals
                families = '; '.join([f"({_to_roman(i+1)}) {h['family']}" for i, h in enumerate(hits)])
                
                # Get family names from substrate map
                family_names = []
                for i, h in enumerate(hits):
                    names = substrate_map.get(h['family'], {}).get('names', [])
                    if names:
                        family_names.append(f"({_to_roman(i+1)}) {names[0]}")
                    else:
                        family_names.append(f"({_to_roman(i+1)}) ")
                family_names_str = '; '.join(family_names)
                
                # Domain positions (alignment ranges)
                positions = '; '.join([f"({_to_roman(i+1)}) {h['ali_from']}-{h['ali_to']}" for i, h in enumerate(hits)])
                
                # Get substrates
                substrates_list = []
                for i, h in enumerate(hits):
                    subs = substrate_map.get(h['family'], {}).get('substrates', [])
                    if subs:
                        substrates_list.append(f"({_to_roman(i+1)}) {'; '.join(subs)}")
                    else:
                        substrates_list.append(f"({_to_roman(i+1)}) ")
                substrates_str = '; '.join(substrates_list)
                
                # Get EC numbers
                ec_list = []
                for i, h in enumerate(hits):
                    ecs = substrate_map.get(h['family'], {}).get('ec_numbers', [])
                    if not ecs:
                        # Try fam_ec_map as backup
                        ecs = fam_ec_map.get(h['family'], [])
                    if ecs:
                        ec_list.append(f"({_to_roman(i+1)}) {'; '.join(ecs)}")
                    else:
                        ec_list.append(f"({_to_roman(i+1)}) ")
                ec_str = '; '.join(ec_list)
                
                domain_evalues = '; '.join([f"({_to_roman(i+1)}) {h['domain_evalue']}" for i, h in enumerate(hits)])
                full_evalues = '; '.join([f"({_to_roman(i+1)}) {h['full_evalue']}" for i, h in enumerate(hits)])
                domain_scores = '; '.join([f"({_to_roman(i+1)}) {h['domain_score']}" for i, h in enumerate(hits)])
                coverages = '; '.join([f"({_to_roman(i+1)}) {h['coverage']}" for i, h in enumerate(hits)])
                
                # Best hit
                best_hit = hits[0]
                best_family = best_hit['family']
                best_evalue = best_hit['domain_evalue']
                best_score = best_hit['domain_score']
                
                row = [
                    organism,
                    feature_id,
                    str(protein_length),
                    str(family_count),
                    families,
                    family_names_str,
                    positions,
                    substrates_str,
                    ec_str,
                    domain_evalues,
                    full_evalues,
                    domain_scores,
                    coverages,
                    best_family,
                    str(best_evalue),
                    str(best_score),
                    '--cut_ga',  # DBCAN_threshold_type: using gathering threshold
                    'GA'         # DBCAN_threshold_value: GA cutoff from HMM
                ]
            else:
                # No CAZyme hits
                row = [
                    organism,
                    feature_id,
                    str(protein_length),
                    '0',
                    '', '', '', '', '', '', '', '', '',
                    '', '', '', '', ''  # Empty best hit + threshold fields
                ]
            
            out.write('\t'.join(row) + '\n')
    
    # Statistics
    total_proteins = len(protein_lengths)
    proteins_with_hits = len(hits_data)
    total_hits = sum(len(hits) for hits in hits_data.values())
    
    # Count unique families and substrates
    all_families = set()
    all_substrates = set()
    for hits in hits_data.values():
        for h in hits:
            all_families.add(h['family'])
            subs = substrate_map.get(h['family'], {}).get('substrates', [])
            all_substrates.update(subs)
    
    print(f"\nStatistics:")
    print(f"  Total proteins: {total_proteins}")
    print(f"  Proteins with CAZyme domains: {proteins_with_hits} ({100*proteins_with_hits/total_proteins:.1f}%)")
    print(f"  Total CAZyme domain hits: {total_hits}")
    print(f"  Unique CAZyme families found: {len(all_families)}")
    print(f"  Unique substrates predicted: {len(all_substrates)}")
    if proteins_with_hits > 0:
        print(f"  Average domains per CAZyme protein: {total_hits/proteins_with_hits:.2f}")
    print(f"\n✓ Consolidation complete!")

if __name__ == '__main__':
    main()
