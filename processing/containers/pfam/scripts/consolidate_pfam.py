#!/usr/bin/env python3
"""
Consolidate and enrich Pfam protein family annotations
Parses HMMER domain table output and enriches with family descriptions and clan information
"""

import sys
import math
from collections import defaultdict
from datetime import datetime


def log(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}]  {message}", flush=True)

def parse_hmmer_domtblout(hmmer_file):
    """
    Parse HMMER domtblout format (one row per domain).
    
    domtblout format (space-separated):
    0: target name (Pfam family)  1: accession  2: tlen  3: query name (protein ID)  
    4: accession  5: qlen  6: E-value  7: score  8: bias  9: #  10: of  
    11: c-Evalue  12: i-Evalue  13: domain score  14: domain bias  
    15: hmm from  16: hmm to  17: ali from  18: ali to  19: env from  20: env to  
    21: acc  22+: description
    """
    hits = defaultdict(list)
    
    with open(hmmer_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) < 22:
                continue
            
            target_name = parts[0]      # Pfam family name (e.g., ABC_tran)
            accession = parts[1]        # Pfam accession (e.g., PF00001.23)
            tlen = int(parts[2]) if parts[2] != '-' else 0
            query_name = parts[3]       # Protein ID
            qlen = int(parts[5]) if parts[5] != '-' else 0
            full_evalue = float(parts[6]) if parts[6] != '-' else 1.0
            full_score = float(parts[7]) if parts[7] != '-' else 0.0
            full_bias = float(parts[8]) if parts[8] != '-' else 0.0
            domain_cevalue = float(parts[11]) if parts[11] != '-' else 1.0
            domain_ievalue = float(parts[12]) if parts[12] != '-' else 1.0
            domain_score = float(parts[13]) if parts[13] != '-' else 0.0
            domain_bias = float(parts[14]) if parts[14] != '-' else 0.0
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            ali_from = int(parts[17])
            ali_to = int(parts[18])
            env_from = int(parts[19])
            env_to = int(parts[20])
            acc = float(parts[21]) if parts[21] != '-' else 0.0
            
            hits[query_name].append({
                'family_name': target_name,  # Family name from HMMER
                'accession': accession,       # Accession for metadata lookup
                'tlen': tlen,
                'qlen': qlen,
                'full_evalue': full_evalue,
                'full_score': full_score,
                'full_bias': full_bias,
                'domain_cevalue': domain_cevalue,
                'domain_ievalue': domain_ievalue,
                'domain_score': domain_score,
                'domain_bias': domain_bias,
                'hmm_from': hmm_from,
                'hmm_to': hmm_to,
                'ali_from': ali_from,
                'ali_to': ali_to,
                'env_from': env_from,
                'env_to': env_to,
                'acc': acc
            })
    
    return hits

def load_pfam_definitions(pfam_def_file):
    """
    Load Pfam family descriptions and clan information.
    
    Format: accession \t family_name \t description \t clan \t type \t model_length
    Example: PF00001.23 \t ABC_tran \t ABC transporter \t CL0192 \t Domain \t 125
    """
    pfam_desc = {}
    
    try:
        with open(pfam_def_file) as f:
            header = f.readline()  # Skip header
            for line in f:
                if not line.strip():
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 6:
                    accession = parts[0]
                    family_name = parts[1]
                    description = parts[2]
                    clan = parts[3] if parts[3] else ''
                    family_type = parts[4]
                    model_length = parts[5]
                    
                    pfam_desc[accession] = {
                        'family_name': family_name,
                        'description': description,
                        'clan': clan,
                        'type': family_type,
                        'model_length': model_length
                    }
    except FileNotFoundError:
        print(f"Warning: Pfam definitions file not found: {pfam_def_file}", file=sys.stderr)
        print("Continuing without descriptions...", file=sys.stderr)
    
    return pfam_desc

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

def calculate_confidence_scores(hits_data, protein_lengths):
    """
    Calculate confidence scores for Pfam assignments.
    
    Based on:
    - Domain i-Evalue (statistical strength)
    - Coverage/completeness (alignment vs HMM and query)
    - Alignment reliability (HMMER acc score)
    
    Weights: 60% E-value, 30% coverage, 10% alignment accuracy
    """
    for feature_id, hits in hits_data.items():
        for hit in hits:
            # s1: Statistical strength from i-Evalue
            try:
                S = -math.log10(max(hit['domain_ievalue'], 1e-200))
                s1 = min(1.0, max(0.0, S / 50.0))
            except:
                s1 = 0.0
            
            # s2: Coverage/completeness
            try:
                prot_len = protein_lengths.get(feature_id, hit['qlen'])
                Lq = hit['ali_to'] - hit['ali_from'] + 1
                Lh = hit['hmm_to'] - hit['hmm_from'] + 1
                qcov = Lq / prot_len if prot_len > 0 else 0
                hcov = Lh / hit['tlen'] if hit['tlen'] > 0 else 0
                cov_min = min(qcov, hcov)
                s2 = min(1.0, max(0.0, (cov_min - 0.3) / 0.7))
            except:
                s2 = 0.0
            
            # s3: Alignment reliability from HMMER acc
            try:
                s3 = min(1.0, max(0.0, hit['acc']))
            except:
                s3 = 0.0
            
            # Composite confidence: 60% E-value, 30% coverage, 10% accuracy
            confidence = 0.6 * s1 + 0.3 * s2 + 0.1 * s3
            
            hit['evalue_score'] = round(s1, 6)
            hit['coverage_score'] = round(s2, 6)
            hit['alignment_score'] = round(s3, 6)
            hit['confidence'] = round(confidence, 6)

def main():
    if len(sys.argv) != 6:
        log("Usage: consolidate_pfam.py <organism> <hmmer_domains.tsv> <pfam_definitions.txt> <proteins.faa> <output.tsv>")
        sys.exit(1)
    
    organism = sys.argv[1]
    hmmer_output = sys.argv[2]
    pfam_def_file = sys.argv[3]
    protein_file = sys.argv[4]
    output_file = sys.argv[5]
    
    log("="*70)
    log("Pfam Annotation Consolidation")
    log("="*70)
    log("")
    log(f"Organism:              {organism}")
    log(f"HMMER domains file:    {hmmer_output}")
    log(f"Pfam definitions:      {pfam_def_file}")
    log(f"Protein FASTA:         {protein_file}")
    log(f"Output file:           {output_file}")
    log("")
    
    log("Loading reference data...")
    log("")
    
    log("  Loading protein lengths...")
    protein_lengths = get_protein_lengths(protein_file)
    log(f"    Loaded {len(protein_lengths)} proteins")
    
    log("  Loading Pfam definitions...")
    pfam_desc = load_pfam_definitions(pfam_def_file)
    log(f"    Loaded {len(pfam_desc)} family definitions")
    log("")
    
    log("Parsing HMMER domain output...")
    hits_data = parse_hmmer_domtblout(hmmer_output)
    log(f"  Found domain hits for {len(hits_data)} proteins")
    log("")
    
    log("Writing consolidated output...")
    log(f"  Output file: {output_file}")
    log("")    
    # Define output columns (one row per protein, semicolon-separated domains)
    header = [
        'organism', 'feature_id', 'protein_length', 'PFAM_domain_count',
        'PFAM_accession', 'PFAM_family_name', 'PFAM_description', 'PFAM_type', 'PFAM_clan',
        'PFAM_hmm_coords', 'PFAM_ali_coords', 'PFAM_env_coords',
        'PFAM_full_evalue', 'PFAM_full_score',
        'PFAM_domain_ievalue', 'PFAM_domain_score'
    ]
    
    total_domains = 0
    with open(output_file, 'w') as out:
        out.write('\t'.join(header) + '\n')
        
        # Write one row per protein with all domains semicolon-separated
        for feature_id in sorted(hits_data.keys()):
            protein_length = protein_lengths.get(feature_id, 0)
            
            # Sort domains by i-Evalue (best first)
            domains = sorted(hits_data[feature_id], key=lambda x: x['domain_ievalue'])
            
            # Collect all domain values
            accessions = []
            family_names = []
            descriptions = []
            types = []
            clans = []
            hmm_coords = []
            ali_coords = []
            env_coords = []
            full_evalues = []
            full_scores = []
            domain_ievalues = []
            domain_scores = []
            
            for hit in domains:
                # Get metadata from pfam_metadata.tsv using accession
                info = pfam_desc.get(hit['accession'], {})
                family_name = info.get('family_name', hit['family_name'])
                description = info.get('description', '')
                clan = info.get('clan', '')
                family_type = info.get('type', '')
                
                accessions.append(hit['accession'])
                family_names.append(family_name)
                descriptions.append(description if description else 'N/A')
                types.append(family_type if family_type else 'N/A')
                clans.append(clan if clan else 'N/A')
                hmm_coords.append(f"{hit['hmm_from']}-{hit['hmm_to']}")
                ali_coords.append(f"{hit['ali_from']}-{hit['ali_to']}")
                env_coords.append(f"{hit['env_from']}-{hit['env_to']}")
                full_evalues.append(f"{hit['full_evalue']:.2e}")
                full_scores.append(f"{hit['full_score']:.1f}")
                domain_ievalues.append(f"{hit['domain_ievalue']:.2e}")
                domain_scores.append(f"{hit['domain_score']:.1f}")
                total_domains += 1
            
            row = [
                organism,
                feature_id,
                str(protein_length),
                str(len(domains)),
                ';'.join(accessions),
                ';'.join(family_names),
                ';'.join(descriptions),
                ';'.join(types),
                ';'.join(clans),
                ';'.join(hmm_coords),
                ';'.join(ali_coords),
                ';'.join(env_coords),
                ';'.join(full_evalues),
                ';'.join(full_scores),
                ';'.join(domain_ievalues),
                ';'.join(domain_scores)
            ]
            
            out.write('\t'.join(row) + '\n')
    
    # Statistics
    total_proteins = len(protein_lengths)
    proteins_with_hits = len(hits_data)
    unique_families = set(hit['accession'] for hits in hits_data.values() for hit in hits)
    
    log("Output file structure:")
    log("  Columns: 16 (organism, feature_id, domain details, scores)")
    log("  Format: One row per protein, multiple domains semicolon-separated")
    log("")
    log("Results summary:")
    log(f"  Total proteins:             {total_proteins}")
    log(f"  Proteins with Pfam domains: {proteins_with_hits}")
    if total_proteins > 0:
        coverage = 100 * proteins_with_hits / total_proteins
        log(f"  Coverage:                   {coverage:.1f}%")
    log(f"  Total domain hits:          {total_domains}")
    log(f"  Unique Pfam families:       {len(unique_families)}")
    if proteins_with_hits > 0:
        avg_domains = total_domains / proteins_with_hits
        log(f"  Avg domains per protein:    {avg_domains:.2f}")
    log("")
    log("✓ Pfam consolidation complete!")


if __name__ == '__main__':
    main()
