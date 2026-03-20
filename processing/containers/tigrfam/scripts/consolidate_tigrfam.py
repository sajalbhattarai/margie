#!/usr/bin/env python3
"""
TIGRfam Domain Consolidation Script

Parses HMMER domtblout format and creates consolidated TSV with:
- One row per domain hit
- TIGRfam family descriptions 
- Confidence scoring based on E-value, coverage, and bit score
- Genome Properties mapping for pathway/system context

Output columns (16):
  1. organism
  2. feature_id
  3. protein_length
  4. TIGRFAM_domain_count
  5. TIGRFAM_accession
  6. TIGRFAM_description
  7. TIGRFAM_isotype
  8. TIGRFAM_tlen
  9. TIGRFAM_hmm_from
  10. TIGRFAM_hmm_to
  11. TIGRFAM_ali_from
  12. TIGRFAM_ali_to
  13. TIGRFAM_env_from
  14. TIGRFAM_env_to
  15. TIGRFAM_full_evalue
  16. TIGRFAM_full_score
  17. TIGRFAM_full_bias
  18. TIGRFAM_domain_cevalue
  19. TIGRFAM_domain_ievalue
  20. TIGRFAM_domain_score
  21. TIGRFAM_domain_bias
  22. TIGRFAM_acc
  23. TIGRFAM_evalue_score
  24. TIGRFAM_coverage_score
  25. TIGRFAM_alignment_score
  26. TIGRFAM_confidence
  27. TIGRFAM_genprop_ids (semicolon-separated)
  28. TIGRFAM_genprop_names (semicolon-separated)
  29. TIGRFAM_genprop_types (SYSTEM/PATHWAY, semicolon-separated)
  30. TIGRFAM_genprop_steps (step names, semicolon-separated)
  31. TIGRFAM_genprop_evidence (necessary/sufficient, semicolon-separated)
"""

import sys
import math
from collections import defaultdict
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
                # Save previous sequence
                if current_id:
                    lengths[current_id] = len(''.join(current_seq))
                # Start new sequence
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_id:
            lengths[current_id] = len(''.join(current_seq))
    
    return lengths

def load_tigrfam_definitions(def_file):
    """Load TIGRfam family definitions from complete_tigrfam_definitions.txt
    
    Format: Accession\tName\tFunction\tIsology Type\tEC Number
    Example: TIGR00001\trpmI_bact\tribosomal protein bL35\tequivalog\t
    """
    definitions = {}
    
    with open(def_file, 'r') as f:
        # Skip header line
        header = f.readline()
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                accession = parts[0]  # TIGR00001
                name = parts[1]       # rpmI_bact
                function = parts[2]   # ribosomal protein bL35
                isotype = parts[3]    # equivalog, subfamily, etc.
                ec = parts[4] if len(parts) > 4 else ""
                
                # Use function as description, include name in parentheses
                description = function
                if name and name != accession and name != "TIGR" + accession.replace("TIGR", ""):
                    description = f"{function} ({name})"
                
                definitions[accession] = {
                    'description': description,
                    'isotype': isotype,
                    'ec': ec
                }
    
    return definitions

def load_genome_properties(mapping_file):
    """Load TIGRfam to Genome Properties mapping
    
    Format: TIGR_accession\tGenProp_ID\tGenProp_name\tGenProp_type\tGenProp_category\tStep_name\tEvidence_type
    Returns: dict mapping TIGR accession to list of genome properties
    """
    tigr_to_genprop = defaultdict(list)
    
    try:
        with open(mapping_file, 'r') as f:
            # Skip header
            header = f.readline()
            
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 7:
                    tigr_acc = parts[0]      # TIGR00006
                    genprop_id = parts[1]    # GenProp0802
                    genprop_name = parts[2]  # Ribosome biogenesis proteins, bacteria
                    genprop_type = parts[3]  # SYSTEM or PATHWAY
                    step_name = parts[5]     # 16S rRNA methyltransferase RsmH
                    evidence = parts[6]      # necessary or sufficient
                    
                    tigr_to_genprop[tigr_acc].append({
                        'id': genprop_id,
                        'name': genprop_name,
                        'type': genprop_type,
                        'step': step_name,
                        'evidence': evidence
                    })
    except FileNotFoundError:
        print("  Warning: Genome properties mapping file not found, skipping...")
        return {}
    
    return tigr_to_genprop

def calculate_confidence(evalue, hmm_coverage, ali_score):
    """
    Calculate confidence score (0-1) based on three components:
    1. E-value score (60% weight): Strong statistical significance
    2. Coverage score (30% weight): How much of HMM is covered
    3. Alignment score (10% weight): Quality of alignment
    """
    # E-value score: -log10(evalue) normalized
    if evalue == 0:
        evalue_score = 1.0
    else:
        log_eval = -math.log10(evalue)
        evalue_score = min(1.0, log_eval / 50.0)  # Cap at 50
    
    # Coverage score: direct percentage
    coverage_score = min(1.0, hmm_coverage)
    
    # Alignment score: normalized bit score
    alignment_score = min(1.0, ali_score / 100.0)
    
    # Weighted combination
    confidence = (0.6 * evalue_score) + (0.3 * coverage_score) + (0.1 * alignment_score)
    
    return round(confidence, 2), round(evalue_score, 6), round(coverage_score, 6), round(alignment_score, 2)

def parse_hmmer_domtblout(domtbl_file, protein_lengths, definitions, genome_properties):
    """Parse HMMER domtblout format"""
    
    # Group domains by protein
    protein_domains = defaultdict(list)
    
    with open(domtbl_file, 'r') as f:
        for line in f:
            # Skip comments
            if line.startswith('#'):
                continue
            
            fields = line.strip().split()
            if len(fields) < 23:
                continue
            
            # Extract fields (0-indexed)
            target_name = fields[0]      # TIGRfam family
            target_acc = fields[1]       # Accession
            query_name = fields[3]       # Protein ID
            full_evalue = float(fields[6])
            full_score = float(fields[7])
            full_bias = float(fields[8])
            dom_cevalue = float(fields[11])
            dom_ievalue = float(fields[12])
            dom_score = float(fields[13])
            dom_bias = float(fields[14])
            hmm_from = int(fields[15])
            hmm_to = int(fields[16])
            ali_from = int(fields[17])
            ali_to = int(fields[18])
            env_from = int(fields[19])
            env_to = int(fields[20])
            tlen = int(fields[2])        # Target length
            
            # Get protein length
            prot_len = protein_lengths.get(query_name, 0)
            
            # Get description and isotype
            fam_info = definitions.get(target_name, {'description': '', 'isotype': '', 'ec': ''})
            description = fam_info['description']
            isotype = fam_info['isotype']
            ec_number = fam_info.get('ec', '')
            
            # Get genome properties for this TIGR family
            genprops = genome_properties.get(target_name, [])
            genprop_ids = ';'.join(gp['id'] for gp in genprops) if genprops else ''
            genprop_names = ';'.join(gp['name'] for gp in genprops) if genprops else ''
            genprop_types = ';'.join(gp['type'] for gp in genprops) if genprops else ''
            genprop_steps = ';'.join(gp['step'] for gp in genprops) if genprops else ''
            genprop_evidence = ';'.join(gp['evidence'] for gp in genprops) if genprops else ''
            
            # Calculate coverage
            hmm_coverage = (hmm_to - hmm_from + 1) / tlen if tlen > 0 else 0
            
            domain = {
                'protein_id': query_name,
                'protein_length': prot_len,
                'family': target_name,
                'description': description,
                'isotype': isotype,
                'tlen': tlen,
                'hmm_from': hmm_from,
                'hmm_to': hmm_to,
                'ali_from': ali_from,
                'ali_to': ali_to,
                'env_from': env_from,
                'env_to': env_to,
                'full_evalue': full_evalue,
                'full_score': full_score,
                'full_bias': full_bias,
                'dom_cevalue': dom_cevalue,
                'dom_ievalue': dom_ievalue,
                'dom_score': dom_score,
                'dom_bias': dom_bias,
                'acc': target_acc,
                'genprop_ids': genprop_ids,
                'genprop_names': genprop_names,
                'genprop_types': genprop_types,
                'genprop_steps': genprop_steps,
                'genprop_evidence': genprop_evidence
            }
            
            protein_domains[query_name].append(domain)
    
    return protein_domains

def main():
    if len(sys.argv) != 6:
        log("Usage: consolidate_tigrfam.py <genome_name> <domtbl_file> <def_file> <faa_file> <output_tsv>")
        sys.exit(1)
    
    genome_name = sys.argv[1]
    domtbl_file = sys.argv[2]
    def_file = sys.argv[3]
    faa_file = sys.argv[4]
    output_file = sys.argv[5]
    
    log("="*70)
    log("TIGRfam Annotation Consolidation")
    log("="*70)
    log("")
    log(f"Organism:              {genome_name}")
    log(f"HMMER domains file:    {domtbl_file}")
    log(f"TIGRfam definitions:   {def_file}")
    log(f"Protein FASTA:         {faa_file}")
    log(f"Output file:           {output_file}")
    log("")
    
    log("Loading reference data...")
    log("")
    
    log("  Loading protein lengths...")
    protein_lengths = load_protein_lengths(faa_file)
    log(f"    Loaded {len(protein_lengths)} proteins")
    
    log("  Loading TIGRfam definitions...")
    definitions = load_tigrfam_definitions(def_file)
    log(f"    Loaded {len(definitions)} family definitions")
    
    # Load genome properties mapping (optional)
    import os
    genprop_file = os.path.join(os.path.dirname(def_file), 'tigr_to_genprop_mapping.tsv')
    log("  Loading genome properties...")
    genome_properties = load_genome_properties(genprop_file)
    if genome_properties:
        log(f"    Loaded {len(genome_properties)} genome property mappings")
    log("")
    
    log("Parsing HMMER domain output...")
    protein_domains = parse_hmmer_domtblout(domtbl_file, protein_lengths, definitions, genome_properties)
    log(f"  Found domain hits for {len(protein_domains)} proteins")
    log("")
    
    log("Writing consolidated output...")
    log(f"  Output file: {output_file}")
    log("")
    
    # Write output TSV (one row per protein, semicolon-separated domains)
    with open(output_file, 'w') as out:
        # Header
        # In TIGRfam HMM files the NAME and ACC fields are both the TIGR ID
        # (e.g. TIGR00006), so storing both produced redundant accession numbers
        # in the old TIGRFAM_id and TIGRFAM_family columns.  The new layout uses
        # a single TIGRFAM_accession column for the TIGR ID and
        # TIGRFAM_description for the human-readable function name.
        header = [
            'organism', 'feature_id', 'protein_length', 'TIGRFAM_domain_count',
            'TIGRFAM_accession', 'TIGRFAM_description', 'TIGRFAM_isotype',
            'TIGRFAM_hmm_coords', 'TIGRFAM_ali_coords', 'TIGRFAM_env_coords',
            'TIGRFAM_full_evalue', 'TIGRFAM_full_score',
            'TIGRFAM_domain_ievalue', 'TIGRFAM_domain_score',
            'TIGRFAM_genprop_ids', 'TIGRFAM_genprop_names'
        ]
        out.write('\t'.join(header) + '\n')

        # Write domains
        total_domains = 0
        unique_families = set()

        for protein_id in sorted(protein_domains.keys()):
            domains = protein_domains[protein_id]

            # Sort domains by position
            domains.sort(key=lambda d: d['env_from'])

            # Collect all domain values
            accessions = []
            descriptions = []
            isotypes = []
            hmm_coords = []
            ali_coords = []
            env_coords = []
            full_evalues = []
            full_scores = []
            domain_ievalues = []
            domain_scores = []
            all_genprop_ids = []
            all_genprop_names = []

            protein_length = domains[0]['protein_length'] if domains else 0

            for domain in domains:
                total_domains += 1
                unique_families.add(domain['family'])

                # Use 'family' (target name) as the canonical TIGR accession;
                # 'acc' (fields[1]) is identical for TIGRfam HMMs.
                accessions.append(domain['family'])
                descriptions.append(domain['description'] if domain['description'] else 'N/A')
                isotypes.append(domain['isotype'] if domain['isotype'] else 'N/A')
                hmm_coords.append(f"{domain['hmm_from']}-{domain['hmm_to']}")
                ali_coords.append(f"{domain['ali_from']}-{domain['ali_to']}")
                env_coords.append(f"{domain['env_from']}-{domain['env_to']}")
                full_evalues.append(f"{domain['full_evalue']:.2e}")
                full_scores.append(f"{domain['full_score']:.1f}")
                domain_ievalues.append(f"{domain['dom_ievalue']:.2e}")
                domain_scores.append(f"{domain['dom_score']:.1f}")
                all_genprop_ids.append(domain['genprop_ids'] if domain['genprop_ids'] else 'N/A')
                all_genprop_names.append(domain['genprop_names'] if domain['genprop_names'] else 'N/A')

            row = [
                genome_name,
                protein_id,
                str(protein_length),
                str(len(domains)),
                ';'.join(accessions),
                ';'.join(descriptions),
                ';'.join(isotypes),
                ';'.join(hmm_coords),
                ';'.join(ali_coords),
                ';'.join(env_coords),
                ';'.join(full_evalues),
                ';'.join(full_scores),
                ';'.join(domain_ievalues),
                ';'.join(domain_scores),
                ';'.join(all_genprop_ids),
                ';'.join(all_genprop_names)
            ]
            out.write('\t'.join(row) + '\n')
    
    # Print statistics
    log("Output file structure:")
    log("  Columns: 17 (organism, feature_id, domain details, genome properties)")
    log("  Format: One row per protein, multiple domains semicolon-separated")
    log("")
    log("Results summary:")
    log(f"  Total proteins:               {len(protein_lengths)}")
    log(f"  Proteins with TIGRfam domains: {len(protein_domains)}")
    if len(protein_lengths) > 0:
        coverage = 100 * len(protein_domains) / len(protein_lengths)
        log(f"  Coverage:                     {coverage:.1f}%")
    log(f"  Total domain hits:            {total_domains}")
    log(f"  Unique TIGRfam families:      {len(unique_families)}")
    if len(protein_domains) > 0:
        avg_domains = total_domains / len(protein_domains)
        log(f"  Avg domains per protein:      {avg_domains:.2f}")
    log("")
    log("✓ TIGRfam consolidation complete!")

if __name__ == '__main__':
    main()
