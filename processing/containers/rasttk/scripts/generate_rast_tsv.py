#!/usr/bin/env python3
"""
Consolidate RASTtk annotation outputs into a single rast.tsv file.

This script merges CDS, RNA, and prophage annotations with sequences,
extracting EC numbers and parsing gene locations.
It also enriches the output with complete SEED Subsystem variant definitions,
including curator metadata, variant descriptions, notes, and versions.

Usage:
    python consolidate_rast_output.py <rast_output_dir> <organism_name>

Example:
    python consolidate_rast_output.py output/rasttk/theta "Bacteroides thetaiotaomicron"
"""

import sys
import os
import re
import csv
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import defaultdict

# -----------------------------------------------------------------------------
# Subsystem Enrichment Logic with SEED Variant Definitions
# -----------------------------------------------------------------------------

def load_subsystem_mapping(mapping_file: str) -> Dict[str, List[Dict[str, str]]]:
    """
    Loads subsystem mapping from a TSV file.
    Format: subsystem_id, roles (:: sep), subsystem_name, superclass, class, subclass
    Returns: Dict { role_name: [ {subsystem_name, superclass, class, subclass} ] }
    """
    role_to_subsystems = {}
    
    if not os.path.exists(mapping_file):
        print(f"Warning: Subsystem mapping file not found at {mapping_file}. Skipping enrichment.")
        return role_to_subsystems

    print(f"Loading subsystem mapping from {mapping_file}...")
    try:
        with open(mapping_file, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            # Check header
            headers = next(reader, None)
            
            for row in reader:
                if len(row) < 6:
                    continue
                
                # Parse row based on p3-all-subsystems output
                # 0: subsystem_id
                # 1: roles (:: separated)
                # 2: subsystem_name
                # 3: superclass
                # 4: class
                # 5: subclass
                
                roles_str = row[1]
                subsystem_info = {
                    'name': row[2],
                    'superclass': row[3],
                    'class': row[4],
                    'subclass': row[5]
                }
                
                # Split roles by '::'
                roles = roles_str.split('::')
                for role in roles:
                    role_clean = role.strip()
                    if not role_clean:
                        continue
                    
                    if role_clean not in role_to_subsystems:
                        role_to_subsystems[role_clean] = []
                    role_to_subsystems[role_clean].append(subsystem_info)
                    
    except Exception as e:
        print(f"Error loading subsystem mapping: {e}")

    print(f"Loaded subsystem mapping for {len(role_to_subsystems)} roles.")
    return role_to_subsystems


def load_subsystem_metadata(metadata_file: str) -> Dict[str, Dict[str, str]]:
    """
    Load subsystem metadata (curator, description, notes, version)
    Returns: Dict[subsystem_name] -> {curator, description, notes, version}
    """
    metadata = {}
    
    if not os.path.exists(metadata_file):
        return metadata
    
    try:
        with open(metadata_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                subsystem = row.get('subsystem', '')
                if subsystem:
                    metadata[subsystem] = {
                        'curator': row.get('curator', ''),
                        'description': row.get('description', '').replace('\n', ' ').strip(),
                        'notes': row.get('notes', '').replace('\n', ' ').strip(),
                        'version': row.get('version', '')
                    }
    except Exception as e:
        print(f"Error loading subsystem metadata: {e}")
    
    return metadata


def load_variant_roles_mapping(roles_file: str) -> Dict[Tuple[str, str], List[str]]:
    """
    Load role assignments for each variant
    Returns: Dict[(subsystem, variant)] -> [list of role names]
    """
    variant_roles = defaultdict(list)
    
    if not os.path.exists(roles_file):
        return variant_roles
    
    try:
        with open(roles_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                subsystem = row.get('subsystem', '')
                variant = row.get('variant', '')
                present = row.get('present', 'False')
                role_name = row.get('role_name', '')
                
                if subsystem and variant and present == 'True' and role_name:
                    key = (subsystem, variant)
                    variant_roles[key].append(role_name)
    except Exception as e:
        print(f"Error loading variant roles: {e}")
    
    return dict(variant_roles)


def load_subsystem_variants(variants_file: str) -> Dict[Tuple[str, str], Dict]:
    """
    Load variant definitions
    Returns: Dict[(subsystem, variant)] -> {description, is_complete, ...}
    """
    variants = {}
    
    if not os.path.exists(variants_file):
        return variants
    
    try:
        with open(variants_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                subsystem = row.get('subsystem', '')
                variant = row.get('variant', '')
                if subsystem and variant:
                    key = (subsystem, variant)
                    variants[key] = {
                        'description': row.get('description', ''),
                        'is_complete': row.get('is_complete', ''),
                        'total_roles': row.get('total_roles', ''),
                        'present_roles': row.get('present_roles', '')
                    }
    except Exception as e:
        print(f"Error loading subsystem variants: {e}")
    
    return variants


def match_variant_for_role(
    role: str,
    subsystem_name: str,
    variant_roles: Dict[Tuple[str, str], List[str]]
) -> Tuple[str, str]:
    """
    Match a role to a specific variant
    First tries the specified subsystem, then searches all subsystems
    Returns: (variant_number, actual_subsystem_used)
    """
    # First try: exact subsystem match (preferred)
    for (sub, variant), roles in variant_roles.items():
        if sub == subsystem_name and role in roles:
            return variant, sub
    
    # Second try: search all subsystems for this role
    # (handles cases where subsystem name differs between mapping sources)
    for (sub, variant), roles in variant_roles.items():
        if role in roles:
            return variant, sub
    
    return '', ''

# -----------------------------------------------------------------------------
# Existing Logic Helper Functions
# -----------------------------------------------------------------------------

def parse_gene_location(gene_id: str) -> Tuple[int, int, str]:
    """
    Parse gene location from gene_id string.
    
    Examples:
        NODE_1_length_527240_cov_20.279845_15807+1122 -> (15807, 16929, '+')
        NODE_1_length_527240_cov_20.279845_15807-1122 -> (14685, 15807, '-')
    
    Returns:
        (start, end, strand)
    """
    # Extract the location part (last component after final underscore)
    match = re.search(r'_(\d+)([+-])(\d+)$', gene_id)
    if match:
        start = int(match.group(1))
        strand = match.group(2)
        length = int(match.group(3))
        
        if strand == '+':
            end = start + length
        else:  # strand == '-'
            end = start
            start = start - length
        
        return start, end, strand
    
    return 0, 0, '?'


def extract_ec_numbers(description: str) -> str:
    """
    Extract EC numbers from description string.
    
    Examples:
        "Aerotolerance protein BatE (EC 4.3.3.7)" -> "EC 4.3.3.7"
        "Some protein (EC 1.1.1.1) and (EC 2.2.2.2)" -> "EC 1.1.1.1; EC 2.2.2.2"
    """
    ec_pattern = r'\(EC[\s]+[\d\.\-]+\)'
    ec_matches = re.findall(ec_pattern, description)
    
    if ec_matches:
        # Clean up and join multiple EC numbers
        ec_numbers = [ec.strip('()').strip() for ec in ec_matches]
        return '; '.join(ec_numbers)
    
    return ''


def load_sequences(fasta_file: str) -> Dict[str, str]:
    """
    Load sequences from FASTA file.
    
    Returns:
        Dictionary mapping feature_id to sequence
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    if not os.path.exists(fasta_file):
        return sequences
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                
                # Start new sequence
                # Header format: >fig|6666666.1476522.peg.1 function
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences


def load_tsv_data(tsv_file: str) -> List[List[str]]:
    """Load TSV file data."""
    data = []
    if not os.path.exists(tsv_file):
        return data
    
    with open(tsv_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                data.append(line.split('\t'))
    
    return data


def consolidate_rast_annotations(rast_dir: str, organism_name: str, output_file: str):
    """
    Consolidate all RASTtk outputs into a single rast.tsv file with complete SEED enrichment.
    """
    rast_path = Path(rast_dir)
    
    # -------------------------------------------------------------------------
    # 1. Setup Enrichment - Load all SEED data
    # -------------------------------------------------------------------------
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.abspath(os.path.join(script_dir, '../'))
    
    # Paths for subsystem mapping and variant definitions
    # Priority: 1) Container path, 2) Local bvbrc_container/database, 3) Legacy path
    mapping_path_container = "/container/db/subsystem_mapping.tsv"
    mapping_path_bvbrc = os.path.join(base_dir, "database/subsystem_mapping.tsv")
    mapping_path_legacy = os.path.join(base_dir, "../processing/database/subsystem_mapping.tsv")
    
    variant_defs_container = "/container/db/variant_definitions"
    variant_defs_bvbrc = os.path.join(base_dir, "database/variant_definitions")
    variant_defs_legacy = os.path.join(base_dir, "../processing/database/variant_definitions")
    
    # Check which paths exist
    if os.path.exists(mapping_path_container):
        mapping_file = mapping_path_container
        variant_defs_dir = variant_defs_container
        print(f"Using container database paths")
    elif os.path.exists(mapping_path_bvbrc):
        mapping_file = mapping_path_bvbrc
        variant_defs_dir = variant_defs_bvbrc
        print(f"Using bvbrc_container database paths: {base_dir}/database/")
    elif os.path.exists(mapping_path_legacy):
        mapping_file = mapping_path_legacy
        variant_defs_dir = variant_defs_legacy
        print(f"Using legacy database paths")
    else:
        mapping_file = None
        variant_defs_dir = None
        print("Warning: No database paths found!")
    
    # Load subsystem mapping (role -> subsystem hierarchy)
    role_to_subsystems = {}
    if mapping_file:
        role_to_subsystems = load_subsystem_mapping(mapping_file)
    else:
        print("No subsystem mapping file found. Proceeding without enrichment.")
    
    # Load SEED variant data
    subsystem_metadata = {}
    variant_roles = {}
    subsystem_variants = {}
    
    if variant_defs_dir and os.path.exists(variant_defs_dir):
        print("Loading SEED variant definitions...")
        
        metadata_file = os.path.join(variant_defs_dir, 'local_references', 'subsystem_metadata.tsv')
        if os.path.exists(metadata_file):
            subsystem_metadata = load_subsystem_metadata(metadata_file)
            print(f"  Loaded metadata for {len(subsystem_metadata)} subsystems")
        
        roles_file = os.path.join(variant_defs_dir, 'structured_database', 'tsv_tables', 'variant_roles.tsv')
        if os.path.exists(roles_file):
            variant_roles = load_variant_roles_mapping(roles_file)
            print(f"  Loaded role assignments for {len(variant_roles)} variants")
        
        variants_file = os.path.join(variant_defs_dir, 'structured_database', 'tsv_tables', 'subsystem_variants.tsv')
        if os.path.exists(variants_file):
            subsystem_variants = load_subsystem_variants(variants_file)
            print(f"  Loaded {len(subsystem_variants)} variant definitions")
    else:
        print("SEED variant definitions not found. Will only include basic subsystem info.")

    # -------------------------------------------------------------------------
    # 2. Identify Input Files
    # -------------------------------------------------------------------------
    # Check if files are in 'native' subdirectory (standard pipeline structure)
    search_path = rast_path / "native"
    if not search_path.exists():
        search_path = rast_path

    # Find the genome prefix (e.g., "theta" from "theta_CDS.tsv")
    cds_files = list(search_path.glob('*_CDS.tsv'))
    if not cds_files:
        print(f"Error: No CDS.tsv file found in {search_path}", file=sys.stderr)
        return
    
    genome_prefix = cds_files[0].stem.replace('_CDS', '')
    
    # Define file paths in the correct directory
    cds_file = search_path / f"{genome_prefix}_CDS.tsv"
    rna_file = search_path / f"{genome_prefix}_RNA.tsv"
    prophage_file = search_path / f"{genome_prefix}_prophage.tsv"
    
    # Sequences might be in 'gene_calls' or 'native' depending on when this is run
    # Pipeline moves them to gene_calls at the end.
    # Check gene_calls first for sequences
    gene_calls_path = rast_path / "gene_calls"
    if (gene_calls_path / f"{genome_prefix}.faa").exists():
        faa_file = gene_calls_path / f"{genome_prefix}.faa"
        ffn_file = gene_calls_path / f"{genome_prefix}.ffn"
    else:
        # Fallback to search path (native or root)
        faa_file = search_path / f"{genome_prefix}.faa"
        ffn_file = search_path / f"{genome_prefix}.ffn"
    
    print(f"Loading sequences from {faa_file}...")
    aa_sequences = load_sequences(str(faa_file))
    
    print(f"Loading sequences from {ffn_file}...")
    na_sequences = load_sequences(str(ffn_file))
    
    # Prepare output data
    all_rows = []
    
    # Helper to process features
    def process_features(file_path, feature_type_override=None):
        if not os.path.exists(file_path):
            return 0
            
        print(f"Processing features from {file_path}...")
        data = load_tsv_data(str(file_path))
        count = 0
        
        for row in data:
            if len(row) < 4:
                continue
            
            feature_id = row[0]
            gene_id = row[1] # Contains location
            gene_type = feature_type_override if feature_type_override else (row[2] if len(row) > 2 else '')
            description = row[3] if len(row) > 3 else ''
            feature_hash = row[5] if len(row) > 5 else ''
            
            # Parse location
            gene_start, gene_end, strand = parse_gene_location(gene_id)
            
            # Extract EC numbers
            ec_numbers = extract_ec_numbers(description)
            
            # Get sequences
            aa_seq = aa_sequences.get(feature_id, '')
            na_seq = na_sequences.get(feature_id, '')
            
            # Calculate lengths
            aa_length = len(aa_seq)
            na_length = len(na_seq)
            
            # --------------------------------
            # Subsystem Enrichment with SEED Variant Data
            # --------------------------------
            sub_superclass = []
            sub_class = []
            sub_subclass = []
            sub_name = []
            
            # SEED variant columns
            seed_variant_num = ''
            seed_variant_def = ''
            seed_curator = ''
            seed_description = ''
            seed_notes = ''
            seed_version = ''
            seed_subsystem_notes = ''
            
            # Lookup role (description) in dictionary
            if description in role_to_subsystems:
                for sub in role_to_subsystems[description]:
                    if sub['superclass']: sub_superclass.append(sub['superclass'])
                    if sub['class']: sub_class.append(sub['class'])
                    if sub['subclass']: sub_subclass.append(sub['subclass'])
                    if sub['name']: sub_name.append(sub['name'])
                
                # Get SEED variant data for the first (primary) subsystem
                if sub_name:
                    primary_subsystem = sub_name[0]
                    
                    # Match to variant (may find in different subsystem than mapping)
                    variant_num, matched_subsystem = match_variant_for_role(description, primary_subsystem, variant_roles)
                    
                    # Use the matched subsystem for metadata lookup
                    metadata_subsystem = matched_subsystem if matched_subsystem else primary_subsystem
                    
                    # Get metadata
                    if metadata_subsystem in subsystem_metadata:
                        meta = subsystem_metadata[metadata_subsystem]
                        seed_curator = meta['curator']
                        seed_description = meta['description']
                        seed_notes = meta['notes']
                        seed_version = meta['version']
                    
                    if variant_num:
                        seed_variant_num = variant_num
                        
                        # Get variant definition using the matched subsystem
                        var_key = (metadata_subsystem, variant_num)
                        if var_key in subsystem_variants:
                            var_info = subsystem_variants[var_key]
                            seed_variant_def = var_info['description']
                            seed_subsystem_notes = f"Roles: {var_info['present_roles']}/{var_info['total_roles']}"
            
            # Join unique values
            def join_unique(items):
                return "; ".join(sorted(list(set(items))))
            
            row_data = [
                organism_name,
                feature_id,
                gene_id,
                str(gene_start),
                str(gene_end),
                str(na_length),
                str(aa_length),
                na_seq,
                aa_seq,
                gene_type,
                description,
                ec_numbers,
                strand,
                feature_hash,
                os.path.basename(str(file_path)),
                json.dumps(row, ensure_ascii=False),
                join_unique(sub_superclass),
                join_unique(sub_class),
                join_unique(sub_subclass),
                join_unique(sub_name),
                seed_variant_num,
                seed_variant_def,
                seed_curator,
                seed_description,
                seed_notes,
                seed_version,
                seed_subsystem_notes
            ]
            all_rows.append(row_data)
            count += 1
            
        return count

    # Process all types
    n_cds = process_features(cds_file)
    n_rna = process_features(rna_file)
    n_prophage = process_features(prophage_file, feature_type_override='prophage')
    
    # Write output file
    print(f"\nWriting consolidated output to {output_file}...")
    with open(output_file, 'w', newline='') as f:
        # Write header with all columns
        header = [
            'organism',
            'feature_id',
            'gene_id',
            'gene_start',
            'gene_end',
            'na_length',
            'aa_length',
            'na_seq',
            'aa_seq',
            'RAST_feature_type',
            'RAST_description',
            'RAST_EC_numbers',
            'RAST_strand',
            'RAST_feature_hash',
            'RAST_source_file',
            'RAST_genecaller_raw_row',
            'RAST_Subsystem_Superclass',
            'RAST_Subsystem_Class',
            'RAST_Subsystem_Subclass',
            'RAST_Subsystem_Name',
            'SEED_variant_number',
            'SEED_variant_definition',
            'SEED_curator',
            'SEED_description',
            'SEED_notes',
            'SEED_version',
            'SEED_subsystem_notes'
        ]
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(header)
        writer.writerows(all_rows)
            
    print(f"\nSuccessfully created {output_file}")
    print(f"  Total features: {len(all_rows)}")
    print(f"  - CDS: {n_cds}")
    print(f"  - RNA: {n_rna}")
    print(f"  - Prophage: {n_prophage}")


def main():
    if len(sys.argv) != 3:
        print("Usage: python consolidate_rast_output.py <rast_output_dir> <organism_name>", file=sys.stderr)
        print("\nExample:", file=sys.stderr)
        print('  python consolidate_rast_output.py output/rasttk/theta "Bacteroides thetaiotaomicron"', file=sys.stderr)
        sys.exit(1)
    
    rast_dir = sys.argv[1]
    organism_name = sys.argv[2]
    
    # Output file in same directory as input
    output_file = os.path.join(rast_dir, 'rast.tsv')
    
    consolidate_rast_annotations(rast_dir, organism_name, output_file)


if __name__ == '__main__':
    main()
