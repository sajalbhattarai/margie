#!/usr/bin/env python3
"""
Consolidate RASTtk annotation outputs into a single rast.tsv file.

This script merges CDS, RNA, and prophage annotations with sequences,
extracting EC numbers and parsing gene locations.
It also enriches the output with Subsystem information if a mapping file is available.

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

# -----------------------------------------------------------------------------
# Subsystem Enrichment Logic
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
    Consolidate all RASTtk outputs into a single rast.tsv file with enrichment.
    """
    rast_path = Path(rast_dir)
    
    # -------------------------------------------------------------------------
    # 1. Setup Enrichment if available
    # -------------------------------------------------------------------------
    # Assuming mapping file is at ../../../processing/database/subsystem_mapping.tsv relative to this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Adjust path: script is in processing/scripts/rasttk
    # Database is in processing/database
    base_dir = os.path.abspath(os.path.join(script_dir, '../../../'))
    # Also check /container/db/subsystem_mapping.tsv for CONTAINER environment
    
    mapping_path_local = os.path.join(base_dir, "processing/database/subsystem_mapping.tsv")
    mapping_path_container = "/container/db/subsystem_mapping.tsv"
    
    mapping_file = None
    if os.path.exists(mapping_path_container):
        mapping_file = mapping_path_container
    elif os.path.exists(mapping_path_local):
        mapping_file = mapping_path_local
    
    role_to_subsystems = {}
    seed_variant_status = "pending"
    seed_variant_message = "SEED variant enrichment is coming soon; placeholders are emitted and pipeline continues."
    if mapping_file:
        role_to_subsystems = load_subsystem_mapping(mapping_file)
    else:
        print("No subsystem mapping file found. Proceeding without enrichment.")

    # SEED variant enrichment placeholders are intentionally emitted for now.
    print(f"[INFO] {seed_variant_message}")

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
            # Subsystem Enrichment
            # --------------------------------
            sub_superclass = []
            sub_class = []
            sub_subclass = []
            sub_name = []
            
            # Lookup role (description) in dictionary
            # Note: RAST descriptions might vary slightly, exact match required here
            if description in role_to_subsystems:
                for sub in role_to_subsystems[description]:
                    if sub['superclass']: sub_superclass.append(sub['superclass'])
                    if sub['class']: sub_class.append(sub['class'])
                    if sub['subclass']: sub_subclass.append(sub['subclass'])
                    if sub['name']: sub_name.append(sub['name'])
            
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
                # New Columns with RAST_ prefix
                join_unique(sub_superclass),
                join_unique(sub_class),
                join_unique(sub_subclass),
                join_unique(sub_name),
                # Placeholder SEED variant fields (non-blocking until full variant integration)
                'NA',
                'COMING_SOON',
                'COMING_SOON',
                seed_variant_message,
                '',
                '',
                seed_variant_status
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
        # Write header
        header = [
            'organism',
            'feature_id',
            'RAST_gene_id',
            'gene_start',
            'gene_end',
            'na_length',
            'aa_length',
            'na_seq',
            'aa_seq',
            'RAST_gene_type',
            'RAST_description',
            'RAST_EC_numbers',
            'RAST_strand',
            'RAST_feature_hash',
            'RAST_source_file',
            'RAST_genecaller_raw_row',
            # New Headers with RAST_ prefix
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
            'SEED_variant_status'
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
