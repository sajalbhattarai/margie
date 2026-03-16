#!/usr/bin/env python3
"""
Parse Pfam-A.hmm.dat file to extract family metadata
Converts Stockholm format to simple TSV
"""

import sys
import re

def parse_pfam_dat(dat_file, output_file):
    """
    Parse Pfam-A.hmm.dat Stockholm format file.
    
    Fields to extract:
    - AC: Accession (PF00001.23)
    - ID: Family name (ABC_tran)
    - DE: Description
    - CL: Clan (optional)
    - TP: Type (Domain/Family/Repeat)
    - ML: Model length
    """
    
    entries = []
    current_entry = {}
    
    with open(dat_file) as f:
        for line in f:
            line = line.rstrip()
            
            # Start of new entry
            if line == "# STOCKHOLM 1.0":
                if current_entry:
                    entries.append(current_entry)
                current_entry = {}
            
            # Parse metadata fields
            elif line.startswith('#=GF '):
                parts = line.split(None, 3)
                if len(parts) >= 3:
                    field = parts[1]
                    value = parts[2] if len(parts) == 3 else ' '.join(parts[2:])
                    
                    if field == 'AC':
                        current_entry['accession'] = value
                    elif field == 'ID':
                        current_entry['family_name'] = value
                    elif field == 'DE':
                        current_entry['description'] = value
                    elif field == 'CL':
                        current_entry['clan'] = value
                    elif field == 'TP':
                        current_entry['type'] = value
                    elif field == 'ML':
                        current_entry['model_length'] = value
            
            # End of entry
            elif line == "//":
                if current_entry:
                    entries.append(current_entry)
                    current_entry = {}
    
    # Don't forget last entry
    if current_entry:
        entries.append(current_entry)
    
    # Write TSV output
    print(f"Writing {len(entries)} Pfam family entries to {output_file}...")
    
    with open(output_file, 'w') as out:
        # Header
        out.write('\t'.join([
            'accession',
            'family_name',
            'description',
            'clan',
            'type',
            'model_length'
        ]) + '\n')
        
        # Entries
        for entry in entries:
            out.write('\t'.join([
                entry.get('accession', ''),
                entry.get('family_name', ''),
                entry.get('description', ''),
                entry.get('clan', ''),
                entry.get('type', ''),
                entry.get('model_length', '')
            ]) + '\n')
    
    print(f"✓ Parsed {len(entries)} Pfam families")
    
    # Print some stats
    with_clan = sum(1 for e in entries if e.get('clan'))
    types = {}
    for e in entries:
        t = e.get('type', 'Unknown')
        types[t] = types.get(t, 0) + 1
    
    print(f"\nStatistics:")
    print(f"  Total families: {len(entries)}")
    print(f"  Families with clan: {with_clan}")
    print(f"  Types: {dict(types)}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: parse_pfam_metadata.py <Pfam-A.hmm.dat> <output.tsv>")
        sys.exit(1)
    
    parse_pfam_dat(sys.argv[1], sys.argv[2])
