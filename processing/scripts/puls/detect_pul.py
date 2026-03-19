#!/usr/bin/env python3
"""
Detect Polysaccharide Utilization Loci (PULs)
"""

import argparse
import pandas as pd
import re

def detect_pul(consolidated_file, output_tsv, output_html=None):
    """
    Identify PULs based on:
    - CAZyme clustering (≥3 CAZymes within 10kb)
    - SusCD systems (SusC/SusD pairs)
    """

    annotations = pd.read_csv(consolidated_file, sep='\t')

    if 'feature_id' not in annotations.columns:
        raise ValueError("Consolidated file is missing required column: feature_id")

    if 'gene_id' in annotations.columns:
        contig_series = annotations['gene_id'].astype(str)
    elif 'RAST_gene_id' in annotations.columns:
        contig_series = annotations['RAST_gene_id'].astype(str)
    else:
        contig_series = pd.Series(['unknown_contig'] * len(annotations))

    if 'RASTTK_gene_start' in annotations.columns and 'RASTTK_gene_end' in annotations.columns:
        start_series = pd.to_numeric(annotations['RASTTK_gene_start'], errors='coerce').fillna(0).astype(int)
        end_series = pd.to_numeric(annotations['RASTTK_gene_end'], errors='coerce').fillna(0).astype(int)
    elif 'gene_start' in annotations.columns and 'gene_end' in annotations.columns:
        start_series = pd.to_numeric(annotations['gene_start'], errors='coerce').fillna(0).astype(int)
        end_series = pd.to_numeric(annotations['gene_end'], errors='coerce').fillna(0).astype(int)
    else:
        start_series = pd.Series([0] * len(annotations))
        end_series = pd.Series([0] * len(annotations))

    merged = annotations.copy()
    merged['id'] = merged['feature_id'].astype(str)
    merged['contig'] = contig_series
    merged['start'] = start_series
    merged['end'] = end_series

    if 'RASTTK_strand' in merged.columns:
        merged['strand'] = merged['RASTTK_strand']
    elif 'RAST_strand' in merged.columns:
        merged['strand'] = merged['RAST_strand']
    else:
        merged['strand'] = '.'

    dbcan_count = pd.to_numeric(merged.get('DBCAN_family_count', 0), errors='coerce').fillna(0)
    dbcan_families = merged.get('DBCAN_families', pd.Series([''] * len(merged)))
    
    # Properly identify CAZymes: must have family_count > 0
    merged['is_cazyme'] = (dbcan_count > 0)

    desc_cols = [
        'RAST_description',
        'UNIPROT_protein_name',
        'KEGG_KO_definitions',
        'TCDB_description',
        'MEROPS_peptidase_name'
    ]
    existing_desc_cols = [col for col in desc_cols if col in merged.columns]
    if existing_desc_cols:
        merged['combined_description'] = merged[existing_desc_cols].fillna('').astype(str).agg(' '.join, axis=1)
    else:
        merged['combined_description'] = ''

    merged['is_susc'] = merged['combined_description'].str.contains(r'\bSusC\b|TonB', case=False, na=False)
    merged['is_susd'] = merged['combined_description'].str.contains(r'\bSusD\b', case=False, na=False)
    
    # Find PUL clusters
    puls = []
    max_distance = 10000  # bp
    min_cazymes = 3
    
    for contig in merged['contig'].dropna().unique():
        contig_genes = merged[merged['contig'] == contig].sort_values('start')
        if contig_genes.empty:
            continue
        
        i = 0
        while i < len(contig_genes):
            if contig_genes.iloc[i]['is_cazyme']:
                # Start potential PUL
                pul_start = i
                pul_genes = [i]
                cazyme_count = 1
                
                # Extend PUL
                j = i + 1
                while j < len(contig_genes):
                    if contig_genes.iloc[j]['start'] - contig_genes.iloc[j-1]['end'] > max_distance:
                        break
                    
                    if contig_genes.iloc[j]['is_cazyme']:
                        cazyme_count += 1
                    
                    pul_genes.append(j)
                    j += 1
                
                # Check if valid PUL
                if cazyme_count >= min_cazymes:
                    pul_data = contig_genes.iloc[pul_genes]
                    puls.append({
                        'pul_id': f"PUL_{contig}_{pul_start}",
                        'contig': contig,
                        'start': pul_data['start'].min(),
                        'end': pul_data['end'].max(),
                        'num_genes': len(pul_genes),
                        'num_cazymes': cazyme_count,
                        'has_suscd': pul_data['is_susc'].any() and pul_data['is_susd'].any(),
                        'genes': ','.join(pul_data['id'].tolist())
                    })
                
                i = j
            else:
                i += 1
    
    # Save results
    pul_df = pd.DataFrame(puls)
    pul_df.to_csv(output_tsv, sep='\t', index=False)

    if output_html:
        generate_pul_viz(pul_df, merged, output_html)

def generate_pul_viz(pul_df, genes_df, output_html):
    """Generate HTML visualization of PULs"""
    
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>PUL Detection Report</title>
        <style>
            body {{ font-family: Arial; margin: 20px; }}
            .pul-box {{ border: 2px solid #4CAF50; padding: 10px; margin: 10px 0; }}
            .gene {{ display: inline-block; padding: 5px; margin: 2px; border: 1px solid #ddd; }}
            .cazyme {{ background-color: #ffeb3b; }}
            .suscd {{ background-color: #2196f3; color: white; }}
        </style>
    </head>
    <body>
        <h1>Polysaccharide Utilization Loci (PULs)</h1>
        <p>Detected {len(pul_df)} PULs</p>
        
        <h2>PUL Details</h2>
        {pul_df.to_html(index=False)}
    </body>
    </html>
    """
    
    with open(output_html, 'w') as f:
        f.write(html)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Detect Polysaccharide Utilization Loci (PULs)')
    parser.add_argument('--annotations', required=True, help='Consolidated annotation TSV file')
    parser.add_argument('--output', required=True, help='Output PUL TSV file')
    parser.add_argument('--html', required=False, help='Optional HTML report output path')
    args = parser.parse_args()

    detect_pul(args.annotations, args.output, args.html)
