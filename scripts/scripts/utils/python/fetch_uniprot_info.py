#!/usr/bin/env python3
#
# Script Title: Fetch UniProt Information for Protein IDs
#
# Description: Retrieves UniProt annotations (e.g., gene names, recommended protein names) 
#              for a list of protein IDs. Saves the merged information to a CSV file.
#
# Author: Elad Oren
# Created: 2025-03-20
# Last Modified: 2025-04-27
#
# Usage:
#   python fetch_uniprot_info.py
#
# Input Files:
#   - protein_id_list.txt : Plain text file with one UniProt ID per line
#
# Output Files:
#   - uniprot_mapping.csv : CSV file containing retrieved UniProt annotations
#
# Dependencies:
#   - pandas
#   - requests
#
# Notes:
#   - The script queries UniProt in batches (max 500 IDs per request) to avoid overload.
#   - UniProt REST API endpoint used: https://rest.uniprot.org/uniprotkb/search
#   - Handles missing annotations gracefully (empty fields if not found).

import sys
import pandas as pd
import requests

# Usage:
# python fetch_uniprot_info.py diamond_results.tsv final_annotations.csv

diamond_file = sys.argv[1]
output_csv = sys.argv[2]

# Step 1: Parse DIAMOND output
print(f"Reading DIAMOND results from {diamond_file}...")

# Columns: qseqid, sseqid, pident, evalue, bitscore
cols = ['query_id', 'uniprot_full', 'pident', 'evalue', 'bitscore']
diamond_df = pd.read_csv(diamond_file, sep="\t", names=cols)

# Extract clean UniProt IDs from sp|ID|Name format
diamond_df['uniprot_id'] = diamond_df['uniprot_full'].apply(lambda x: x.split('|')[1] if '|' in x else x)

# Keep unique best matches per query
best_hits = diamond_df.sort_values(['query_id', 'bitscore'], ascending=[True, False]).drop_duplicates('query_id')

print(f"Found {len(best_hits)} unique best hits.")

# Step 2: Query UniProt for annotations
print("Fetching annotations from UniProt...")

def fetch_uniprot_full(uniprot_ids):
    import time
    base_url = "https://rest.uniprot.org/uniprotkb/"
    headers = {"Accept": "application/json"}
    results = []

    for uid in uniprot_ids:
        url = base_url + uid + ".json"
        print(f"Fetching {uid}...")
        try:
            r = requests.get(url, headers=headers)
            r.raise_for_status()
            data = r.json()

            # Extract protein name
            name = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', '')

            # Extract and categorize GO terms
            bp_terms = []
            mf_terms = []
            cc_terms = []
            cross_refs = data.get('uniProtKBCrossReferences', [])
            for ref in cross_refs:
                if ref.get('database') == 'GO':
                    props = {prop['key']: prop['value'] for prop in ref.get('properties', [])}
                    go_term = props.get('GoTerm', '')
                    if go_term.startswith("P:"):
                        bp_terms.append(go_term[2:])  # Remove "P:" prefix
                    elif go_term.startswith("F:"):
                        mf_terms.append(go_term[2:])
                    elif go_term.startswith("C:"):
                        cc_terms.append(go_term[2:])

            print(f"  {uid}: {len(bp_terms)} BP, {len(mf_terms)} MF, {len(cc_terms)} CC terms found.")

            results.append({
                "uniprot_id": uid,
                "Description": name,
                "BP_terms": "; ".join(bp_terms),
                "MF_terms": "; ".join(mf_terms),
                "CC_terms": "; ".join(cc_terms)
            })

        except Exception as e:
            print(f"  Error fetching {uid}: {e}")
            results.append({
                "uniprot_id": uid,
                "Description": "",
                "BP_terms": "",
                "MF_terms": "",
                "CC_terms": ""
            })
        
        time.sleep(0.1)

    return pd.DataFrame(results)

# Fetch annotations
uniprot_ids = best_hits['uniprot_id'].dropna().unique().tolist()

anno_df = fetch_uniprot_full(uniprot_ids)

print(f"Fetched annotations for {len(anno_df)} proteins.")

# Step 3: Merge based on CLEAN IDs
final_df = best_hits.merge(anno_df, on='uniprot_id', how='left')
final_df = final_df[['query_id', 'uniprot_id', 'Description', 'BP_terms', 'MF_terms', 'CC_terms']]
    
# Save final table
final_df.to_csv(output_csv, index=False)

# Debug: Show a few lines of the final table
print("\nPreview of saved file:")
print(final_df.head())

print(f"\nSaved final annotations to {output_csv}.")