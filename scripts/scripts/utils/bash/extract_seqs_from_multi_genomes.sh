
#!/bin/bash
#
# Script Title: Extract Specific Protein Sequences from Proteomes and Optionally Rename Headers
#
# Description: Extracts protein sequences matching a list of IDs across multiple proteomes, 
#              consolidates them into a single FASTA file, and optionally renames the headers 
#              using a provided mapping file. Supports flexible input formats and automatic 
#              column detection.
#
# Author: Elad Oren
# Created: 2025-03-20
# Last Modified: 2025-04-27
#
# Input:
#   - <proteomes_dir> - Directory containing species-specific proteome FASTA files
#   - <input_file> - Tab-separated table containing at least 'protein_id' and 'proteome' columns
#   - [mapping_file] (optional) - Two-column table mapping old protein IDs to new names
#
# Output:
#   - <output_fasta> - FASTA file containing all extracted sequences
#   - <output_fasta>_renamed.aa.fa (optional) - FASTA file with headers renamed according to mapping
#
# Dependencies:
#   - seqkit (for sequence extraction and renaming)
#   - bash
#
# Notes:
#   - Columns in the input file are automatically detected by name (case-insensitive).
#   - Supports proteomes organized as separate FASTA files per species.
#   - Mapping file must have 'protein_id' and 'short_name' columns if renaming is requested.


# Usage check
if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 <proteomes_dir> <input_file> <output_fasta> [mapping_file]"
    exit 1
fi

# Assign command-line arguments
PROTEOMES_DIR="$1"
INPUT_FILE="$2"
OUTPUT_FASTA="$3"
MAPPING_FILE="$4"  # Optional mapping file

# Read the header line
HEADER=$(head -n 1 "$INPUT_FILE")

# Find the column numbers
PROTEIN_ID_COL=$(echo "$HEADER" | tr '\t' '\n' | grep -n -i "protein_id" | cut -d: -f1)
PROTEOME_COL=$(echo "$HEADER" | tr '\t' '\n' | grep -n -i "proteome" | cut -d: -f1)

# Check if columns were found
if [ -z "$PROTEIN_ID_COL" ] || [ -z "$PROTEOME_COL" ]; then
    echo "Error: Could not find 'protein_id' or 'proteome' columns in the header."
    exit 1
fi

echo "Detected columns: protein_id=$PROTEIN_ID_COL, proteome=$PROTEOME_COL"

# Remove output file if it exists
rm -f "$OUTPUT_FASTA"

# Extract unique proteomes
for proteome in $(awk -v pc="$PROTEOME_COL" 'BEGIN{FS="\t"} NR>1 {print $pc}' "$INPUT_FILE" | sort -u); do
    echo "Processing proteome: $proteome"
    
    # Get list of protein IDs for this proteome
    awk -v prot="$proteome" -v pidc="$PROTEIN_ID_COL" -v pc="$PROTEOME_COL" 'BEGIN{FS="\t"} NR>1 && $pc==prot {print $pidc}' "$INPUT_FILE" > "${proteome}.ids"
    
    # Extract sequences
    seqkit grep -f "${proteome}.ids" "${PROTEOMES_DIR}/${proteome}" >> "$OUTPUT_FASTA"
    
    # Clean up
    rm "${proteome}.ids"
done

echo "Extraction complete. Sequences combined into $OUTPUT_FASTA."

# If mapping file is provided, rename sequences
if [ ! -z "$MAPPING_FILE" ]; then
    echo "Renaming sequences using mapping file: $MAPPING_FILE"
    
    # Create a temporary file without header for seqkit
    if head -n 1 "$MAPPING_FILE" | grep -q "protein_id"; then
        echo "Detected header in mapping file, removing for seqkit compatibility"
        tail -n +2 "$MAPPING_FILE" > "${MAPPING_FILE}.tmp"
        MAPPING_FILE_FOR_SEQKIT="${MAPPING_FILE}.tmp"
    else
        MAPPING_FILE_FOR_SEQKIT="$MAPPING_FILE"
    fi
    
    # Rename the sequences using seqkit
    seqkit replace -p "^(\S+)" -r '{kv}' --kv-file "$MAPPING_FILE_FOR_SEQKIT" "$OUTPUT_FASTA" > "${OUTPUT_FASTA%.aa.fa}_renamed.aa.fa"
    
    # # Update the output file variable
    # mv "${OUTPUT_FASTA%.fasta}_renamed.fasta" "$OUTPUT_FASTA"
    
    # Clean up
    if [ -f "${MAPPING_FILE}.tmp" ]; then
        rm "${MAPPING_FILE}.tmp"
    fi
    
    echo "Sequence renaming complete."
fi

echo "All processing complete. Final sequences in $OUTPUT_FASTA."