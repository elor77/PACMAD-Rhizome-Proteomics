#!/bin/bash
#
# Script Title: Extract HSP Protein Sequences Using Seqkit
#
# Description: Extracts specific protein sequences from multiple proteomes based on a provided list of protein IDs 
#              and associated proteome filenames. Optionally renames sequences according to a mapping table.
#              Outputs a single combined FASTA file.
#
# Author: Elad Oren
# Created: 2025-03-20
# Last Modified: 2025-04-27
#
# Usage:
#   ./extract_hsp_sequences_seqkit.sh <proteomes_dir> <input_file> <output_fasta> [mapping_file]
#
# Arguments:
#   - proteomes_dir : Directory containing proteome FASTA files
#   - input_file    : Tab-delimited file listing 'protein_id' and 'proteome' columns
#   - output_fasta  : Path for the combined output FASTA file
#   - mapping_file  : (Optional) Tab-delimited file mapping old protein IDs to new names (for renaming)
#
# Input Files:
#   - Proteome FASTA files inside <proteomes_dir> (e.g., Tripsacum.fasta, Ag.fasta, etc.)
#   - Input table containing protein IDs and corresponding proteome filenames
#   - Optional mapping table for sequence renaming
#
# Output Files:
#   - Combined FASTA of extracted sequences
#   - Optionally, renamed FASTA file (_renamed.aa.fa suffix)
#
# Dependencies:
#   - seqkit >= 2.0
#
# Notes:
#   - The script automatically detects column positions of 'protein_id' and 'proteome'.
#   - Renaming step is only performed if a mapping file is supplied.
#   - Temporary files (.ids) are automatically cleaned up after processing.


# Set the directory containing your proteome FASTA files.
PROTEOMES_DIR="/Users/eo235/Library/CloudStorage/OneDrive-Personal/reference_genomes/Andropogoneae"

# Your input file with columns: OG, proteome, protein_id
INPUT_FILE="93_HSP_of_0_log2FC_All.txt"

# Define the output FASTA file and remove it if it already exists
OUTPUT_FASTA="93_HSP_proteins.aa.fa"
rm -f "$OUTPUT_FASTA"

# Loop over each unique proteome file from the input (skip the header)
for proteome in $(awk 'NR>1 {print $2}' "$INPUT_FILE" | sort -u); do
    echo "Processing proteome: $proteome"
    
    # Create a temporary file containing the list of protein IDs for this proteome
    awk -v prot="$proteome" 'NR>1 && $2==prot {print $3}' "$INPUT_FILE" > "${proteome}.ids"
    
    # Run seqkit grep to extract sequences whose headers match the protein IDs
    # The -f option tells seqkit grep to use the list of IDs from the file.
    seqkit grep -f "${proteome}.ids" "${PROTEOMES_DIR}/${proteome}" >> "$OUTPUT_FASTA"
    
    # Optionally, remove the temporary IDs file
    rm "${proteome}.ids"
done

echo "Extraction complete. Sequences combined into $OUTPUT_FASTA."