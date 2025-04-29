# Script Title: LEA3-like Motif Scanner in Hardcoded Sequences
#
# Description:
#   Scans predefined amino acid sequences for characteristic LEA3-like motifs,
#   defined as a hydrophobic stretch followed by a basic residue cluster.
#   Hydropathy scores are manually assigned and simple sliding windows are used
#   to detect motif patterns. Detected motifs and their positions are printed.
#
# Author: Elad Oren
# Created: 2025-04-27
# Last Modified: 2025-04-27
#
# Input:
#   - No external input files. Sequences are defined inside the script.
#
# Output:
#   - Motif matches are printed directly to the console.
#
# Dependencies:
#   - None (base Python only)
#
# Notes:
#   - Sequences must be manually updated within the script if new scans are desired.
#   - Designed for quick scanning of a few LEA3 candidate proteins, not for batch analysis.


hydrophobic = {'T','A','V','I','L','F','M','W','Y'}  # Apolar
positive    = {'K','R'}
negative    = {'D','E'}
def is_dure_motif(mer):
    """Return True if this 11-residue string meets the Dure motif rules:
       Positions 1,2,5,9 (0-based: 0,1,4,8) are hydrophobic,
       Positions 6,7,8 (0-based: 5,6,7) are either + - + or + Q +"""
    if len(mer) != 11:
        return False
    
    # Positions 0,1,4,8 → hydrophobic
    if (mer[0] not in hydrophobic or
        mer[1] not in hydrophobic or
        mer[4] not in hydrophobic or
        mer[8] not in hydrophobic):
        return False
    
    # Positions 5,6,7 → either (+, -, +) or (+, Q, +)
    c5, c6, c7 = mer[5], mer[6], mer[7]
    plus_minus_plus = (
        c5 in positive and
        c6 in negative and
        c7 in positive
    )
    plus_Q_plus = (
        c5 in positive and
        c6 == 'Q' and
        c7 in positive
    )
    
    return plus_minus_plus or plus_Q_plus

def find_motifs(seq):
    """Return a list of (start_index_1based, the_11mer) for every Dure motif match."""
    hits = []
    for i in range(len(seq) - 10):  # up to len(seq) - 11
        window = seq[i:i+11]
        if is_dure_motif(window):
            # store 1-based start index
            hits.append((i+1, window))
    return hits

sequences = {
  "ZmLEA3_Zm00001eb294480_P001":
    "MASHQDKASYQAGETKARTEEKTGQAVGATKDTAQHAKDRAADAAGHAAGKGQDAKEATKQKASDTGSYLGKKTDEAKHKAGETTEATKQKAGETTEAAKQKAADAMEAAKQKAAEAGQYAKDTAVSGKDKSGGVIQQATEQVKSAAAGAKDAVMSTLGMGGDDKQGDANTNKDSSTITRDH",
  "TdFL_v1.0_Td00001aa022594_T002":
    "MASHQDKASYQAGETKARTEEKTGQAVGATKDTAQHAKDRASDAAGHAEATKQRAAETAEATKQRAAETAEATKQRAAETAEATKQKAAETAEATKQKAAEYAKDTAVSGKDKSGGVIQQATEQVKSAAAGAKDAVMNTLGMGGDNKQGDTDTTKDSSTITRDH",
  "Ag00001aa069574_T002":
    "MASHQDKASYQAGETKARTEEKTGQAMGATKDTAQHAKETTKQKASDTSSYLGQKTEEAKQKAGQTTEATKQKTGETTEATKQKAGQTTEATKQKTGQTTEATKQKTGETTEAAKQKAAEAMEATKQKAAEAGQYAKETVDSGKDKSGSVIQQATEQVKSAAAGAKDAVMNTLGMGGDNNQQSDTNTNKDSSTITRDH",
  "Misin17G219700.1.v7.1":
    "MASHQDKASYQAGETKARTEEKTGQAMGATKDTAQHAKETTKQKASDTGSYLGQKTEEAKQKAGETTEAAKQKAGQTTEAAKQKAAETTEAAKQKAAEPTEEKSGGVIQQATEQVKSAAAGAKDAVMNTLGMGGDNSKQGDTNTNNSKDSSTITRDH",
  "Sn00001aa056018_T001":
    "MASHQDKASYQAGETKARTEEKTGQAMGATKDTAEHAKETTKQKASDTGSYLGQKTEEAKHKTGETTEATKQKAGETTEAAKQKTAETTEAAKQKTAETTEAAKQKAAEATEAAKETAVSGKDKSGGVIQQATEQVKSAAVGAKDAVMNTLGMGGDNNNQQQSDTNHKDSSTITRDH",
  "PvWBCH1.3NG170900.1.v1.1":
    "MASHQDKASYQAGETKARTEEKAGQAMGATKDTAQHAKDRASDAAGHAAGKGHDAKEATKQKASDTGSYLGQKTDEATKHKAGETTEATKHKAGETTEAEATKQKAGETTEAAKQKTAEATKAAKQKAAEAGEYAKESAVAGKDKTGSAIQQATEQVKSAAVGAKDAVMSTLGMSGDNKEGGAGNGKEEDHSTITRDQ"
}

for name, seq in sequences.items():
    hits = find_motifs(seq)
    if hits:
        print(f"{name} => found {len(hits)} match(es):")
        for (start, mer) in hits:
            print(f"  - starts at {start}, 11-mer = {mer}")
    else:
        print(f"{name} => no matches found.")