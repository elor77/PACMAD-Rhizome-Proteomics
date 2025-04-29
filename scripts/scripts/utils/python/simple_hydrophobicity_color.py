# Script Title: Color Protein Structure by Hydrophobicity (PyMOL)
#
# Description:
#   Defines residue groups based on hydrophobicity and colors them in PyMOL.
#   Hydrophilic residues are colored blue, special residues (Cys, Lys, Arg) are
#   colored marine, and hydrophobic residues are colored red.
#   The function is registered to the PyMOL command line as 'simple_hydrophobicity'.
#
# Author: Lara Brindisi
# Created: 2025-04-27
# Last Modified: 2025-04-27
#
# Input:
#   - None (applies to the currently loaded structure in PyMOL).
#
# Output:
#   - Colors residues in the 3D view based on hydrophobicity.
#
# Dependencies:
#   - pymol.cmd (PyMOL environment required)
#
# Notes:
#   - Run 'simple_hydrophobicity' from the PyMOL command line after loading your molecule.
#   - Color scheme:
#       Hydrophilic (Ser, Thr, Asn, Gln, His, Tyr, Asp, Glu) → blue
#       Special (Cys, Lys, Arg) → marine
#       Hydrophobic (Gly, Ala, Val, Leu, Ile, Met, Phe, Trp, Pro) → red


from pymol import cmd

def simple_hydrophobicity():
    # Select groups of residues based on hydrophobicity
    cmd.select('hydrophilic', 'resn SER+THR+ASN+GLN+HIS+TYR+ASP+GLU')
    cmd.select('special', 'resn CYS+LYS+ARG')
    cmd.select('hydrophobic', 'resn GLY+ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO')

    # Apply colors to the different groups
    cmd.color('blue', 'hydrophilic')      # Hydrophilic residues in blue
    cmd.color('marine', 'special')        # Special residues in marine blue
    cmd.color('red', 'hydrophobic')       # Hydrophobic residues in red

# The following is to ensure the function works in PyMOL's environment
cmd.extend('simple_hydrophobicity', simple_hydrophobicity)