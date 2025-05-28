## script to be executed in PyMol
# This script will read a list of accessions and look for the corresponding PDB files
# Then, it will read a TSV and color the residues in a given range from column X in the TSV

import os
import pandas as pd
from pymol import cmd
from pymol import util
os.chdir('/Users/estebanlopeztavera/Library/CloudStorage/OneDrive-NorwegianUniversityofLifeSciences/Collab/02-Zarah/2023-08-CBMs_alignment_analysis')  # Change to the directory where the PDB files are located

# Read list of accessions
accession_list = []
listpath = '2025-04-14/Data/03-find_motifs/accessions_MMH_noAA10.txt'  # Path to the file containing the list of accessions
with open(listpath, 'r') as f:
    for line in f:
        accession = line.strip()
        if accession:
            accession_list.append(accession)
print(f"Number of accessions: {len(accession_list)}")

# Read the TSV file
tsv_path = '2025-04-14/Data/06-structure_analysis/CBM2_interpro_ranges.tsv'  # Path to the TSV file
df = pd.read_csv(tsv_path, sep='\t').set_index('uniprot')

# load the PDB files
pdbdir = '2025-04-14/Data/05-structures'  # Directory where the PDB files are located
pdbs_list = os.listdir(pdbdir)
loadpdbs = [pdb for pdb in pdbs_list if pdb.endswith('.pdb') and pdb.split('_')[0] in accession_list]
for pdb in loadpdbs:
    cmd.load(f'{pdbdir}/{pdb}')
print(f"Number of PDB files loaded: {len(loadpdbs)}")

# Color the residues in the given range
cmd.spectrum('b', 'orange_yellow_cyan_blue', minimum=50)
for pdb in loadpdbs:
    # Get the accession number from the PDB filename
    accession = pdb.split('_')[0]
    try: 
        cmd.color('white', f'{pdb.replace(".pdb","")} and resi {df.loc[accession, "CBM2"]}')
        # show methionines and histidines in the cbm2
        cmd.show('sticks', f'{pdb.replace(".pdb","")} and resi {df.loc[accession, "CBM2"]} and (resn MET or resn HIS)')
        cmd.color('green', f'{pdb.replace(".pdb","")} and resi {df.loc[accession, "CBM2"]} and (resn MET or resn HIS)')
        util.cnc('all')
    except Exception as e:
        print(f"Error coloring {pdb}: {e}")
        continue

# superimpose the structures
target = loadpdbs[-1].replace('.pdb', '')
for pdb in loadpdbs[:-1]:
    cmd.super(pdb.replace('.pdb', ''), target)

