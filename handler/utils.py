from pathlib import Path
import os

def join_protein_ligand_pdbs(protein_pdb, ligand_pdb, filename):

    os.system('cat {} {} > {}'.format(protein_pdb, ligand_pdb, filename))

    return filename

    
    
