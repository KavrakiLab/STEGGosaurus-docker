import os
import pandas as pd
from biopandas.pdb import PandasPdb
from anarci import anarci

from utils import *

def renumber_pdb(pdb,chain_letters,out_dir='output/'):
    """
    Renumber the residues in specified chains of a PDB file sequentially,
    starting from 1 for each chain.

    Args:
        pdb (str): Path to the input PDB file.
        chain_letters (list): A list of chain IDs (e.g., ['A', 'B']) whose
                                  residues should be renumbered.
        out_dir (str): Directory where the renumbered PDB will be saved.
                                 Defaults to 'output/'.

    Returns:
        str: The path to the renumbered PDB file.
    """

    ppdb = PandasPdb().read_pdb(pdb)
    atom_df = ppdb.df['ATOM']
    
    alpha_letter,beta_letter = chain_letters

    # Iterate through each specified chain and renumber its residues
    for chain in atom_df['chain_id'].unique():
        temp = atom_df[atom_df['chain_id'] == chain]
        new_nums = {}
        for i, num in enumerate(temp['residue_number'].unique()):
            new_nums[num]=i+1
        atom_df['residue_number'] = atom_df.apply(lambda x: new_nums[x.residue_number] if x.chain_id == chain else x.residue_number ,axis=1)

    atom_df['atom_number'] = list(range(1,atom_df.shape[0]+1))
    # Create a new PandasPdb object and save the modified DataFrame
    ppdb.to_pdb(out_dir+'renumbered_TCR.pdb')
    return out_dir+'renumbered_TCR.pdb'

def find_loops(pdb,chain_letters,out_dir='output/'):
    """
    Identifies the start and end residue numbers (IMGT numbering) for CDR loops
    (CDR1, CDR2, CDR3) of TCR alpha and beta chains in a PDB file.
    It uses ANARCI to get IMGT numbering.

    Args:
        pdb (str): Path to the renumbered input PDB file.
        chain_letters (list): A list of two chain IDs, where the first is
                                  assumed to be alpha and the second beta.
        out_dir (str): Directory where ANARCI output files will be saved.
                                 Defaults to 'output/'.

    Returns:
        tuple: A tuple containing six tuples, each representing the (start, end)
               residue numbers for (a1, a2, a3, b1, b2, b3) CDR loops.
    """
    ppdb = PandasPdb().read_pdb(pdb)
    atom_df = ppdb.df['ATOM']
    
    alpha_letter,beta_letter = chain_letters

    # Extract sequences for ANARCI
    TCR_A_list = atom_df[atom_df['chain_id']==alpha_letter].drop_duplicates(subset=['residue_number'])['residue_name'].apply(lambda x: d[x]).values
    TCR_A = ''
    for aa in TCR_A_list:
        TCR_A += aa

    TCR_B_list = atom_df[atom_df['chain_id']==beta_letter].drop_duplicates(subset=['residue_number'])['residue_name'].apply(lambda x: d[x]).values
    TCR_B = ''
    for aa in TCR_B_list:
        TCR_B += aa

    TCR_seq = [("alpha",TCR_A),("beta",TCR_B)]
    anarci(TCR_seq,scheme='imgt',output=True,outfile=out_dir+'TCR_anarci.txt')
    os.system("tr -s '[:blank:]' ',' < "+out_dir+"TCR_anarci.txt > "+out_dir+"TCR_anarci_stripped.txt")

    #remove cases with additional letters on rows...
    with open(out_dir+"TCR_anarci_stripped.txt", 'r') as infile, open(out_dir+"TCR_anarci_stripped_clean.txt", 'w') as outfile:
        for line in infile:
            if line.count(',') > 2:
                line = line[:-4]+line[-2:]
            outfile.write(line)

    renumbered_TCR = pd.read_csv(out_dir+"TCR_anarci_stripped_clean.txt",names=['chain','idx','acid'],comment='#')
    renumbered_TCR = renumbered_TCR[renumbered_TCR['acid'] != '-'].dropna()
    renumbered_TCR_A = renumbered_TCR[renumbered_TCR['chain']=='A'].reset_index(drop=True)
    renumbered_TCR_B = renumbered_TCR[renumbered_TCR['chain']=='B'].reset_index(drop=True)
    # ANARCI output indices are 1-based, so adjust DataFrame index to match residue numbers directly
    renumbered_TCR_A.index += 1
    renumbered_TCR_B.index += 1

    #Find the (start, end) residue numbers for each CDR loop using IMGT definitions
    a1 = (renumbered_TCR_A[renumbered_TCR_A['idx']==imgt_loops[0][0]].index[0], renumbered_TCR_A[renumbered_TCR_A['idx']==imgt_loops[0][1]].index[0])
    a2 = (renumbered_TCR_A[renumbered_TCR_A['idx']==imgt_loops[1][0]].index[0], renumbered_TCR_A[renumbered_TCR_A['idx']==imgt_loops[1][1]].index[0])
    a3 = (renumbered_TCR_A[renumbered_TCR_A['idx']==imgt_loops[2][0]].index[0], renumbered_TCR_A[renumbered_TCR_A['idx']==imgt_loops[2][1]].index[0])

    b1 = (renumbered_TCR_B[renumbered_TCR_B['idx']==imgt_loops[0][0]].index[0], renumbered_TCR_B[renumbered_TCR_B['idx']==imgt_loops[0][1]].index[0])
    b2 = (renumbered_TCR_B[renumbered_TCR_B['idx']==imgt_loops[1][0]].index[0], renumbered_TCR_B[renumbered_TCR_B['idx']==imgt_loops[1][1]].index[0])
    b3 = (renumbered_TCR_B[renumbered_TCR_B['idx']==imgt_loops[2][0]].index[0], renumbered_TCR_B[renumbered_TCR_B['idx']==imgt_loops[2][1]].index[0])

    return a1,a2,a3,b1,b2,b3
