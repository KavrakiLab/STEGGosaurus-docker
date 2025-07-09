import os
import pandas as pd
from biopandas.pdb import PandasPdb
from anarci import anarci
import numpy as np

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def decouple_chains(ref_merged,out_dir,output_filename_prefix):
    """
    Decouples chains in a merged PDB file, identifies TCR alpha and beta chains,
    renumbers residues based on ANARCI output, and filters for relevant regions.
    The input PDB is expected to have TCR alpha and beta chains merged with
    the pMHC complex, potentially resulting in chains 'D', 'E', 'C', 'A'.

    Args:
        ref_merged (str): Path to the merged PDB file containing all chains.
        out_dir (str): Directory where the decoupled and renumbered PDB
                                 file will be saved.
        output_filename_prefix (str): Prefix for the output PDB filename.
    """
    # TODO: preserve remarks in merged pdbs...

    merged_ppdb = PandasPdb().read_pdb(ref_merged)
    atom_df = merged_ppdb.df['ATOM']

    # Splitting chain 'D' (TCR) into 'D' and 'E' based on a large distance jump
    tcr_chain_d_df = atom_df[atom_df['chain_id'] == 'D']
    if not tcr_chain_d_df.empty:
        tcr_d_coords = np.array(tcr_chain_d_df[['x_coord','y_coord','z_coord']])
        tcr_d_distances = np.linalg.norm(tcr_d_coords[:-1,:] - tcr_d_coords[1:,:], axis=1)
        tcr_d_distances_padded = np.insert(tcr_d_distances, [0, len(tcr_d_distances)], 0)
        split_point_tcr = np.argmax(tcr_d_distances_padded)
        len_tcr_d = tcr_chain_d_df.shape[0]
        new_tcr_chain_ids = ['D'] * split_point_tcr + ['E'] * (len_tcr_d - split_point_tcr)
        atom_df.loc[tcr_chain_d_df.index, 'chain_id'] = new_tcr_chain_ids
    else:
        print(f"Warning: Chain 'D' (expected TCR alpha) not found in {ref_merged}.")

    # Splitting chain 'A' (pMHC) into 'C' (peptide) and 'A' (MHC)
    pmhc_chain_a_df = atom_df[atom_df['chain_id'] == 'A']
    if not pmhc_chain_a_df.empty:
        pmhc_a_coords = np.array(pmhc_chain_a_df[['x_coord','y_coord','z_coord']])
        pmhc_a_distances = np.linalg.norm(pmhc_a_coords[:-1,:] - pmhc_a_coords[1:,:], axis=1)
        pmhc_a_distances_padded = np.insert(pmhc_a_distances, [0, len(pmhc_a_distances)], 0)
        split_point_pmhc = np.argmax(pmhc_a_distances_padded)

        len_pmhc_a = len(pmhc_chain_a_df)
        new_pmhc_chain_ids = ['C'] * split_point_pmhc + ['A'] * (len_pmhc_a - split_point_pmhc)
        atom_df.loc[pmhc_chain_a_df.index, 'chain_id'] = new_pmhc_chain_ids
    else:
        print(f"Warning: Chain 'A' (expected pMHC) not found in {ref_merged}.")

    # --- ANARCI renumbering and sequence extraction for TCR chains ---
    # Extract sequences for ANARCI from the split TCR chains ('D' and 'E')
    tcr_alpha_list = atom_df[atom_df['chain_id']=='D'].drop_duplicates(subset=['residue_number'])['residue_name'].apply(lambda x: d[x]).values
    tcr_alpha_sequence = ''.join(tcr_alpha_list)

    tcr_beta_list = atom_df[atom_df['chain_id']=='E'].drop_duplicates(subset=['residue_number'])['residue_name'].apply(lambda x: d[x]).values
    tcr_beta_sequence = ''.join(tcr_beta_list)

    tcr_sequences_for_anarci = [("alpha", tcr_alpha_sequence), ("beta", tcr_beta_sequence)]

    anarci_output_filepath = os.path.join(out_dir, 'TCR_anarci.txt')
    anarci_stripped_filepath = os.path.join(out_dir, 'TCR_anarci_stripped.txt')
    anarci_clean_filepath = os.path.join(out_dir, 'TCR_anarci_stripped_clean.txt')

    anarci(tcr_sequences_for_anarci, scheme='aho', output=True, outfile=anarci_output_filepath)
    os.system(f"tr -s '[:blank:]' ',' < {anarci_output_filepath} > {anarci_stripped_filepath}")

    # Clean up ANARCI output: remove lines with too many commas (extra letters)
    with open(anarci_stripped_filepath, 'r') as infile, open(anarci_clean_filepath, 'w') as outfile:
        for line in infile:
            if line.count(',') > 2:
                line = line[:-4]+line[-2:]
            outfile.write(line)

    # Read the cleaned ANARCI output into a DataFrame
    renumbered_tcr_df = pd.read_csv(anarci_clean_filepath, names=['chain','idx','acid'], comment='#')
    # Filter out rows where amino acid is '-' (gap) or NaN
    renumbered_tcr_df = renumbered_tcr_df[renumbered_tcr_df['acid'] != '-'].dropna()

    # Extract renumbered sequences from ANARCI output for each chain
    anarci_tcr1_sequence = ''.join(renumbered_tcr_df[renumbered_tcr_df['chain']=='A']['acid'].tolist())
    anarci_tcr2_sequence = ''.join(renumbered_tcr_df[renumbered_tcr_df['chain']=='B']['acid'].tolist())

    # --- Filtering and renumbering atom DataFrame based on ANARCI results ---
    # Process TCR alpha chain ('D')
    full_sequence_alpha = tcr_sequences_for_anarci[0][1] # Original sequence used for ANARCI
    target_sequence_alpha = anarci_tcr1_sequence # Renumbered sequence from ANARCI
    chain_d_atoms = atom_df[atom_df['chain_id'] == 'D'].copy()
    if not chain_d_atoms.empty:
        # Get unique residue numbers and names to find the segment corresponding to `target_sequence_alpha`
        residues_d = chain_d_atoms.groupby(['residue_number', 'residue_name']).first().reset_index()
        try:
            start_idx_alpha = full_sequence_alpha.index(target_sequence_alpha)
            end_idx_alpha = start_idx_alpha + len(target_sequence_alpha)
            keep_residues_alpha = residues_d.iloc[start_idx_alpha:end_idx_alpha]['residue_number'].tolist()
            # Filter the atoms to keep only selected residues
            filtered_atoms_alpha = chain_d_atoms[chain_d_atoms['residue_number'].isin(keep_residues_alpha)]
            # Remove original chain D atoms and concatenate with filtered ones
            atom_df = atom_df[atom_df['chain_id'] != 'D']
            atom_df = pd.concat([atom_df, filtered_atoms_alpha])
        except ValueError:
            print(f"Warning: Target alpha sequence '{target_sequence_alpha}' not found in full alpha sequence.")
            # If not found, keep the original chain D or handle as an error

    # Process TCR beta chain ('E')
    full_sequence_beta = tcr_sequences_for_anarci[1][1] # Original sequence used for ANARCI
    target_sequence_beta = anarci_tcr2_sequence # Renumbered sequence from ANARCI
    chain_e_atoms = atom_df[atom_df['chain_id'] == 'E'].copy()
    if not chain_e_atoms.empty:
        residues_e = chain_e_atoms.groupby(['residue_number', 'residue_name']).first().reset_index()
        try:
            start_idx_beta = full_sequence_beta.index(target_sequence_beta)
            end_idx_beta = start_idx_beta + len(target_sequence_beta)
            keep_residues_beta = residues_e.iloc[start_idx_beta:end_idx_beta]['residue_number'].tolist()
            filtered_atoms_beta = chain_e_atoms[chain_e_atoms['residue_number'].isin(keep_residues_beta)]
            atom_df = atom_df[atom_df['chain_id'] != 'E']
            atom_df = pd.concat([atom_df, filtered_atoms_beta])
        except ValueError:
            print(f"Warning: Target beta sequence '{target_sequence_beta}' not found in full beta sequence.")

    # Re-map residue numbers based on ANARCI output or simply renumber sequentially
    updated_dfs = []
    for chain_id in atom_df['chain_id'].unique():
        chain_df = atom_df[atom_df['chain_id'] == chain_id].copy()

        if chain_id == 'D': # TCR Alpha chain
            anarci_alpha = renumbered_tcr_df[renumbered_tcr_df['chain']=='A'] # ANARCI 'A' corresponds to TCR alpha
            sorted_pdb_residues = sorted(chain_df['residue_number'].unique())
            if len(sorted_pdb_residues) == len(anarci_alpha['idx'].tolist()):
                residue_map = dict(zip(sorted_pdb_residues, anarci_alpha['idx'].tolist()))
            else:
                print(f"Warning: Mismatch in residue count for chain D. Renumbering sequentially.")
                residue_map = {old: new for new, old in enumerate(sorted_pdb_residues, start=1)}

        elif chain_id == 'E': # TCR Beta chain
            anarci_beta = renumbered_tcr_df[renumbered_tcr_df['chain']=='B'] # ANARCI 'B' corresponds to TCR beta
            sorted_pdb_residues = sorted(chain_df['residue_number'].unique())
            if len(sorted_pdb_residues) == len(anarci_beta['idx'].tolist()):
                residue_map = dict(zip(sorted_pdb_residues, anarci_beta['idx'].tolist()))
            else:
                print(f"Warning: Mismatch in residue count for chain E. Renumbering sequentially.")
                residue_map = {old: new for new, old in enumerate(sorted_pdb_residues, start=1)}
        else: # For other chains (e.g., peptide 'C', MHC 'A'), renumber sequentially starting from 1
            residue_map = {old: new for new, old in enumerate(sorted(chain_df['residue_number'].unique()), start=1)}

        chain_df['residue_number'] = chain_df['residue_number'].map(residue_map)
        chain_df['residue_number'] = chain_df['residue_number'].astype(int)
        chain_df = chain_df.reset_index(drop=True)
        updated_dfs.append(chain_df)

    # Concatenate all updated chain data and sort by atom number (or a logical order)
    atom_df = pd.concat(updated_dfs).sort_values(['chain_id', 'residue_number', 'atom_number']) # Sorting by atom_number might re-sort them improperly later.

    # Filter out MHC residues beyond a certain residue number ( >180)
    atom_df = atom_df.drop(atom_df[(atom_df['chain_id'] == 'A') & (atom_df['residue_number'] > 180)].index)

    # Re-index atom numbers sequentially after all filtering and renumbering
    atom_df['atom_number'] = range(1, len(atom_df) + 1)
    atom_df = atom_df.reset_index(drop=True)

    # Create a new PandasPdb object and save the modified DataFrame to a PDB file
    new_ppdb = PandasPdb()
    new_ppdb.df['ATOM'] = atom_df
    output_filepath = os.path.join(out_dir, f'{output_filename_prefix}.pdb')
    new_ppdb.to_pdb(output_filepath)