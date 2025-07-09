import os
import shutil
import json
import pandas as pd
import numpy as np
from biopandas.pdb import PandasPdb
from anarci import anarci
from Bio import PDB
from Bio.Data.IUPACData import protein_letters_3to1

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def safe_three_to_one(resname):
    """Convert 3-letter residue to 1-letter. Returns 'X' if unknown."""
    try:
        return protein_letters_3to1[resname.capitalize()]
    except KeyError:
        return 'X'

def extract_chain_sequence(structure, chain_id):
    """Extract full sequence from a given chain."""
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                seq = ''
                for residue in chain:
                    if PDB.is_aa(residue, standard=True):
                        seq += safe_three_to_one(residue.get_resname())
                return seq
    return ''

def extract_residue_range_sequence(structure, chain_id, start, end):
    """Extract sequence for residues in a specified number range in a chain."""
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                residues = [residue for residue in chain if start <= residue.id[1] <= end]
                return ''.join(
                    safe_three_to_one(res.get_resname()) if PDB.is_aa(res, standard=True) else 'X'
                    for res in residues
                )
    return ''

def process_pdb_files_to_json(input_dir):
    parser = PDB.PDBParser(QUIET=True)

    for filename in os.listdir(input_dir):
        if filename.endswith(".pdb"):
            base_name = os.path.splitext(filename)[0]
            pdb_path = os.path.join(input_dir, filename)
            output_dir = os.path.join(input_dir, base_name)

            os.makedirs(output_dir, exist_ok=True)
            shutil.copy(pdb_path, output_dir)

            structure = parser.get_structure(base_name, pdb_path)

            data = {
                "MHC": extract_chain_sequence(structure, 'A'),
                "B2M": extract_chain_sequence(structure, 'B'),
                "peptide": extract_chain_sequence(structure, 'C'),
                "alpha": extract_chain_sequence(structure, 'D'),
                "beta": extract_chain_sequence(structure, 'E'),
                "CDR3a": extract_residue_range_sequence(structure, 'D', 105, 116),
                "CDR3b": extract_residue_range_sequence(structure, 'E', 105, 116)
            }

            with open(os.path.join(output_dir, f"{base_name}.json"), 'w') as json_file:
                json.dump(data, json_file, indent=4)

def get_ROI(input_pdb,out_dir):
    ppdb = PandasPdb().read_pdb(input_pdb)
    atom_df = ppdb.df['ATOM']

    #drop B2M
    atom_df = atom_df[atom_df['chain_id'] != 'B']

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
    output_filepath = os.path.join(out_dir, out_dir+'_trimmed.pdb')
    new_ppdb.to_pdb(output_filepath)

def trim_pdbs_to_ROI(input_dir):
    for dirname in os.listdir(input_dir):
        if os.path.isdir(os.path.join(input_dir, dirname)):
            print('trimming ',dirname)
            filename = input_dir+'/'+dirname+'/'+dirname+'.pdb'
            get_ROI(filename,os.path.join(input_dir, dirname))


# Example usage:
# process_pdb_files_to_json("imgt")
trim_pdbs_to_ROI('/home/jared/Downloads/STEGG_benchmark_dataset-selected/processed_structures/')

