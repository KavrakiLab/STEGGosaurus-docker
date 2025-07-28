import pandas as pd
from biopandas.pdb import PandasPdb
from anarci import anarci
import os
import shutil
import sys

from utils import *

def prepare_pdb_files(pdb_filepath, molecule_type, output_directory='output/'):
    """
    Prepares PDB files by selecting chains, removing heteroatoms, fixing insertions,
    selecting alternate locations, keeping coordinates, and tidying the structure.
    Specifically handles TCR (T-cell Receptor) and pMHC (peptide-MHC) molecules.

    Args:
        pdb_filepath (str): Path to the input PDB file.
        molecule_type (str): Type of molecule ('TCR' or 'pMHC').
        output_directory (str): Directory where cleaned PDB files will be saved.
                                 Defaults to 'output/'.
    """
    # Ensure the output directory exists
    os.makedirs(output_directory, exist_ok=True)

    if molecule_type == 'TCR':
        # Define output filenames for TCR chains and the merged clean TCR
        tcr_chain0_path = os.path.join(output_directory, 'TCR_chain0.pdb')
        tcr_chain1_path = os.path.join(output_directory, 'TCR_chain1.pdb')
        tcr_clean_path = os.path.join(output_directory, 'TCR_clean.pdb')

        # Clean up previous runs' temporary files if they exist
        for f in [tcr_chain0_path, tcr_chain1_path, tcr_clean_path]:
            if os.path.exists(f):
                os.remove(f)

        # Process TCR chains: A and B
        os.system(f'pdb_selchain -A {pdb_filepath} | pdb_delhetatm | pdb_fixinsert | pdb_selaltloc | pdb_keepcoord | pdb_tidy -strict > {tcr_chain0_path}')
        os.system(f'pdb_selchain -B {pdb_filepath} | pdb_delhetatm | pdb_fixinsert | pdb_selaltloc | pdb_keepcoord | pdb_tidy -strict > {tcr_chain1_path}')

        # Merge TCR chains, renumber residues, change chain IDs, and tidy
        os.system(f'pdb_merge {tcr_chain0_path} {tcr_chain1_path} | pdb_reres -1 | pdb_chain -D | pdb_chainxseg | pdb_tidy -strict > {tcr_clean_path}')

    elif molecule_type == 'pMHC':
        # Define output filenames for peptide, MHC chain, and the merged clean pMHC
        peptide_chain_path = os.path.join(output_directory, 'pep_chain.pdb')
        mhc_chain_path = os.path.join(output_directory, 'MHC_chain.pdb')
        pmhc_clean_path = os.path.join(output_directory, 'pMHC_clean.pdb')

        # Clean up previous runs' temporary files if they exist
        for f in [peptide_chain_path, mhc_chain_path, pmhc_clean_path]:
            if os.path.exists(f):
                os.remove(f)

        # Process pMHC chains: C (peptide) and A (MHC)
        os.system(f'pdb_selchain -C {pdb_filepath} | pdb_delhetatm | pdb_fixinsert | pdb_selaltloc | pdb_keepcoord | pdb_tidy -strict > {peptide_chain_path}')
        os.system(f'pdb_selchain -A {pdb_filepath} | pdb_delhetatm | pdb_fixinsert | pdb_selaltloc | pdb_keepcoord | pdb_tidy -strict > {mhc_chain_path}')

        # Merge peptide and MHC chains, renumber residues, change chain IDs, and tidy
        os.system(f'pdb_merge {peptide_chain_path} {mhc_chain_path} | pdb_reres -1 | pdb_chain -A | pdb_chainxseg | pdb_tidy -strict > {pmhc_clean_path}')
    else:
        raise ValueError("Invalid molecule_type. Must be 'TCR' or 'pMHC'.")


def generate_haddock_air_file(tcr_pdb_path, pmhc_pdb_path, sequence_dictionary, output_directory='output/'):
    """
    Generates an Active Interface Residues (AIR) restraint file for HADDOCK docking.
    This file specifies distance restraints between CDR3 loops of the TCR and the
    peptide of the pMHC.

    Args:
        tcr_pdb_path (str): Path to the cleaned TCR PDB file.
        pmhc_pdb_path (str): Path to the cleaned pMHC PDB file.
        sequence_dictionary (dict): A dictionary containing the sequences of 'a3' (alpha CDR3),
                                    'b3' (beta CDR3), and 'p' (peptide).
        output_directory (str): Directory where the restraint files will be saved.
                                 Defaults to 'output/'.
    """
    # Ensure the output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # Read TCR PDB and extract atom information
    tcr_ppdb = PandasPdb().read_pdb(tcr_pdb_path)
    tcr_atom_df = tcr_ppdb.df['ATOM']

    # Read pMHC PDB and extract atom information
    pmhc_ppdb = PandasPdb().read_pdb(pmhc_pdb_path)
    pmhc_atom_df = pmhc_ppdb.df['ATOM']

    # Extract TCR sequence from residue names, converting 3-letter codes to 1-letter codes
    tcr_sequence_list = tcr_atom_df.drop_duplicates(subset=['residue_number'])['residue_name'].apply(lambda x: d[x]).values
    tcr_sequence = ''.join(tcr_sequence_list)

    # Extract pMHC sequence from residue names, converting 3-letter codes to 1-letter codes
    pmhc_sequence_list = pmhc_atom_df.drop_duplicates(subset=['residue_number'])['residue_name'].apply(lambda x: d[x]).values
    pmhc_sequence = ''.join(pmhc_sequence_list)

    # Find the start and end indices of CDR3 alpha (a3), CDR3 beta (b3) in TCR sequence,
    # and peptide (p) in pMHC sequence.
    cdr3a_indices = [tcr_sequence.find(sequence_dictionary['a3']), tcr_sequence.find(sequence_dictionary['a3']) + len(sequence_dictionary['a3'])]
    cdr3b_indices = [tcr_sequence.find(sequence_dictionary['b3']), tcr_sequence.find(sequence_dictionary['b3']) + len(sequence_dictionary['b3'])]
    peptide_indices = [pmhc_sequence.find(sequence_dictionary['p']), pmhc_sequence.find(sequence_dictionary['p']) + len(sequence_dictionary['p'])]

    # Generate the HADDOCK AIR restraint file
    air_restraints_filepath = os.path.join(output_directory, 'tcr_phla_restraints.tbl')
    with open(air_restraints_filepath, 'w') as outfile:
        outfile.write('! HADDOCK AIR restraints\n')
        outfile.write('! generated automatically by STEGG\n')
        outfile.write('!\n')
        outfile.write('! alpha restraints\n')
        # Iterate through residues of CDR3 alpha and define restraints with peptide residues
        for i in range(cdr3a_indices[0], cdr3a_indices[1] + 1):
            outfile.write(f'assign ( resid {i} and segid D)\n') # Assuming TCR alpha chain is 'D'
            outfile.write('(\n')
            for j in range(peptide_indices[0], peptide_indices[1]):
                outfile.write(f'( resid {j} and segid A)\n') # Assuming MHC and peptide chain is 'A'
                outfile.write('or\n')
            outfile.write(f'( resid {peptide_indices[1]} and segid A)\n')
            outfile.write(') 2.0 2.0 0.0\n\n') # Distance restraint parameters (distance, lower bound, upper bound)

        outfile.write('!\n')
        outfile.write('! beta restraints\n')
        # Iterate through residues of CDR3 beta and define restraints with peptide residues
        for i in range(cdr3b_indices[0], cdr3b_indices[1] + 1):
            outfile.write(f'assign ( resid {i} and segid D)\n') # Assuming TCR beta chain is also 'D' for now
            outfile.write('(\n')
            for j in range(peptide_indices[0], peptide_indices[1]):
                outfile.write(f'( resid {j} and segid A)\n')
                outfile.write('or\n')
            outfile.write(f'( resid {peptide_indices[1]} and segid A)\n')
            outfile.write(') 2.0 2.0 0.0\n\n')

    # Generate HADDOCK body restraints for TCR and pMHC
    chain_restraint0_path = os.path.join(output_directory, 'chain_restraint0.tbl')
    chain_restraint1_path = os.path.join(output_directory, 'chain_restraint1.tbl')
    combined_chain_restraint_path = os.path.join(output_directory, 'chain_restraint.tbl')

    os.system(f'haddock3-restraints restrain_bodies {tcr_pdb_path} > {chain_restraint0_path}')
    os.system(f'haddock3-restraints restrain_bodies {pmhc_pdb_path} > {chain_restraint1_path}')

    # Combine the two chain restraint files into a single file
    with open(combined_chain_restraint_path, 'w') as outfile:
        for fname in [chain_restraint0_path, chain_restraint1_path]:
            with open(fname, 'r') as infile:
                shutil.copyfileobj(infile, outfile) # Efficiently copy file content

def generate_cfg(sample_id,template='template_config.txt'):

    with open(template, 'r') as f:
        config_template = f.read()

    config_modified = config_template.replace('?', sample_id)

    with open(sample_id+'.cfg', 'w') as f:
        f.write(config_modified)

    print(f"Config written to: "+sample_id+'.cfg')

if __name__ == "__main__":
    # Command-line arguments:
    # sys.argv[1]: Directory containing TCR PDB files
    # sys.argv[2]: Directory containing pMHC PDB files
    # sys.argv[3]: CDR3 alpha sequence
    # sys.argv[4]: CDR3 beta sequence
    # sys.argv[5]: Peptide sequence

    if len(sys.argv) != 6:
        print("Usage: python dock_TCR_pHLA.py <TCR_PDB_dir> <PMHC_PDB_dir> <CDR3a_sequence> <CDR3b_sequence> <Peptide_sequence>")
        sys.exit(1)

    tcr_pdb_directory = sys.argv[1]
    pmhc_pdb_directory = sys.argv[2]
    cdr3a_sequence = sys.argv[3]
    cdr3b_sequence = sys.argv[4]
    peptide_sequence = sys.argv[5]

    # Store sequences in a dictionary for easy access
    sequences = {'a3': cdr3a_sequence, 'b3': cdr3b_sequence, 'p': peptide_sequence}

    # Get list of PDB files from the input directories
    tcr_pdb_files = [f for f in os.listdir(tcr_pdb_directory) if f.endswith('.pdb')]
    pmhc_pdb_files = [f for f in os.listdir(pmhc_pdb_directory) if f.endswith('.pdb')]

    cleaned_tcr_paths = []
    cleaned_pmhc_paths = []
    output_base_directory = 'output/'+tcr_pdb_directory.rsplit('/', 1)[-1]+'_'+pmhc_pdb_directory.rsplit('/', 1)[-1] # Define base output directory

    # delete base output directory and remake it
    if os.path.exists(output_base_directory) and os.path.isdir(output_base_directory):
        shutil.rmtree(output_base_directory)
    os.makedirs(output_base_directory, exist_ok=True)

    # Prepare TCR PDB files
    for i, tcr_filename in enumerate(tcr_pdb_files):
        full_tcr_path = os.path.join(tcr_pdb_directory, tcr_filename)
        prepare_pdb_files(full_tcr_path, 'TCR', output_base_directory)
        # Rename the cleaned TCR file to a unique name and store its path
        temp_tcr_name = f'{tcr_filename[:-4]}_{i}.pdb'
        final_tcr_clean_path = os.path.join(output_base_directory, temp_tcr_name)
        shutil.move(os.path.join(output_base_directory, 'TCR_clean.pdb'), final_tcr_clean_path)
        cleaned_tcr_paths.append(final_tcr_clean_path)

    # Prepare pMHC PDB files
    for i, pmhc_filename in enumerate(pmhc_pdb_files):
        full_pmhc_path = os.path.join(pmhc_pdb_directory, pmhc_filename)
        prepare_pdb_files(full_pmhc_path, 'pMHC', output_base_directory)
        # Rename the cleaned pMHC file to a unique name and store its path
        temp_pmhc_name = f'{pmhc_filename[:-4]}_{i}.pdb'
        final_pmhc_clean_path = os.path.join(output_base_directory, temp_pmhc_name)
        shutil.move(os.path.join(output_base_directory, 'pMHC_clean.pdb'), final_pmhc_clean_path)
        cleaned_pmhc_paths.append(final_pmhc_clean_path)

    # Create ensembles of cleaned TCR and pMHC structures
    # pdb_mkensemble: combine multiple PDB files into a single ensemble PDB file
    tcr_ensemble_path = os.path.join(output_base_directory, 'TCR_ensemble.pdb')
    pmhc_ensemble_path = os.path.join(output_base_directory, 'pMHC_ensemble.pdb')

    tcr_ensemble_command = f'pdb_mkensemble {" ".join(cleaned_tcr_paths)} > {tcr_ensemble_path}'
    os.system(tcr_ensemble_command)

    pmhc_ensemble_command = f'pdb_mkensemble {" ".join(cleaned_pmhc_paths)} > {pmhc_ensemble_path}'
    os.system(pmhc_ensemble_command)

    # Generate HADDOCK AIR restraints using the first cleaned TCR and pMHC (representative structures)
    if cleaned_tcr_paths and cleaned_pmhc_paths:
        generate_haddock_air_file(cleaned_tcr_paths[0], cleaned_pmhc_paths[0], sequences, output_base_directory)
    else:
        print("Error: No cleaned TCR or pMHC files available for restraint generation.")
        sys.exit(1)

    generate_cfg(tcr_pdb_directory.rsplit('/', 1)[-1]+'_'+pmhc_pdb_directory.rsplit('/', 1)[-1])
    os.system('haddock3 '+tcr_pdb_directory.rsplit('/', 1)[-1]+'_'+pmhc_pdb_directory.rsplit('/', 1)[-1]+'.cfg')