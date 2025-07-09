import os
import shutil
import sys
from biopandas.pdb import PandasPdb

from numbering import *
from construction import *
from utils import *

if __name__ == "__main__":
    # Command-line arguments:
    # sys.argv[1]: Path to the input PDB file
    # sys.argv[2]: Number of conformations to generate

    if len(sys.argv) != 3:
        print("Usage: python T-RECS.py <input_pdb_file> <number_of_conformations>")
        sys.exit(1)

    pdb = sys.argv[1]
    num_confs = int(sys.argv[2])
    
    ppdb = PandasPdb().read_pdb(pdb)
    atom_df = ppdb.df['ATOM']
    
    pdb_chain_ids = set(atom_df['chain_id'].unique())
    if len(pdb_chain_ids) != 2:
        raise Exception("Expected a PBD file with two chains")
    if pdb_chain_ids == {'D','E'}:
        alpha_letter = 'D'
        beta_letter = 'E'
    elif pdb_chain_ids == {'A','B'}:
        alpha_letter = 'A'
        beta_letter = 'B'
    else:
        raise Exception("TCR chain letters not well defined. Please use (A,B) or (D,E)")
    
    dir_to_keep = 'output/'+pdb[:-4]+'/'
    if not os.path.exists(dir_to_keep):
        os.makedirs(dir_to_keep)
    output_dir = 'output/'+pdb[:-4]+'_temp/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Step 1: Renumber the PDB file
    renumbered_pdb_path = renumber_pdb(pdb,[alpha_letter,beta_letter],output_dir)

    # Step 2: Find CDR loop regions (a1, a2, a3, b1, b2, b3) using ANARCI and IMGT scheme
    a1,a2,a3,b1,b2,b3 = find_loops(renumbered_pdb_path,[alpha_letter,beta_letter],output_dir)
    working_pdb = output_dir+'TCR0.pdb'
    shutil.copyfile(renumbered_pdb_path,working_pdb)

    # Step 3: Generate loop conformations for each CDR loop for the specified number of conformations
    print(f"Generating {num_confs} TCR conformations...")
    for i in range(num_confs):
        print(f"Generating conformation {i+1}/{num_confs}...")
        #a1
        fname = 'a1'+str(i)
        working_pdb = generate_loop_conf(working_pdb,fname,alpha_letter,output_dir,a1[0],a1[1],200)
        #a2
        fname = 'a2'+str(i)
        working_pdb = generate_loop_conf(working_pdb,fname,alpha_letter,output_dir,a2[0],a2[1],200)
        #a3
        fname = 'a3'+str(i)
        working_pdb = generate_loop_conf(working_pdb,fname,alpha_letter,output_dir,a3[0],a3[1],500)
        #b1
        fname = 'a1'+str(i)
        working_pdb = generate_loop_conf(working_pdb,fname,beta_letter,output_dir,b1[0],b1[1],200)
        #b2
        fname = 'a1'+str(i)
        working_pdb = generate_loop_conf(working_pdb,fname,beta_letter,output_dir,b2[0],b2[1],200)
        #b3
        fname = 'a1'+str(i)
        working_pdb = generate_loop_conf(working_pdb,fname,beta_letter,output_dir,b3[0],b3[1],500)

        shutil.copyfile(working_pdb,output_dir+'TCR'+str(i+1)+'.pdb')
    
    # Step 4: Copy the generated conformations from the temporary directory to the final output directory
    for i in range(num_confs+1):
        shutil.copyfile(output_dir+'TCR'+str(i)+'.pdb',dir_to_keep+'TCR'+str(i)+'.pdb',)
    shutil.rmtree(output_dir)
