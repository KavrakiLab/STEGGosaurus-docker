import os
import shutil
import sys
import random
import multiprocessing
from biopandas.pdb import PandasPdb
from mpire import WorkerPool

# Import existing modules
from numbering import *
from construction import *
from utils import *

def generate_single_conformation(idx, initial_pdb, loop_definitions, chain_letters, temp_base_dir, final_output_dir):
    """
    Worker function to generate a single TCR conformation with randomized loop order.
    
    Args:
        idx (int): The index of the conformation (0-based).
        initial_pdb (str): Path to the starting renumbered PDB (TCR0.pdb).
        loop_definitions (dict): Dictionary containing (start, end) tuples for loops a1-b3.
        chain_letters (tuple): (alpha_chain_id, beta_chain_id).
        temp_base_dir (str): The main temporary directory.
        final_output_dir (str): Directory to save the final PDB.
    """
    
    # 1. Create a unique subdirectory for this worker to prevent file conflicts
    # (Since generate_loop_conf uses fixed filenames like 'rcd_targets.txt')
    worker_dir = os.path.join(temp_base_dir, f"worker_{idx}")
    if not os.path.exists(worker_dir):
        os.makedirs(worker_dir)
        
    # 2. Copy the initial PDB to the worker directory
    current_pdb = os.path.join(worker_dir, "working_structure.pdb")
    shutil.copyfile(initial_pdb, current_pdb)
    
    alpha, beta = chain_letters
    
    # 3. Define the tasks for this conformation
    # Format: (loop_name, chain_id, start_residue, end_residue, n_samples)
    tasks = [
        ('a1', alpha, loop_definitions['a1'][0], loop_definitions['a1'][1], 200),
        ('a2', alpha, loop_definitions['a2'][0], loop_definitions['a2'][1], 200),
        ('a3', alpha, loop_definitions['a3'][0], loop_definitions['a3'][1], 500),
        ('b1', beta,  loop_definitions['b1'][0], loop_definitions['b1'][1], 200),
        ('b2', beta,  loop_definitions['b2'][0], loop_definitions['b2'][1], 200),
        ('b3', beta,  loop_definitions['b3'][0], loop_definitions['b3'][1], 500)
    ]
    
    # 4. Randomize the order of loop modeling
    random.shuffle(tasks)
    
    # 5. Process loops sequentially for this specific conformation
    for loop_name, chain, start, end, n_samples in tasks:
        fname = f"{loop_name}_{idx}"
        # generate_loop_conf takes the current PDB, modifies it, and returns the path to the new packed PDB
        current_pdb = generate_loop_conf(
            working_pdb=current_pdb,
            fname=fname,
            chain=chain,
            output_dir=worker_dir + '/', # Must end with slash for some internal concatenations
            loop_start=start,
            loop_end=end,
            n_samples=n_samples
        )
        
    # 6. Save the final result
    final_pdb_name = f"TCR{idx + 1}.pdb"
    destination = os.path.join(final_output_dir, final_pdb_name)
    shutil.copyfile(current_pdb, destination)
    
    # 7. Clean up worker directory
    shutil.rmtree(worker_dir)
    
    return final_pdb_name

if __name__ == "__main__":
    # Command-line arguments:
    # sys.argv[1]: Path to the input PDB file
    # sys.argv[2]: Number of conformations to generate
    # sys.argv[3]: (Optional) Number of cores to use. Defaults to CPU count - 1.

    if len(sys.argv) < 3:
        print("Usage: python T-RECS.py <input_pdb_file> <number_of_conformations> [num_cores]")
        sys.exit(1)

    pdb = sys.argv[1]
    num_confs = int(sys.argv[2])
    
    # Determine number of cores
    if len(sys.argv) > 3:
        num_cores = int(sys.argv[3])
    else:
        num_cores = max(1, multiprocessing.cpu_count() - 1)

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
    
    # Step 1: Renumber the PDB file (Single threaded setup)
    print("Renumbering PDB and identifying loops...")
    renumbered_pdb_path = renumber_pdb(pdb, [alpha_letter, beta_letter], output_dir)

    # Step 2: Find CDR loop regions (a1, a2, a3, b1, b2, b3) using ANARCI and IMGT scheme
    a1, a2, a3, b1, b2, b3 = find_loops(renumbered_pdb_path, [alpha_letter, beta_letter], output_dir)
    
    # Create the base template that all workers will use
    base_template_pdb = output_dir + 'TCR0.pdb'
    shutil.copyfile(renumbered_pdb_path, base_template_pdb)

    # Bundle loop definitions for the workers
    loop_definitions = {
        'a1': a1, 'a2': a2, 'a3': a3,
        'b1': b1, 'b2': b2, 'b3': b3
    }

    # Step 3: Generate loop conformations in parallel
    print(f"Generating {num_confs} TCR conformations using {num_cores} cores...")
    
    # Prepare arguments for map
    # We create a list of arguments where each element corresponds to one conformation task
    tasks_args = []
    for i in range(num_confs):
        tasks_args.append((
            i,                  # index
            base_template_pdb,  # initial structure
            loop_definitions,   # loop coordinates
            (alpha_letter, beta_letter), # chains
            output_dir,         # temp storage
            dir_to_keep         # final storage
        ))

    # Execute with mpire
    with WorkerPool(n_jobs=num_cores) as pool:
        # We assume verbose=True to see progress bars if mpire supports it in this env, 
        # otherwise a simple print helps.
        results = pool.map(generate_single_conformation, tasks_args, progress_bar=True)

    # Step 4: Cleanup
    # (Copying to dir_to_keep is done inside the worker now)
    print("Cleaning up temporary files...")
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
        
    print(f"Done! Generated {len(results)} conformations in {dir_to_keep}")
