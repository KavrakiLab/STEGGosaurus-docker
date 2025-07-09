import mdtraj as md 
from pathlib import Path

def align_all_pdbs(input_dir):
    """
    Aligns all PDB files within a specified directory to the first PDB file
    in the directory, using only heavy atoms (non-hydrogen) of chain A.

    Args:
        input_dir (str): The path to the directory containing PDB files to align.

    Raises:
        ValueError: If no PDB files are found or if Chain A is not found
                    in the reference structure or if atom counts for Chain A
                    mismatch between structures.
    """
    input_path = Path(input_dir)

    # Get a sorted list of all PDB files in the directory
    pdb_files = sorted(input_path.glob("*.pdb"))
    if not pdb_files:
        raise ValueError(f"No PDB files found in {input_dir}.")

    # Load reference structure
    ref_full = md.load_pdb(str(pdb_files[0]))
    print(f"Using {pdb_files[0].name} as the reference structure for alignment.")

    # Identify heavy atoms in chain A (chain index 3)
    ref_chain_index = 3
    ref_chain_atom_indices = [
        atom.index for atom in ref_full.topology.atoms
        if atom.residue.chain.index == ref_chain_index and atom.element.symbol != 'H'
    ]

    if not ref_chain_atom_indices:
        raise ValueError("No heavy atoms found in chain 'A' of the reference structure.")

    for pdb_file in pdb_files:
        traj_full = md.load_pdb(str(pdb_file))

        # Get heavy atoms in chain A for the current structure
        traj_chain_atom_indices = [
            atom.index for atom in traj_full.topology.atoms
            if atom.residue.chain.index == ref_chain_index and atom.element.symbol != 'H'
        ]

        if len(traj_chain_atom_indices) != len(ref_chain_atom_indices):
            raise ValueError(f"Chain A heavy atom count mismatch in {pdb_file.name}")

        # Align the full structure using only heavy atoms in chain A
        traj_full.superpose(ref_full, atom_indices=traj_chain_atom_indices, ref_atom_indices=ref_chain_atom_indices)

        # Save the aligned structure
        output_file = input_path / pdb_file.name
        traj_full.save_pdb(str(output_file))
        print(f"Aligned {pdb_file.name} based on chain A heavy atoms.")
