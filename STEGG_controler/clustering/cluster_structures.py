import os
import sys
import numpy as np
from scipy.spatial import distance_matrix
from biopandas.pdb import PandasPdb
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering

def get_irmsd_vector(pdb_filepath, atom_type='all'):
    """
    Loads atom coordinates from a PDB file and calculates a vector of
    inter-residue distances (internal RMSD components). This vector can be
    used for structural comparisons and clustering.

    Args:
        pdb_filepath (str): The path to the PDB file.
        atom_type (str): The type of atoms to consider for distance calculation.
                         'all': all heavy atoms (default).
                         'CA': only C-alpha atoms.
                         'heavy': all atoms except Hydrogen and Deuterium.

    Returns:
        numpy.ndarray: A 1D array of unique inter-atomic distances (upper triangle
                       of the distance matrix), which serves as a structural fingerprint.
    """
    ppdb = PandasPdb().read_pdb(pdb_filepath)
    atom_df = ppdb.df['ATOM']

    # Filter atoms based on the specified type
    if atom_type == 'CA':
        atom_df = atom_df[atom_df['atom_name'] == 'CA']
    elif atom_type == 'heavy':
        atom_df = atom_df[~atom_df['element_symbol'].isin(['H', 'D'])]
    # If 'all', no filtering is applied (all atoms in ATOM section are used)

    # Extract coordinates as a NumPy array
    coordinates = atom_df[['x_coord', 'y_coord', 'z_coord']].to_numpy()

    # Calculate the pairwise distance matrix
    distance_matrix_full = distance_matrix(coordinates, coordinates)

    # Extract the upper triangular part of the distance matrix (excluding diagonal)
    # This generates a unique vector of all pairwise distances without redundancy.
    row_indices, col_indices = np.triu_indices_from(distance_matrix_full, k=1)
    irmsd_vector = distance_matrix_full[row_indices, col_indices]
    return irmsd_vector

def cluster_main(pdb_directory, atom_type='all', pca_components=3, n_clusters=3):
    """
    Performs clustering on a set of PDB structures based on their internal RMSD (iRMSD) vectors.
    The process involves:
    1. Calculating iRMSD vectors for each PDB file.
    2. Reducing the dimensionality of these vectors using Principal Component Analysis (PCA).
    3. Applying Agglomerative Clustering on the PCA-reduced data.

    Args:
        pdb_directory (str): The path to the directory containing the PDB files.
        atom_type (str): The type of atoms to use for iRMSD calculation ('all', 'CA', 'heavy').
                         Defaults to 'all'.
        pca_components (int): The number of principal components to keep after PCA.
                              Defaults to 3.
        n_clusters (int): The desired number of clusters for Agglomerative Clustering.
                          Defaults to 3.

    Returns:
        dict: A dictionary mapping each PDB filename to its assigned cluster label.
    """
    # Get a list of all PDB files in the specified directory
    pdb_filenames = [f for f in os.listdir(pdb_directory) if f.endswith('.pdb')]
    pdb_filepaths = [os.path.join(pdb_directory, f) for f in pdb_filenames]

    irmsd_vectors_list = []
    valid_pdb_filenames = []

    # Iterate through PDB files, calculate iRMSD vectors, and handle potential errors
    for filepath, filename in zip(pdb_filepaths, pdb_filenames):
        try:
            vector = get_irmsd_vector(filepath, atom_type=atom_type)
            irmsd_vectors_list.append(vector)
            valid_pdb_filenames.append(filename)
        except Exception as e:
            print(f"Skipping {filename} due to error during iRMSD calculation: {e}")

    # Convert the list of vectors into a NumPy array for PCA
    if not irmsd_vectors_list:
        print("No valid PDB files found or processed for clustering.")
        return {} # Return an empty dictionary if no vectors were generated

    structure_vectors = np.array(irmsd_vectors_list)

    # Perform Principal Component Analysis (PCA) for dimensionality reduction
    pca = PCA(n_components=pca_components, svd_solver='randomized')
    reduced_vectors = pca.fit_transform(structure_vectors)

    # Perform Agglomerative Clustering
    clustering_model = AgglomerativeClustering(n_clusters=n_clusters)
    cluster_labels = clustering_model.fit_predict(reduced_vectors)

    # Map filenames to their assigned cluster labels
    cluster_assignments = dict(zip(valid_pdb_filenames, cluster_labels))

    print("\nCluster assignments:")
    for filename, label in cluster_assignments.items():
        print(f"{filename} => Cluster {label}")

    return cluster_assignments