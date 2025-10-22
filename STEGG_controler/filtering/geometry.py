
from glob import glob
import itertools
import os
from typing import Dict, List, Tuple

import numpy as np
from numpy.typing import NDArray
import pandas as pd
from sklearn.decomposition import PCA

from biopandas.pdb import PandasPdb

from .criteria import CUTOFFS, DIRECTIONS, APPLIED


def CoM(pdb_path: str, chains: List[str]) -> np.ndarray:

    # Load pdb
    ppdb = PandasPdb().read_pdb(pdb_path)

    # Extract only CA atoms
    df = ppdb.df['ATOM']
    df = df[df['atom_name'] == 'CA']
    df = df[df['chain_id'].isin(chains)]
    coords = df[['x_coord', 'y_coord', 'z_coord']].to_numpy()
    #get CoM
    com = np.mean(coords, axis=0)

    return com


def orient_vector(vector: NDArray, direction_vector: NDArray) -> NDArray:

    # Calculate the dot product of the direction vector and the vector
    dot_product = np.dot(direction_vector, vector)

    # Change the direction of the vector based on direction_vector
    oriented_vector = vector if dot_product >= 0 else -vector

    return oriented_vector


def get_end_ref(pdb_path: str, chains: List[str]) -> NDArray:

    # Load pdb
    ppdb = PandasPdb().read_pdb(pdb_path)

    # Extract only CA atoms
    df = ppdb.df['ATOM']
    df = df[df['atom_name'] == 'CA']
    df = df[df['chain_id'].isin(chains)]
    coords = df[['x_coord', 'y_coord', 'z_coord']].to_numpy()

    #return last CA atom
    return coords[-1,:]


def get_MHC_vectors(pdb_path: str, chains: List[str] = ['A']) -> NDArray:

    # get elements we need
    v_mhc = get_end_ref(pdb_path,['C']) - CoM(pdb_path,['A','C'])

    # Load pdb
    ppdb = PandasPdb().read_pdb(pdb_path)

    # Extract only CA atoms
    df = ppdb.df['ATOM']
    df = df[df['atom_name'] == 'CA']
    df = df[df['chain_id'].isin(chains)]
    coords = df[['x_coord', 'y_coord', 'z_coord']].to_numpy()

    pcaMHC = PCA()
    pcaMHC.fit(coords)
    pca1MHC = orient_vector(pcaMHC.components_[0], v_mhc) #MHC binding groove vector
    pca2MHC = orient_vector(pcaMHC.components_[1], v_mhc) #the waist of the MHC binding groove
    pca3MHC = np.cross(pca1MHC, pca2MHC) #the normal vector of MHC plane

    return pca1MHC, pca2MHC, pca3MHC


def get_TCR_vectors(pdb_path: str, chains: List[str] = ['D','E']) -> NDArray:

    # get elements we need
    v_alpha = get_end_ref(pdb_path,['D']) - CoM(pdb_path,['D'])
    v_beta = get_end_ref(pdb_path,['E']) - CoM(pdb_path,['E'])

    # Load pdb
    ppdb = PandasPdb().read_pdb(pdb_path)

    # Extract only CA atoms
    df = ppdb.df['ATOM']
    df = df[df['atom_name'] == 'CA']
    df = df[df['chain_id'].isin(chains)]
    coords_D = df[df['chain_id'] == 'D'][['x_coord', 'y_coord', 'z_coord']].to_numpy()
    coords_E = df[df['chain_id'] == 'E'][['x_coord', 'y_coord', 'z_coord']].to_numpy()

    pcaD = PCA()
    pcaD.fit(coords_D)
    pca1D = orient_vector(pcaD.components_[0], v_alpha)
    pca2D = orient_vector(pcaD.components_[1], v_alpha)
    pca3D = np.cross(pca1D, pca2D)

    pcaE = PCA()
    pcaE.fit(coords_E)
    pca1E = orient_vector(pcaE.components_[0], v_beta)
    pca2E = orient_vector(pcaE.components_[1], v_beta)
    pca3E = np.cross(pca1E, pca2E)

    coma_D = CoM(pdb_path,['D','E']) -  CoM(pdb_path,['D'])
    coma_D = coma_D / np.linalg.norm(coma_D)
    coma_E = CoM(pdb_path,['D','E']) -  CoM(pdb_path,['E'])
    coma_E = coma_E / np.linalg.norm(coma_E)

    return coma_D, pca1D, pca2D, pca3D, coma_E, pca1E, pca2E, pca3E


def calculate_incident_angle(vector1: NDArray, vector2: NDArray) -> float:

    dot_product = np.dot(vector1, vector2)
    magnitude_v1 = np.linalg.norm(vector1)
    magnitude_v2 = np.linalg.norm(vector2)
    cosine_theta = dot_product / (magnitude_v1 * magnitude_v2)
    cosine_theta = np.clip(cosine_theta, -1.0, 1.0)
    angle_radians = np.arccos(cosine_theta)

    return angle_radians


def min_ca_distances(pdb_path: str, chains: Tuple[str] = ('A', 'C', 'D', 'E')):
    """
    Compute minimum distance between CA atoms for each pair of chains.
    Parameters
    ----------
    pdb_path : str
        Path to the PDB file.
    chains : tuple
        Chain IDs to compare.
    Returns
    -------
    dict
        Dictionary with chain-pair keys and minimum CA–CA distances as values.
    """

    # Load pdb
    ppdb = PandasPdb().read_pdb(pdb_path)
    # Extract only CA atoms
    df = ppdb.df['ATOM']
    ca_atoms = df[df['atom_name'] == 'CA']
    # Store results
    min_dists = {}
    # Iterate over chain pairs
    for chain1, chain2 in itertools.combinations(chains, 2):
        coords1 = ca_atoms[ca_atoms['chain_id'] == chain1][['x_coord', 'y_coord', 'z_coord']].to_numpy()
        coords2 = ca_atoms[ca_atoms['chain_id'] == chain2][['x_coord', 'y_coord', 'z_coord']].to_numpy()
        if coords1.size == 0 or coords2.size == 0:
            min_dists[f"{chain1}{chain2}"] = np.nan
            continue
        # Compute pairwise distances efficiently using broadcasting
        diff = coords1[:, None, :] - coords2[None, :, :]
        dists = np.linalg.norm(diff, axis=2)
        min_dists[f"{chain1}{chain2}"] = np.min(dists)

    return min_dists


def get_centers_and_angles(pdb_path: str) -> Tuple[List[float], List[float]]:

    com_AC = CoM(pdb_path, ['A','C'])
    com_DE = CoM(pdb_path, ['D','E'])
    com_D = CoM(pdb_path, ['D'])
    com_E = CoM(pdb_path, ['E'])

    dists = [
        np.linalg.norm(com_AC - com_DE), np.linalg.norm(com_AC - com_D),
        np.linalg.norm(com_AC - com_E), np.linalg.norm(com_DE - com_D),
        np.linalg.norm(com_DE - com_E), np.linalg.norm(com_D - com_E)
    ]

    MHC_vecs = get_MHC_vectors(pdb_path)
    TCR_vecs = get_TCR_vectors(pdb_path)

    vecs = []
    for mhc_vec in MHC_vecs:
        vecs.append(mhc_vec)
    for tcr_vec in TCR_vecs:
        vecs.append(tcr_vec)

    angles = []
    for v1 in vecs:
        for v2 in vecs:
            angles.append(calculate_incident_angle(v1, v2))

    return angles, dists


def calc_geometric_features(folder_path: str, chains: Tuple[str] = ('A', 'C', 'D', 'E')):
    """
    For every pdb in a folder:
    - Compute min CA-CA distances between chain pairs
    - Compute geometric angles and distances (get_centers_and_angles)
    and append all results to scores.csv
    """

    # Chain pair labels
    pair_cols = [f"{c1}{c2}" for i, c1 in enumerate(chains) for c2 in chains[i+1:]]

    # Geometry feature column names
    angle_cols = [f"angle_{i+1}" for i in range(81)]  # 9 vectors × 9 vectors = 81 pairwise angles
    dist_cols = [f"geom_dist_{i+1}" for i in range(6)]

    df = pd.DataFrame(columns=["pdb_name"] + pair_cols + angle_cols + dist_cols)

    # Ensure columns exist
    # for col in pair_cols + angle_cols + dist_cols:
    #     if col not in scores_df.columns:
    #         scores_df[col] = None

    # Loop over pdbs
    for pdb_path in glob(os.path.join(folder_path, "*.pdb")):

        pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]

        # Compute CA-CA distances
        ca_distances = min_ca_distances(pdb_path, chains)

        # Compute geometry features
        try:
            angles, dists = get_centers_and_angles(pdb_path)
        except Exception as e:
            print(f"⚠️ Skipping {pdb_name} due to error: {e}")
            continue

        # # Find corresponding row
        # row_idx = df.index[df.iloc[:, 0].astype(str).str.contains(pdb_name, case=False, regex=False)]
        # if len(row_idx) == 0:
        #     print(f"⚠️ No matching row for {pdb_name}")
        #     continue

        # idx = row_idx[0]
        idx = len(df)
        df.loc[idx, "pdb_name"] = pdb_name

        # Update CA distance columns
        for pair, val in ca_distances.items():
            df.loc[idx, pair] = val

        # Update geometric distance columns
        for i, val in enumerate(dists):
            df.loc[idx, f"geom_dist_{i+1}"] = val

        # Update angle columns
        for i, val in enumerate(angles):
            df.loc[idx, f"angle_{i+1}"] = val

    # Save updated csv
    df_path = os.path.join(folder_path, "geo.csv")
    df.to_csv(df_path, index=False)
    print(f"✅ Updated geometry and distance features saved to {df_path}")

    return


def filter_pdbs(
    folder_path: str, cutoffs: Dict = CUTOFFS, directions: Dict = DIRECTIONS,
    applied: Dict = APPLIED
):

    df = pd.read_csv(os.path.join(folder_path, "geo.csv"))
    feature_cols = list(cutoffs.keys())

    # Default directions/applied flags
    if directions is None:
        directions = {col: True for col in feature_cols}  # default ≤
    if applied is None:
        applied = {col: True for col in feature_cols}

    # Apply cutoffs
    mask = np.ones(len(df), dtype=bool)
    for col in feature_cols:
        # Check if the feature column exists in the DataFrame before applying the cutoff
        if col in df.columns and applied[col]:
            if directions[col]:
                mask &= df[col] <= cutoffs[col]
            else:
                mask &= df[col] >= cutoffs[col]

    eliminated = df[~mask]

    for pdb_name in eliminated.pdb_name.values:
        os.remove(os.path.join(folder_path, pdb_name + ".pdb"))

    print(f"✅ {len(df) - len(eliminated)} / {len(df)} conformations passed geometric filtering")

    return
