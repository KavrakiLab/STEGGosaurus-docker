
from glob import glob
from multiprocessing import Pool
import os
import sys
from tqdm import tqdm

import pandas as pd

import pyrosetta
pyrosetta.init("-out:levels basic:error core:error")

from pyrosetta import get_fa_scorefxn, pose_from_pdb


def score_ensemble(pdb_dir: str):

    pdb_fps = glob(os.path.join(pdb_dir, "*_complex_*_roi.pdb"))
    model_ids = [os.path.basename(pdb_fp) for pdb_fp in pdb_fps]
    models = list(zip(model_ids, pdb_fps))

    fxn = get_fa_scorefxn()
    score_types = fxn.get_nonzero_weighted_scoretypes()
    terms = [score_type.name for score_type in score_types]

    vectors = dict()
    for model_id, pdb_fp in tqdm(models):
        pose = pose_from_pdb(pdb_fp)
        fxn(pose)
        vectors[model_id] = [
            pose.energies().total_energies()[score_type]
            for score_type in score_types
        ]

    data = dict()
    for i, term in enumerate(terms):
        data[term] = {model_id: vectors[model_id][i] for model_id in model_ids}

    df = pd.DataFrame.from_dict(data)
    df.to_csv(os.path.join(pdb_dir, "energy_terms.csv"))

    return


if __name__ == "__main__":

    ensemble_dirs = glob(os.path.join(sys.argv[1], "*"))

    with Pool(processes=8) as p:
        p.map(score_ensemble, ensemble_dirs)
