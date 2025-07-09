# ****************************************************************
# *                       DockQ                                  *
# *   Docking scoring for biomolecular models                    *
# *   DockQ score legend:                                        *
# *    0.00 <= DockQ <  0.23 - Incorrect                         *
# *    0.23 <= DockQ <  0.49 - Acceptable quality                *
# *    0.49 <= DockQ <  0.80 - Medium quality                    *
# *            DockQ >= 0.80 - High quality                      *
# *   Ref: Mirabello and Wallner, 'DockQ v2: Improved automatic  *
# *   quality measure for protein multimers, nucleic acids       *
# *   and small molecules'                                       *
# *                                                              *
# *   For comments, please email: bjorn.wallner@.liu.se          *
# ****************************************************************


import os
import matplotlib.pyplot as plt
import numpy as np
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces

def merge_chains(model, chains_to_merge):
    for chain in chains_to_merge[1:]:
        for i, res in enumerate(model[chain]):
            res.id = (chain, res.id[1], res.id[2])
            model[chains_to_merge[0]].add(res)
        model.detach_child(chain)
    model[chains_to_merge[0]].id = "".join(chains_to_merge)
    return model

def get_dock_q_score(native_pdb,model_pdb):
    # native = load_PDB("../../dockq_testing/1ao7_trimmed.pdb")
    # model = load_PDB("../../dockq_testing/1ao7_initial_models/final_complex_complex_1.pdb")
    native = load_PDB(native_pdb)
    model = load_PDB(model_pdb)

    model = merge_chains(model, ["A", "C"])
    model = merge_chains(model, ["D", "E"])
    native = merge_chains(native, ["A", "C"])
    native = merge_chains(native, ["D", "E"])

    chain_map = {"AC":"AC", "DE":"DE"}

    # returns a dictionary containing the results and the total DockQ score
    res = run_on_all_native_interfaces(model, native,chain_map)
    complex_dockq_score = res[0]['ACDE']
    return complex_dockq_score

def process_all_models(native_pdb, model_dir):
    dockq_scores = []
    scores_with_names = []

    for fname in sorted(os.listdir(model_dir)):
        if fname.endswith(".pdb"):
            model_path = os.path.join(model_dir, fname)
            try:
                score = get_dock_q_score(native_pdb, model_path)
                dockq_scores.append(score['DockQ'])
                scores_with_names.append((fname, score))
                print(f"{fname}: DockQ = {score['DockQ']:.3f}")
            except Exception as e:
                print(f"Error processing {fname}: {e}")

    if not dockq_scores:
        print("No scores computed.")
        return

    # Compute and print percentiles
    percentiles = np.percentile(dockq_scores, [0, 25, 50, 75, 100])
    print("\nDockQ Score Percentiles:")
    for label, val in zip(["Min", "25%", "Median", "75%", "Max"], percentiles):
        print(f"{label}: {val:.3f}")

    # Find highest scoring complex
    best_model, best_score = max(scores_with_names, key=lambda x: x[1]['DockQ'])
    print(f"\nHighest scoring model: {best_model} with the following metrics:")
    for key, value in best_score.items():
        print(f"Key: {key}, Value: {value}")

# Example usage:
native_pdb = "../../dockq_testing/1ao7_trimmed.pdb"
model_dir = "../../dockq_testing/1ao7_initial_models"
process_all_models(native_pdb, model_dir)
