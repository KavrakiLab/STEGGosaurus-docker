
from glob import glob
import json
import os
import shutil
import sys
from tqdm import tqdm

from clustering.get_best_diverse import get_ape_gen_diverse
from filtering.geometry import calc_geometric_features, filter_pdbs
from post_processing.align_pdbs import *
from post_processing.get_ROI import *
from utils import *


WORKING_DIR = os.path.abspath(os.getcwd())
APEGEN_WORKING_DIR = os.path.abspath("../Ape-Gen2.0-main")
TFOLD_WORKING_DIR = os.path.abspath("../tfold/")
TRECS_WORKING_DIR = os.path.abspath("../T-RECS/")
HADDOCK3_WORKING_DIR = os.path.abspath("../haddock3_TCR/")

PMHC_DB_DIR = "pMHC_ensemble_DB"
TCR_DB_DIR = "TCR_ensemble_DB"
TCRPMHC_DB_DIR = "STEGG_complex_DB"


if __name__ == "__main__":

    input_json_fp = sys.argv[1]
    output_dir = sys.argv[2]

    os.makedirs(output_dir, exist_ok=True)

    with open(input_json_fp, "r") as f:
        input_data = json.load(f)

    job_id = input_data["jobID"]
    peptide_seq = input_data["peptide"]
    mhc_seq = input_data["MHC"]
    tcr_alpha_seq = input_data["alpha"]
    tcr_beta_seq = input_data["beta"]
    cdr3a_seq = input_data["CDR3a"]
    cdr3b_seq = input_data["CDR3b"]

    pmhc_db_dir = os.path.join(output_dir, PMHC_DB_DIR, job_id)
    tcr_db_dir = os.path.join(output_dir, TCR_DB_DIR, job_id)
    tcrpmhc_db_dir = os.path.join(output_dir, TCRPMHC_DB_DIR, job_id)

    # --- Model pMHC with tFold + Ape-GEN 2.0 ---

    print("Modeling pMHC ensemble")

    apegen_output_file = os.path.join(APEGEN_WORKING_DIR, f"{job_id}_output.txt")
    apegen_output_dir = os.path.join(APEGEN_WORKING_DIR, "output/", job_id)

    os.chdir(APEGEN_WORKING_DIR)
    os.system(f"python3 New_APE-Gen.py {peptide_seq} {mhc_seq} --dir {os.path.relpath(apegen_output_dir)} --num_cores 50 --num_loops_for_optimization 200 --verbose > {apegen_output_file}")
    os.chdir(WORKING_DIR)

    os.makedirs(pmhc_db_dir)
    for pdb_fp in tqdm(glob(os.path.join(apegen_output_dir, "results/5_final_conformations/*.pdb"))):
        shutil.move(pdb_fp, pmhc_db_dir)

    shutil.rmtree(apegen_output_dir)

    best_diverse_pmhc = get_ape_gen_diverse(pmhc_db_dir, apegen_output_file, 10)
    print(best_diverse_pmhc)
    if best_diverse_pmhc is not None:
        for pdb_fp in tqdm(glob(os.path.join(pmhc_db_dir, "*.pdb"))):
            if os.path.basename(pdb_fp) not in list(best_diverse_pmhc["Peptide_Index"]):
                os.remove(pdb_fp)

    os.remove(apegen_output_file)

    metadata = dict(hash=job_id, peptide=peptide_seq, MHC=mhc_seq)
    with open(os.path.join(pmhc_db_dir, "meta.json"), "w") as f:
        json.dump(metadata, f)
    new_df = pd.DataFrame([metadata])
    new_df.to_csv(os.path.join(output_dir, PMHC_DB_DIR, "pMHC_ensemble_DB.csv"), mode="a", header=False, index=False)

    # --- Model TCR with tFold + T-RECS ---

    print("Modeling TCR ensemble")

    tcr_fasta_fp = os.path.join(TFOLD_WORKING_DIR, f"{job_id}_TCR.fasta")
    with open(tcr_fasta_fp, "w") as f:
        f.write(f">A\n{tcr_alpha_seq}\n")
        f.write(f">B\n{tcr_beta_seq}\n")

    trecs_input_dir = os.path.join(TRECS_WORKING_DIR, "input/")
    os.makedirs(trecs_input_dir, exist_ok=True)
    trecs_input_fp = os.path.join(trecs_input_dir, f"{job_id}.pdb")

    os.chdir(TFOLD_WORKING_DIR)
    os.system(f"python3 projects/tfold_tcr/predict.py --fasta {tcr_fasta_fp} --output {trecs_input_fp} --model_version TCR")
    os.chdir(WORKING_DIR)

    os.chdir(TRECS_WORKING_DIR)
    os.system(f"python3 T-RECS.py {trecs_input_fp} 25")
    os.chdir(WORKING_DIR)

    os.makedirs(tcr_db_dir)
    trecs_output_dir = os.path.join(TRECS_WORKING_DIR, "output/", job_id)
    for pdb_fp in glob(os.path.join(trecs_output_dir, "*.pdb")):
        shutil.move(pdb_fp, tcr_db_dir)

    os.remove(trecs_input_fp)
    shutil.rmtree(trecs_output_dir)

    metadata = dict(hash=job_id, alpha=tcr_alpha_seq, beta=tcr_beta_seq)
    with open(os.path.join(tcr_db_dir, "meta.json"), "w") as f:
        json.dump(metadata, f)
    new_df = pd.DataFrame([metadata])
    new_df.to_csv(os.path.join(output_dir, TCR_DB_DIR, "TCR_ensemble_DB.csv"), mode="a", header=False, index=False)

    # --- Model complexes with HADDOCK3 --- #

    print("Running ensemble docking")

    os.chdir(HADDOCK3_WORKING_DIR)
    os.system(f"python3 dock_TCR_pHLA.py {tcr_db_dir} {pmhc_db_dir} {tcr_alpha_seq} {tcr_beta_seq} {peptide_seq}")
    os.chdir(WORKING_DIR)

    haddock3_output_dir = os.path.join(HADDOCK3_WORKING_DIR, "output/", f"{job_id}_{job_id}")
    os.makedirs(tcrpmhc_db_dir)
    unzip_gz_files(os.path.join(haddock3_output_dir, "run/9_emref/"), tcrpmhc_db_dir)
    shutil.rmtree(haddock3_output_dir)

    for pdb_fp in glob(os.path.join(tcrpmhc_db_dir, "*.pdb")):
        decouple_chains(pdb_fp, tcrpmhc_db_dir, f"{os.path.basename(pdb_fp)[:-4]}_roi")
        os.remove(pdb_fp)

    for fp in glob(os.path.join(tcrpmhc_db_dir, "*anarci*")):
        os.remove(fp)

    align_all_pdbs(tcrpmhc_db_dir)

    metadata = dict(hash=job_id, alpha=tcr_alpha_seq, beta=tcr_beta_seq, peptide=peptide_seq, MHC=mhc_seq)
    with open(os.path.join(tcrpmhc_db_dir, "meta.json"), "w") as f:
        json.dump(metadata, f)
    new_df = pd.DataFrame([metadata])
    new_df.to_csv(os.path.join(output_dir, TCRPMHC_DB_DIR, "STEGG_complex_DB.csv"), mode="a", header=False, index=False)
