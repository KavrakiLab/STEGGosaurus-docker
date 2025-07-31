import json
import os
import shutil
import sys
import pandas as pd
from biopandas.pdb import PandasPdb

from utils import *
sys.path.insert(1,'/home/STEGG_controler/clustering/')
from get_best_diverse import *
sys.path.insert(1,'/home/STEGG_controler/post_processing/')
from get_ROI import *
from align_pdbs import *

if __name__ == "__main__":
    # Command-line arguments:
    # sys.argv[1]: Path to the input JSON file
    input_json = sys.argv[1]


    # Read JSON data from file
    with open(input_json, 'r') as f:
        data = json.load(f)

    # Assign sequences to variables
    peptide = data["peptide"]
    MHC = data["MHC"]
    alpha = data["alpha"]
    beta = data["beta"]
    CDR3a = data["CDR3a"]
    CDR3b = data["CDR3b"]

    TCR_pMHC_pair_id = data["jobID"]

    print('Modeling pMHC ensemble')
    # --- Model pMHC with tfold + Ape-GEN 2.0 ---
    os.chdir('../Ape-Gen2.0-main')
    os.system('python3 New_APE-Gen.py '+peptide+' '+MHC+' --dir '+TCR_pMHC_pair_id+' --num_cores 50 --num_loops_for_optimization 200 --verbose > ape-gen_output.txt')
    os.chdir('../STEGG_controler')

    shutil.move('../Ape-Gen2.0-main/'+TCR_pMHC_pair_id+'/results/5_final_conformations','pMHC_ensemble_DB/'+TCR_pMHC_pair_id)

    # Get the best diverse pMHC conformations based on Ape-GEN scores
    best_diverse_pmhc = get_ape_gen_diverse('pMHC_ensemble_DB/'+TCR_pMHC_pair_id,'../Ape-Gen2.0-main/ape-gen_output.txt',10)

    if best_diverse_pmhc is not None:
        for file in os.listdir('pMHC_ensemble_DB/'+TCR_pMHC_pair_id):
            if file not in list(best_diverse_pmhc['Peptide_Index']):
                os.remove('pMHC_ensemble_DB/'+TCR_pMHC_pair_id+'/'+file)

    # TODO: repeat modeling if too few conformations

    # Save metadata for the new pMHC ensemble
    meta = {'hash':TCR_pMHC_pair_id,'peptide':peptide,'MHC':MHC}
    with open('pMHC_ensemble_DB/'+TCR_pMHC_pair_id+'/meta.json','w') as f:
        json.dump(meta, f)
    new_df = pd.DataFrame([meta])
    new_df.to_csv('pMHC_ensemble_DB/pMHC_ensemble_DB.csv', mode='a', header=False, index=False)

    print('modeling TCR ensemble')
    os.makedirs('../T-RECS/output/'+TCR_pMHC_pair_id, exist_ok=True)
    # --- Model TCR ensemble with tfold + T-RECS ---
    os.chdir('../tfold')
    with open('TCR.fasta','w') as f:
        f.write('>A\n'+alpha+'\n')
        f.write('>B\n'+beta+'\n')
    # Run tfold to predict TCR structure
    os.system('python3 projects/tfold_tcr/predict.py --fasta TCR.fasta --output ../T-RECS/output/'+TCR_pMHC_pair_id+'/'+TCR_pMHC_pair_id+'TCR.pdb --model_version TCR')
    os.chdir('../STEGG_controler')

    # # Move tfold output to T-RECS directory and run T-RECS for conformational sampling
    # shutil.move('../tfold/output/TCR.pdb','../T-RECS/TCR.pdb')
    os.chdir('../T-RECS')
    os.system('python3 T-RECS.py output/'+TCR_pMHC_pair_id+'/'+TCR_pMHC_pair_id+'TCR.pdb 25')# start with 40 TCR conformations if filtering
    os.chdir('../STEGG_controler')

    shutil.move('../T-RECS/output/output/'+TCR_pMHC_pair_id+'/'+TCR_pMHC_pair_id+'TCR','TCR_ensemble_DB/'+TCR_pMHC_pair_id)
    # os.remove('../T-RECS/TCR.pdb')

    # # clustering and filtering of TCR conformations,
    # best_diverse_tcr = get_t_recs_diverse('TCR_ensemble_DB/'+TCR_pMHC_pair_id,25)
    #
    # for file in os.listdir('TCR_ensemble_DB/'+TCR_pMHC_pair_id):
    #     if file not in list(best_diverse_tcr['tcr_pdb_name']):
    #         os.remove('TCR_ensemble_DB/'+TCR_pMHC_pair_id+'/'+file)

    # TODO: re-modeling if needed.

    # Save metadata for the new TCR ensemble
    meta = {'hash':TCR_pMHC_pair_id,'alpha':alpha,'beta':beta}
    with open('TCR_ensemble_DB/'+TCR_pMHC_pair_id+'/meta.json','w') as f:
        json.dump(meta, f)
    new_df = pd.DataFrame([meta])
    new_df.to_csv('TCR_ensemble_DB/TCR_ensemble_DB.csv', mode='a', header=False, index=False)

    print('running ensemble docking')
    # --- Dock TCRs and pHLAs using HADDOCK3 ---
    os.chdir('../haddock3_TCR')
    cmd = 'python3 dock_TCR_pHLA.py ../STEGG_controler/TCR_ensemble_DB/'+TCR_pMHC_pair_id+' ../STEGG_controler/pMHC_ensemble_DB/'+TCR_pMHC_pair_id+' '+CDR3a+' '+CDR3b+' '+peptide
    print(cmd)
    os.system(cmd)
    os.chdir('../STEGG_controler')

    # --- Post-processing of docked complexes ---
    unzip_gz_files('../haddock3_TCR/output/'+TCR_pMHC_pair_id+'_'+TCR_pMHC_pair_id+'/run/9_emref', 'STEGG_complex_DB/'+TCR_pMHC_pair_id)
    for new_file in os.listdir('STEGG_complex_DB/'+TCR_pMHC_pair_id):
        decouple_chains('STEGG_complex_DB/'+TCR_pMHC_pair_id+'/'+new_file,'STEGG_complex_DB/'+TCR_pMHC_pair_id,new_file[:-4]+'_roi')
        os.remove('STEGG_complex_DB/'+TCR_pMHC_pair_id+'/'+new_file)
    for trash_file in os.listdir('STEGG_complex_DB'):
        if 'anarci' in trash_file:
            os.remove('STEGG_complex_DB/'+trash_file)

    #TODO energy min... ?

    # Align all structures based on a reference chain (MHC chain)
    align_all_pdbs('STEGG_complex_DB/'+TCR_pMHC_pair_id)

    # Save metadata for the new STEGG complex ensemble
    meta = {'hash':TCR_pMHC_pair_id,'alpha':alpha,'beta':beta,'peptide':peptide,'MHC':MHC}
    with open('STEGG_complex_DB/'+TCR_pMHC_pair_id+'/meta.json','w') as f:
        json.dump(meta, f)
    new_df = pd.DataFrame([meta])
    new_df.to_csv('STEGG_complex_DB/STEGG_complex_DB.csv', mode='a', header=False, index=False)

    # print("cleaning files")
    # #TODO: clean excess files

    print("\nComplex modeling and docking complete!")

    print('cleaning files')
    files_to_clean = [
        '../haddock3_TCR/'+TCR_pMHC_pair_id+'_'+TCR_pMHC_pair_id+'.cfg',
    ]
    directories_to_clean = [
        '../Ape-Gen2.0-main/'+TCR_pMHC_pair_id,
        '../T-RECS/output/'+TCR_pMHC_pair_id,
        '../T-RECS/output/output/'+TCR_pMHC_pair_id,
        '../haddock3_TCR/output/'+TCR_pMHC_pair_id+'_'+TCR_pMHC_pair_id,
        '../tfold/output/'+TCR_pMHC_pair_id+'/MODELLER_outputMHC.pdb',
    ]

    for file in files_to_clean:
        os.remove(file)
    
    for directory in directories_to_clean:
        shutil.rmtree(directory)
