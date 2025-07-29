import sys
import os
import pandas as pd

var = int(sys.argv[1])-1

df = pd.read_csv('stegg_new_benchmark_dataset.csv')

job_id = df.iloc[var]['job_id']
alpha = df.iloc[var]['TCR_A']
beta = df.iloc[var]['TCR_B']
pep = df.iloc[var]['peptide']
MHC = df.iloc[var]['MHC_A']

os.system('python handle_input.py --jobID '+job_id+' --TCRalpha '+alpha+' --TCRbeta '+beta+' --peptide '+pep+' --MHC '+MHC)
os.system('python model_complex.py input/'+job_id+'.json')