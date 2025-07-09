import pandas as pd
from biopandas.pdb import PandasPdb
import os
import shutil

from utils import *

def make_RCD_target(pdb,chain_id,loop_start,loop_end,filename='rcd_targets.txt',out_dir='output/'):
    """
    Generates a target file for RCD (Randomized Coordinate Decent) loop modeling.

    This function extracts the amino acid sequence of a specified loop region
    from a given PDB file and saves it in a format suitable for RCD.

    Args:
        pdb (str): Path to the input PDB file.
        chain_id (str): The chain identifier of the loop.
        loop_start (int): The starting residue number of the loop.
        loop_end (int): The ending residue number of the loop.
        filename (str, optional): The name of the output target file. Defaults to 'rcd_targets.txt'.
        out_dir (str, optional): The output directory for the target file. Defaults to 'output/'.

    Returns:
        str: The full path to the generated RCD target file.
    """
    ppdb = PandasPdb().read_pdb(pdb)
    atom_df = ppdb.df['ATOM']
    atom_df = atom_df[atom_df['chain_id']==chain_id].drop_duplicates(subset=['residue_number'])
    atom_df = atom_df[atom_df['residue_number'].isin([x for x in range(loop_start,loop_end+1)])]
    aa_list = atom_df['residue_name'].apply(lambda x: d[x]).values
    aa_seq = ''
    for aa in aa_list:
        aa_seq += aa
        
    with open(out_dir+filename,'w') as out_file:
        out_file.write(os.path.basename(pdb)+' '+str(loop_start)+' '+str(loop_end)+' '+chain_id+' '+aa_seq+'\n')
    return out_dir+filename
    
def run_RCD(n_samples=200,n_keep=1,dd=0.5,t=0.90,target_dir='output/',rcd_file_location='RCD_required_files',target_file_name='rcd_targets.txt'):
    """
    Executes the RCD (Randomized Coordinate Decent) program for loop modeling.

    This function changes the current working directory to the target directory,
    runs the RCD executable with specified parameters, and then reverts to the
    original working directory.

    Args:
        n_samples (int, optional): Number of samples to generate. Defaults to 200.
        n_keep (int, optional): Number of best conformations to keep. Defaults to 1.
        dd (float, optional): Distance deviation threshold. Defaults to 0.5.
        t (float, optional): Temperature parameter. Defaults to 0.90.
        target_dir (str, optional): Directory containing the RCD target file. Defaults to 'output/'.
        rcd_file_location (str, optional): Path to the RCD required files directory. Defaults to 'RCD_required_files'.
        target_file_name (str, optional): Name of the RCD target file. Defaults to 'rcd_targets.txt'.
    """
    home = os.getcwd()
    os.chdir(target_dir)
    rcd_file_location = home+'/'+rcd_file_location
    os.system(rcd_file_location+'/bin/rcd '+target_file_name+' -n '+str(n_samples)+' -r -t '+str(t)+' -d '+str(dd)+' --linear -x '+rcd_file_location+'/dunbrack.bin --energy_file '+rcd_file_location+'/korp6Dv1.bin --loco_best '+str(n_keep)+' --bench -o rcd_files')
    os.chdir(home)

def update_loop(loop_pdb,pdb,postfix='_rcd'):
    """
    Updates the main PDB file with the coordinates of a newly modeled loop.

    This function replaces the original loop region in the main PDB file with
    the atoms from the provided loop PDB file and renumbers the atoms.

    Args:
        loop_pdb (str): Path to the PDB file containing the modeled loop coordinates.
        pdb (str): Path to the original main PDB file.
        postfix (str, optional): Postfix to add to the output PDB filename. Defaults to '_rcd'.

    Returns:
        str: The full path to the updated PDB file.
    """
    ppdb = PandasPdb().read_pdb(pdb)
    atom_df = ppdb.df['ATOM']
    
    looppdb = PandasPdb().read_pdb(loop_pdb)
    loop_atom_df = looppdb.df['ATOM']
    
    chain = loop_atom_df['chain_id'][0]
    residues = loop_atom_df['residue_number'].unique()
    
    atom_df = atom_df[~((atom_df['chain_id'] == chain) & (atom_df['residue_number'].isin(residues)))]

    pos = atom_df.loc[(atom_df['chain_id'] == chain) & (atom_df['residue_number'] ==residues[0]-1)].index[-1]+1
  
    atom_df = pd.concat([atom_df.iloc[:pos],loop_atom_df,atom_df.iloc[pos:]])
    
    atom_df['atom_number'] = list(range(1,atom_df.shape[0]+1))
    atom_df['line_idx'] = list(range(1,atom_df.shape[0]+1))
    
    opdb = PandasPdb()
    opdb.df['ATOM'] = atom_df
    opdb.to_pdb(pdb[:-4]+postfix+'.pdb')
    return pdb[:-4]+postfix+'.pdb'

def pack_sidechains(pdb, chain_id, loop_start, loop_end, output_dir='output/'):
    """
    Repacks sidechains around a specified loop region using FASPR.

    This function generates a sequence file for FASPR, highlighting the loop
    region (uppercase) and other regions (lowercase), then runs FASPR to
    optimize sidechain conformations.

    Args:
        pdb (str): Path to the input PDB file.
        chain_id (str): The chain identifier of the loop.
        loop_start (int): The starting residue number of the loop.
        loop_end (int): The ending residue number of the loop.
        output_dir (str, optional): The output directory for temporary files and the packed PDB. Defaults to 'output/'.

    Returns:
        str: The full path to the PDB file with packed sidechains.
    """
    ppdb = PandasPdb().read_pdb(pdb)
    sequence = ppdb.amino3to1()
    sequence_list = list(sequence['residue_name'])
    print('seq:',sequence.head())
    atom_df = ppdb.df['ATOM']
    atom_df = atom_df.drop_duplicates(subset=['chain_id','residue_number'])
    print('df_shape:',atom_df.shape)
    aa_seq = ''
    for i,row in atom_df.iterrows():
        if row['chain_id'] == chain_id and row['residue_number'] >= loop_start and row['residue_number'] <= loop_end+1:
            aa_seq += d[row['residue_name']]
        else:
            aa_seq += d[row['residue_name']].lower()
        
    with open(output_dir+'faspr_seq.txt','w') as out_file:
        out_file.write(aa_seq+'\n')
    print('seq_len:',len(aa_seq))
        
    os.system('FASPR_required_files/FASPR -i '+pdb+' -o '+pdb[:-4]+'_packed.pdb -s '+output_dir+'faspr_seq.txt')
    return pdb[:-4]+'_packed.pdb'

def generate_loop_conf(working_pdb,fname,chain,output_dir,loop_start,loop_end,n_samples):
    """
    Generates a loop conformation using RCD and packs sidechains with FASPR.

    This is a comprehensive function that orchestrates the loop modeling process:
    1. Copies the input PDB to a working directory.
    2. Generates an RCD target file for the specified loop.
    3. Runs RCD to generate loop conformations.
    4. Updates the main PDB with the best RCD-modeled loop.
    5. Packs sidechains around the new loop using FASPR.

    Args:
        working_pdb (str): Path to the input PDB file.
        fname (str): Base filename for intermediate and output files.
        chain (str): The chain identifier of the loop.
        output_dir (str): The directory for all output and intermediate files.
        loop_start (int): The starting residue number of the loop.
        loop_end (int): The ending residue number of the loop.
        n_samples (int): Number of samples to generate for RCD.

    Returns:
        str: The full path to the final PDB file with the modeled loop and packed sidechains.
    """
    pdb_a1 = output_dir+fname+'.pdb'
    shutil.copyfile(working_pdb,pdb_a1)
    rcd_target = make_RCD_target(pdb_a1,chain,loop_start,loop_end,filename='rcd_targets.txt',out_dir=output_dir)
    run_RCD(n_samples=n_samples,n_keep=1,dd=0.5,t=0.90,target_dir=output_dir,rcd_file_location='RCD_required_files',target_file_name='rcd_targets.txt')
    rcd_loop = output_dir+'rcd_files/'+fname+'_closed.pdb'
    rcd_pdb = update_loop(rcd_loop,pdb_a1)
    packed_pdb = pack_sidechains(rcd_pdb,chain,loop_start,loop_end,output_dir)
    return packed_pdb
