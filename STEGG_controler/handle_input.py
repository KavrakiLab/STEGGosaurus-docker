import os
from anarci import anarci
import argparse
import json

def main():
    parser = argparse.ArgumentParser(description="Write TCR-related sequences to a JSON file.")
    parser.add_argument('--jobID', required=True, help="unique ID for job")
    parser.add_argument('--TCRalpha', required=True, help="TCR alpha chain sequence")
    parser.add_argument('--TCRbeta', required=True, help="TCR beta chain sequence")
    parser.add_argument('--peptide', required=True, help="Peptide sequence")
    parser.add_argument('--MHC', required=True, help="MHC sequence")
    
    args = parser.parse_args()

    CDR3a = get_CDR3_seq(args.TCRalpha)
    CDR3b = get_CDR3_seq(args.TCRbeta)

    data = {
        'alpha': args.TCRalpha,
        'beta': args.TCRbeta,
        'peptide': args.peptide,
        'MHC': args.MHC,
        'CDR3a': CDR3a,
        'CDR3b': CDR3b,
        'jobID': args.jobID
    }

    with open('input/'+args.jobID+'.json', 'w') as f:
        json.dump(data, f, indent=4)

    print("JSON file 'input/'+args.jobID+'.json' has been written successfully.")

def get_CDR3_seq(TCR_chain_seq):
    """
    """
    TCR_seq = [("input_seq",TCR_chain_seq)]
    anarci_out = anarci(TCR_seq,scheme='imgt',output=False)
    # Extract the list of tuples
    try:
        sequence_data = anarci_out[0][0][0][0]
    except:
        print("Not a valid TCR sequence!")
        raise
    
    result = ''.join(char for ((index, _), char) in sequence_data if 105 <= index <= 117)

    return result.replace('-','')

# TODO: verify peptide input using MHCflurry (already verified compatabililty with docker setup)

if __name__ == "__main__":
    main()