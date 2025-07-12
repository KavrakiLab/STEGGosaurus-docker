import string
import random
import gzip
import shutil
from pathlib import Path
import pandas as pd

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def generate_unique_string(existing_strings, postfix='', length=16):
    """
    Generates a unique random string of a specified length that is not present
    in a given list of existing strings.

    Args:
        existing_strings (list): A list of strings to check for uniqueness.
        length (int): The desired length of the unique string. Defaults to 16.

    Returns:
        str: A unique random string.
    """
    # Define the character set for string generation (alphanumeric)
    charset = string.ascii_letters + string.digits
    while True:
        # Generate a random string
        new_string = ''.join(random.choices(charset, k=length))
        # Check if the generated string is unique
        if new_string+postfix not in existing_strings:
            return new_string+postfix

def unzip_gz_files(source_directory, destination_directory, pair_identifier=''):
    """
    Unzips all .gz files from a source directory to a destination directory.
    Optionally prepends a pair identifier to the unzipped filenames.

    Args:
        source_directory (str): The path to the directory containing .gz files.
        destination_directory (str): The path to the directory where unzipped files
                                     will be saved.
        pair_identifier (str): An optional string to prepend to the unzipped filenames.
                               Defaults to ''.
    """
    source_path = Path(source_directory)
    destination_path = Path(destination_directory)
    # Create the destination directory if it does not exist
    destination_path.mkdir(parents=True, exist_ok=True)

    # Iterate over all .gz files in the source directory
    for gz_file_path in source_path.glob("*.gz"):
        # Construct the output filename: remove '.gz' and replace 'emref' with pair_identifier_complex
        output_filename = gz_file_path.stem.replace("emref", pair_identifier + "_complex")
        output_filepath = destination_path / output_filename

        # Open the gzipped file for reading in binary mode and the output file for writing in binary mode
        with gzip.open(gz_file_path, 'rb') as f_in, open(output_filepath, 'wb') as f_out:
            # Copy the content from the gzipped file to the output file
            shutil.copyfileobj(f_in, f_out)
        print(f"Unzipped {gz_file_path.name} -> {output_filepath.name}")
