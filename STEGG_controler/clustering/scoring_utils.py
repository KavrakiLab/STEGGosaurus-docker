import pandas as pd
import re
from io import StringIO # Used for reading string as a file

def get_ape_gen_scores(log_filepath):
    """
    Parses the Ape-GEN 2.0 log file to extract a table of peptide
    modeling results, specifically focusing on 'Peptide Index' and 'Affinity'.

    Args:
        log_filepath (str): The path to the Ape-GEN 2.0 output log file.

    Returns:
        pandas.DataFrame: A DataFrame containing 'Peptide_Index' (as a padded string)
                          and 'Affinity' (as a numeric value) for each modeled peptide.

    Raises:
        ValueError: If the table header is not found in the log file.
    """
    with open(log_filepath, 'r') as f:
        lines = f.readlines()

    table_start_index = None
    # Find the line where the table header starts
    for i, line in enumerate(lines):
        if line.strip().startswith('|') and 'Peptide index' in line:
            table_start_index = i
            break

    if table_start_index is None:
        raise ValueError(f"Table header 'Peptide index' not found in log file: {log_filepath}")

    # Collect all lines that are part of the table
    table_lines = []
    for line in lines[table_start_index:]:
        if line.strip().startswith('|'):
            table_lines.append(line.strip())
        else:
            break # Stop when a non-table line is encountered

    # Filter for lines that match the expected data format (e.g., | 0001 | Successfully Modeled | -123.45 |)
    relevant_table_lines = [
        line for line in table_lines
        if re.match(r"^\|\s{12}(\d{4})\s\|\sSuccessfully Modeled\s\|\s*(-?\d+\.\d+)\s*\|$", line)
    ]

    # Clean the table lines
    # Example transformation: "|    0001   | Successfully Modeled |   -123.45   |" becomes "0001,Successfully Modeled,-123.45"
    table_text_for_csv = '\n'.join(
        line.replace('|', ',').replace(" ", "").strip()[1:-1] # Remove leading/trailing '|' and spaces
        for line in relevant_table_lines
    )

    # Read the cleaned text into a DataFrame
    df = pd.read_csv(StringIO(table_text_for_csv), usecols=[0, 2], header=None)

    # Rename columns for clarity
    df.columns = ['Peptide_Index', 'Affinity']

    # Convert 'Peptide_Index' to string and pad with leading zeros to 4 digits (e.g., 1 -> '0001')
    df['Peptide_Index'] = df['Peptide_Index'].astype(str).str.zfill(4)
    # Convert 'Affinity' to a numeric type
    df['Affinity'] = pd.to_numeric(df['Affinity'])

    return df
