"""
Script to create a filemap for running pong from ADMIXTURE Q files.

This script generates a filemap by scanning a directory for Q files 
matching a specific prefix. The resulting filemap includes run IDs, 
K values, and file paths, formatted for pong.

Usage:
    python make_pong_filemap.py <path_to_qfiles> <nameprefix> <outfile>

Arguments:
    path_to_qfiles (str): Path to the directory containing Q files.
    nameprefix (str): Prefix for the names of Q files to include.
    outfile (str): Path to save the output filemap.

Requirements:
    - Python 3.6 or higher
    - pandas library

Example:
    python make-pong-filemap.py ./qfiles admix output_filemap.txt
"""

import sys
from pathlib import Path

import pandas as pd


def make_pong_filemap(path_to_qfiles, nameprefix, outfile):
    """
    Create a filemap for running pong.

    Args:
        path_to_qfiles (str): Path to the directory containing Q files.
        nameprefix (str): Prefix for the names of Q files to include.
        outfile (str): Path to save the output filemap, formatted for pong.

    Returns:
        pandas.DataFrame: A DataFrame representing the filemap, saved to the outfile.
    """
    path_to_qfiles = Path(path_to_qfiles)

    # Validate input path
    if not path_to_qfiles.is_dir():
        raise ValueError(f"The path '{path_to_qfiles}' is not a valid directory.")

    # Collect Q files matching the prefix
    qfiles = [x for x in path_to_qfiles.glob("*Q") if x.name.startswith(nameprefix)]
    if not qfiles:
        raise ValueError(
            f"No Q files found with prefix '{nameprefix}' in '{path_to_qfiles}'."
        )

    # Extract K values and construct run IDs
    try:
        k_values = [int(x.name.split(".")[1]) for x in qfiles]
    except (IndexError, ValueError) as e:
        raise ValueError(
            "Failed to parse K values from Q file names. Ensure file names are in the format 'prefix.K.Q'."
        ) from e

    run_ids = [f"{qf.stem.replace('.', 'x')}-{k}" for qf, k in zip(qfiles, k_values)]

    # Create the filemap DataFrame
    filemap = pd.DataFrame(
        {
            "runID": run_ids,
            "Kvalue": k_values,
            "filepath": [str(qf) for qf in qfiles],  # Ensure paths are strings
        }
    ).sort_values("Kvalue")

    # Save the filemap to the output file
    try:
        filemap.to_csv(outfile, sep="\t", index=False, header=False)
    except Exception as e:
        raise IOError(f"Failed to write filemap to '{outfile}'.") from e

    return filemap


def main():
    """
    Main function to parse command-line arguments and create the filemap.
    Usage:
        python make_pong_filemap.py <path_to_qfiles> <nameprefix> <outfile>
    """
    if len(sys.argv) != 4:
        print(
            "Usage: python make_pong_filemap.py <path_to_qfiles> <nameprefix> <outfile>"
        )
        sys.exit(1)

    path_to_qfiles = sys.argv[1]
    nameprefix = sys.argv[2]
    outfile = sys.argv[3]

    try:
        make_pong_filemap(path_to_qfiles, nameprefix, outfile)
        print(f"Filemap successfully created and saved to '{outfile}'")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
