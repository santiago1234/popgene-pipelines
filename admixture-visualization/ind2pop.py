"""
Generate files for population mapping from sample metadata.

This script processes a population information file and a sample order file 
to create two output files:
1. 'ind2pop.txt': A mapping of individuals to populations (or other covariates).
2. 'pop_order.txt': A list of unique populations (or covariates) in the order they appear.

Usage:
    python script.py <popinfo> <sample_order> <covariate> <sample_id_col>

Arguments:
    popinfo (str): Path to the population info file (CSV format).
    sample_order (str): Path to the sample order file (space-separated format).
    covariate (str): Column name in popinfo to use as population/covariate.
    sample_id_col (str): Column name in popinfo containing sample IDs.
"""
import sys

import pandas as pd


def main():
    if len(sys.argv) != 5:
        print(
            "Usage: python script.py <popinfo> <sample_order> <covariate> <sample_id_col>"
        )
        sys.exit(1)

    # Command-line arguments
    popinfo = sys.argv[1]
    sample_order = sys.argv[2]
    covariate = sys.argv[3]
    sample_id_col = sys.argv[4]

    # Fixed output file names
    outfile = "ind2pop.txt"
    pop_order_file = "pop_order.txt"  # File to save unique covariates

    # Load the population info and sample order
    popinfo = pd.read_csv(popinfo)
    sample_order = pd.read_csv(sample_order, header=None, sep=" ")
    # The 2nd column contains the sample IDs
    sample_order = sample_order.iloc[:, 1].values.tolist()

    # Ensure each sample in sample_order exists in popinfo
    assert all(
        sample in popinfo[sample_id_col].values for sample in sample_order
    ), "Some samples in sample_order are missing in popinfo"

    # Filter and reorder popinfo to match the sample_order
    popinfo_filtered = popinfo[popinfo[sample_id_col].isin(sample_order)]
    popinfo_filtered = (
        popinfo_filtered.set_index(sample_id_col).loc[sample_order].reset_index()
    )

    # Extract the desired covariate (e.g., population label) for each sample
    ind2pop_data = popinfo_filtered[covariate]

    # Replace spaces by underscores
    ind2pop_data = ind2pop_data.str.replace(" ", "_")

    # Save the ind2pop file
    with open(outfile, "w") as f:
        for label in ind2pop_data:
            f.write(f"{label}\n")

    print(f"ind2pop file successfully created: {outfile}")

    # Save unique covariates to the pop order file
    unique_covariates = ind2pop_data.unique()

    with open(pop_order_file, "w") as f:
        for covariate in unique_covariates:
            f.write(f"{covariate}\n")

    print(f"Population order file successfully created: {pop_order_file}")


if __name__ == "__main__":
    main()
