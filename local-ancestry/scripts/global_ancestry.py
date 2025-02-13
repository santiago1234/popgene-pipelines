import sys
from pathlib import Path
from collections import defaultdict
import re

import pandas as pd


def load_bed_data(bed_hap0, bed_hap1):
    # Convert paths to Path objects
    bed_hap0 = Path(bed_hap0)
    bed_hap1 = Path(bed_hap1)

    # Extract sample names by removing '_0.bed' and '_1.bed'
    sample_0 = re.sub(r"_0\.bed$", "", bed_hap0.name)
    sample_1 = re.sub(r"_1\.bed$", "", bed_hap1.name)

    # Ensure they refer to the same sample
    assert (
        sample_0 == sample_1
    ), f"Sample names do not match: '{sample_0}' vs '{sample_1}'"

    # Load BED files
    df_hap0 = pd.read_csv(bed_hap0, sep="\t", na_values=[".", "NA"])
    df_hap1 = pd.read_csv(bed_hap1, sep="\t", na_values=[".", "NA"])

    # Concatenate and return
    bed_data = pd.concat([df_hap0, df_hap1], ignore_index=True)
    bed_data["sample"] = sample_0
    return bed_data


def compute_ancestry_proportions(bed_data):
    """
    Compute ancestry proportions based on weighted genomic segments.

    Parameters:
    bed_data (pd.DataFrame): A DataFrame containing columns ['sample', 'ancestry', 'spos', 'epos'].

    Returns:
    pd.DataFrame: A DataFrame where each row is a sample and each column is an ancestry proportion.
    """
    # Compute weight (genomic segment length)
    bed_data["weight"] = bed_data["epos"] - bed_data["spos"]

    # Group by sample and ancestry, summing weights
    ancestry_weights = bed_data.groupby(["sample", "ancestry"])["weight"].sum()

    # Compute ancestry proportions
    ancestry_proportions = ancestry_weights / ancestry_weights.groupby(level=0).sum()

    # Reshape: pivot to have one row per sample and one column per ancestry
    ancestry_proportions_df = ancestry_proportions.unstack().reset_index()

    return ancestry_proportions_df


def list_bed_files(dir_with_bed_files):
    """
    Lists and pairs BED files from a given directory.

    Parameters:
    dir_with_bed_files (str or Path): Directory containing the BED files.

    Returns:
    dict: A dictionary where keys are sample names and values are lists containing two BED file paths.
    """
    # Convert to Path object if needed
    dir_with_bed_files = Path(dir_with_bed_files)

    # Get all .bed files in the directory
    bed_files = list(dir_with_bed_files.glob("*.bed"))

    # Dictionary to store sample -> [bed_hap0, bed_hap1]
    bed_dict = defaultdict(list)

    # Populate dictionary
    for bed_file in bed_files:
        # Extract sample name by removing "_0.bed" and "_1.bed"
        sample_name = bed_file.name.replace("_0.bed", "").replace("_1.bed", "")

        # Append the file to the corresponding sample
        bed_dict[sample_name].append(bed_file)

    # Ensure each sample has exactly two files (one for hap0 and one for hap1), and sort them
    bed_dict = {
        sample: sorted(files) for sample, files in bed_dict.items() if len(files) == 2
    }

    return bed_dict


def main():
    # Ensure the correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python script.py <bed_directory> <output_file>")
        sys.exit(1)

    # Parse command-line arguments
    bed_dir, out_file = sys.argv[1:3]

    # Get all paired BED files
    bed_files = list_bed_files(bed_dir)

    # Load and process each sample
    bed_data = [load_bed_data(*files) for files in bed_files.values()]
    ancestry_proportions = [compute_ancestry_proportions(data) for data in bed_data]

    # Concatenate results and handle NaNs
    ancestry_proportions = pd.concat(ancestry_proportions, ignore_index=True)
    ancestry_proportions.fillna(0, inplace=True)

    # Save results to CSV
    ancestry_proportions.to_csv(out_file, index=False)
    print(f"Ancestry proportions saved to {out_file}")


# Ensure script runs only when executed directly
if __name__ == "__main__":
    main()
