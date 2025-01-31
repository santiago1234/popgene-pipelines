"""
authors: S.G.M.M & C.B.J
Mask VCF by local ancestry. Masked genotypes are encoded as missing values.

Usage:
    python mask-vcf-by-ancestry.py <ANC> <vcf> <lai> <output>

Args:
    ANC: Ancestry code, see the first line of input <lai> file.
    vcf: VCF file containing genotype data.
    lai: Local ancestry file (gnomix *.msp output).
    output: Output name for the masked VCF file.

NOTES:
    - Missing genotypes are encoded as -1.
    - The script assumes the VCF and MSP files contain data for a single chromosome.
"""

import sys
import os
import logging
import argparse
from typing import Tuple, Set

import pandas as pd
from cyvcf2 import VCF, Writer
import allel

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


def validate_input_files(vcf_file: str, msp_file: str) -> None:
    """
    Validate that the input files exist and are accessible.

    Args:
        vcf_file (str): Path to the VCF file.
        msp_file (str): Path to the MSP file.

    Raises:
        FileNotFoundError: If either file does not exist.
    """
    if not os.path.exists(vcf_file):
        raise FileNotFoundError(f"VCF file {vcf_file} does not exist.")
    if not os.path.exists(msp_file):
        raise FileNotFoundError(f"MSP file {msp_file} does not exist.")


def get_vcf_range(vcf_file: str) -> Tuple[int, int, str]:
    """
    Get the range of variants and the chromosome present in the VCF file.

    Args:
        vcf_file (str): Path to the VCF file.

    Returns:
        Tuple[int, int, str]: Minimum position, maximum position, and chromosome.

    Raises:
        ValueError: If the VCF file contains more than one chromosome.
    """
    callset = allel.read_vcf(vcf_file, fields=["variants/CHROM", "variants/POS"])
    min_position = callset["variants/POS"].min()
    max_position = callset["variants/POS"].max()

    chromosomes = set(callset["variants/CHROM"])
    if len(chromosomes) != 1:
        raise ValueError("The VCF file must contain data for exactly one chromosome.")

    return min_position, max_position, chromosomes.pop()


def load_and_restrict_msp(msp_file: str, vcf_file: str) -> pd.DataFrame:
    """
    Load the MSP file and restrict it to the range of variants present in the VCF file.

    Args:
        msp_file (str): Path to the MSP file.
        vcf_file (str): Path to the VCF file.

    Returns:
        pd.DataFrame: Filtered local ancestry data.

    Raises:
        ValueError: If the MSP file contains more than one chromosome or if the chromosomes in the VCF and MSP files do not match.
    """
    min_position, max_position, vcf_chromosome = get_vcf_range(vcf_file)
    ancestry_data = pd.read_csv(msp_file, sep="\t", skiprows=[0])

    chromosomes_in_ancestry = set(ancestry_data["#chm"].unique())
    if len(chromosomes_in_ancestry) != 1:
        raise ValueError("The MSP file must contain data for exactly one chromosome.")

    msp_chromosome = chromosomes_in_ancestry.pop()
    if str(vcf_chromosome) != str(msp_chromosome):
        raise ValueError(
            f"Chromosome mismatch: VCF has {vcf_chromosome}, MSP has {msp_chromosome}."
        )

    # Filter MSP data to the range of the VCF file
    ancestry_range = ancestry_data[
        (ancestry_data["spos"] <= max_position)
        & (ancestry_data["epos"] >= min_position)
    ]
    return ancestry_range.drop_duplicates()


def retrieve_lai_at(pos: int, chrom: str, lai_data: pd.DataFrame) -> pd.DataFrame:
    """
    Retrieve the local ancestry information for a specific position and chromosome.

    Args:
        pos (int): Genomic position.
        chrom (str): Chromosome.
        lai_data (pd.DataFrame): Local ancestry data.

    Returns:
        pd.DataFrame: Local ancestry information for the specified position.

    Raises:
        ValueError: If the position is not found or is in multiple ranges.
    """
    lai_pos = lai_data[(lai_data["spos"] <= pos) & (lai_data["epos"] >= pos)]
    lai_pos = lai_pos[lai_pos["#chm"].map(str) == str(chrom)]

    if lai_pos.shape[0] < 1:
        raise ValueError(f"Position {pos} not found in local ancestry data.")
    if lai_pos.shape[0] > 1:
        raise ValueError(f"Position {pos} is in multiple ranges in the local ancestry data.")

    return lai_pos


def sample_ancestry_at(lai_pos: pd.DataFrame, sample_index: int, sample_name_map: dict) -> Tuple[str, str]:
    """
    Get the ancestry of a specific sample at a given position.

    Args:
        lai_pos (pd.DataFrame): Local ancestry information for the position.
        sample_index (int): Index of the sample in the VCF file.
        sample_name_map (dict): Mapping of sample indices to sample names.

    Returns:
        Tuple[str, str]: Ancestry codes for the two haplotypes of the sample.
    """
    sample_name = sample_name_map[sample_index]
    sample_hap_0 = f"{sample_name}.0"
    sample_hap_1 = f"{sample_name}.1"

    anc_hap_0 = lai_pos[sample_hap_0].to_list()[0]
    anc_hap_1 = lai_pos[sample_hap_1].to_list()[0]

    return str(anc_hap_0), str(anc_hap_1)


def main():
    """Main function to mask VCF by local ancestry."""
    parser = argparse.ArgumentParser(description="Mask VCF by local ancestry.")
    parser.add_argument("ANC", type=str, help="Ancestry code to mask against.")
    parser.add_argument("vcf", type=str, help="Input VCF file.")
    parser.add_argument("lai", type=str, help="Local ancestry file (gnomix *.msp output).")
    parser.add_argument("output", type=str, help="Output VCF file.")
    args = parser.parse_args()

    # Validate input files
    validate_input_files(args.vcf, args.lai)

    # Load and restrict local ancestry data
    logger.info("Loading and restricting local ancestry data...")
    lai_data = load_and_restrict_msp(args.lai, args.vcf)

    # Open VCF file and initialize writer
    vcf_reader = VCF(args.vcf)
    vcf_writer = Writer(args.output, vcf_reader)
    sample_name_map = {x: y for (x, y) in enumerate(vcf_reader.samples)}

    # Process each variant
    logger.info("Masking genotypes...")
    for i, variant in enumerate(vcf_reader):
        chrom, pos = variant.CHROM, variant.POS

        if i % 10000 == 0:
            logger.info(f"Processing variant at {chrom}:{pos}...")

        try:
            lai_pos = retrieve_lai_at(pos, chrom, lai_data)
            for sample_index, _ in enumerate(variant.genotypes):
                anc_h0, anc_h1 = sample_ancestry_at(lai_pos, sample_index, sample_name_map)
                if anc_h0 != args.ANC:
                    variant.genotypes[sample_index][0] = -1
                if anc_h1 != args.ANC:
                    variant.genotypes[sample_index][1] = -1
            vcf_writer.write_record(variant)
        except ValueError as e:
            logger.warning(f"Skipping variant at {chrom}:{pos}: {e}")

    # Close files
    vcf_writer.close()
    vcf_reader.close()
    logger.info(f"Masked VCF saved to {args.output}.")


if __name__ == "__main__":
    main()
