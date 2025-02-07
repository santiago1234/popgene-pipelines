"""
Test ancestry masking
"""

import pandas as pd
import allel
import numpy as np


def sampler(callset):
    """
    Sample a random POS index, its actual POS value, and a random SAMPLE index from the VCF dataset.

    Returns:
        pos_idx (int): Index of the randomly chosen variant in callset['variants/POS'].
        pos_value (int): Actual POS value at that index.
        sample_idx (int): Index of the randomly chosen sample.
    """
    positions = callset['variants/POS']
    num_samples = len(callset['samples'])

    pos_idx = np.random.randint(0, len(positions))
    pos_value = positions[pos_idx]
    sample_idx = np.random.randint(0, num_samples)

    return pos_idx, pos_value, sample_idx


def get_alleles(rpos_idx, rsample_idx):
    """
    Get the genotype data from callset['calldata/GT'] for a given variant position index and sample index.

    Returns:
        tuple: (allele_1, allele_2) or (-1, -1) for missing data.
    """
    return tuple(callset['calldata/GT'][rpos_idx, rsample_idx, :])


def locate_ancestry_region(pos, msp):
    """
    Find the ancestry region for a given genomic position.

    Returns:
        pd.Series or None: The row containing the ancestry region if found, None otherwise.
    """
    region = msp[(msp['spos'] <= pos) & (msp['epos'] >= pos)]
    
    if len(region) > 1:
        raise ValueError(f"Multiple ancestry regions found for position {pos}. Check MSP file integrity.")

    return region.iloc[0] if not region.empty else None


def get_ancestry(pos_idx, sample_idx, msp, callset):
    """
    Return the ancestry region for a given variant position and sample index.

    Returns:
        tuple(str, str): Ancestry of haplotype 0 and haplotype 1.
    """
    pos_value = callset['variants/POS'][pos_idx]
    region = locate_ancestry_region(pos_value, msp)

    if region is None:
        return "Unknown", "Unknown"

    sample_name = callset['samples'][sample_idx]
    hap_0, hap_1 = f"{sample_name}.0", f"{sample_name}.1"

    if hap_0 not in msp.columns or hap_1 not in msp.columns:
        raise KeyError(f"Haplotype columns '{hap_0}' or '{hap_1}' not found in MSP file.")

    return region[hap_0], region[hap_1]


def gather_test_data():
    """
    Gather data for testing by sampling a variant position, genotype, and ancestry.

    Returns:
        dict: A dictionary containing sampled test data.
    """
    pos_idx, pos_value, sample_idx = sampler(callset)
    sample = callset['samples'][sample_idx]
    alleles = get_alleles(pos_idx, sample_idx)
    hap_0_ancestry, hap_1_ancestry = get_ancestry(pos_idx, sample_idx, msp, callset)

    return {
        "variant_position": pos_value,
        "sample_index": sample_idx,
        "sample_name": sample,
        "alleles": alleles,  # Keeping it as a tuple for easy processing
        "ancestry": (hap_0_ancestry, hap_1_ancestry)
    }


def validate_haplotype_ancestry(ANC, hap_ancestry, genotype):
    """
    Validate whether a given genotype should be encoded as missing based on ancestry.

    Returns:
        bool: True if the genotype is correctly assigned, False if there is a mismatch.
    """
    if ANC != hap_ancestry:
        return genotype == -1  # Expect genotype to be missing (-1)

    return genotype in (0, 1)  # Otherwise, genotype must be valid (0 or 1)


def run_test_iteration(ANC):
    """
    Run a single test iteration to validate ancestry-based genotype masking.
    If the test fails, the script will immediately exit.
    
    Args:
        ANC (str): The expected ancestry to validate against.
    """
    # Step 1: Gather test data
    test_data = gather_test_data()

    # Extract relevant values
    pos_value = test_data["variant_position"]
    sample = test_data["sample_name"]
    sample_idx = test_data["sample_index"]
    alleles = test_data["alleles"]
    hap_0_ancestry, hap_1_ancestry = test_data["ancestry"]

    # Step 2: Validate both haplotypes using the expected ancestry (ANC)
    valid_0 = validate_haplotype_ancestry(ANC, hap_0_ancestry, alleles[0])
    valid_1 = validate_haplotype_ancestry(ANC, hap_1_ancestry, alleles[1])

    # Step 3: Print results
    print("=" * 50)
    print(" TEST ITERATION RESULTS")
    print("=" * 50)
    print(f"  Expected Ancestry   : {ANC}")
    print(f"  Variant Position    : {pos_value}")
    print(f"  Sample              : {sample} (Index: {sample_idx})")
    print(f"  Genotype (GT)       : {alleles[0]} | {alleles[1]}")
    print(f"  Ancestry            : {hap_0_ancestry} | {hap_1_ancestry}")
    print(f"  Validation (GT Masking)")
    print(f"    - Haplotype 0     : {'✔️ Passed' if valid_0 else '❌ Failed'}")
    print(f"    - Haplotype 1     : {'✔️ Passed' if valid_1 else '❌ Failed'}")

    # Step 4: Check for failure and exit if necessary
    if valid_0 and valid_1:
        print("\n✅ SUCCESS: Both haplotypes passed validation ✅")
    else:
        print("\n❌ FAILURE: At least one haplotype failed validation ❌")
        print("Aborting script due to test failure...\n")
        sys.exit(1)  # Exit script with error code 1

    print("=" * 50)



if __name__ == "__main__":
    # Load the VCF dataset and the MSP file
    callset = allel.read_vcf("data/mask_10.vcf.gz")
    msp = pd.read_csv("data/query_results.msp", sep="\t", skiprows=1)

    expected_ancestry = 1

    n_iterations = 1000
    for i in range(n_iterations):
        run_test_iteration(expected_ancestry)
