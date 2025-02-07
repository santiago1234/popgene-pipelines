"""
Authors: S.G.M.M & C.B.J
Mask VCF by local ancestry. Masked genotypes are encoded as missing values.
Usage:
    python mask-vcf-by-ancestry.py <ANC> <vcf> <lai> <output>

Args:
    ANC: Ancestry code, see the first line of input <lai> file.
    vcf: VCF file
    lai: Local ancestry file (gnomix *.msp output)
    output: output name for masked vcf file
    
NOTES:

    * -1 indicates a missing genotype
    * See imports for required pacakges
"""
import sys

import pandas as pd
import numpy as np
from cyvcf2 import VCF, Writer
import allel

ANC, vcf, lai, output = sys.argv[1:]



def get_vcf_range(vcf_file):
    """
    Get the range of variants
    present in the vcf file
    """
    callset = allel.read_vcf(vcf_file, fields=['variants/CHROM', 'variants/POS'])
    min_position = callset['variants/POS'].min()
    max_position = callset['variants/POS'].max()

    chromosomes = set(callset['variants/CHROM'])
    assert len(chromosomes) == 1, "There should be one and only one chromosome in the vcf file"

    return min_position, max_position, chromosomes.pop()


def load_and_restrict_msp(msp_file, vcf_file):
    """
    Load the vcf file and restrict it to the
    range of variants present in the vcf
    """
    min_position, max_position, vcf_chromosome = get_vcf_range(vcf_file)
    ancestry_data = pd.read_csv(msp_file, sep='\t', skiprows=[0])

    chromosomes_in_ancestry = set(ancestry_data['#chm'].unique())

    assert len(chromosomes_in_ancestry) == 1, "There should be one and only one chromosome in the msp file"
    msp_chromosome = chromosomes_in_ancestry.pop()

    assert str(vcf_chromosome) == str(msp_chromosome), "Chromosome in the vcf file and msp file should be the same"

    ranges = [ancestry_data[(ancestry_data['spos'] <= pos) & (ancestry_data['epos'] >= pos)] for pos in [min_position, max_position]]
    ranges.append(ancestry_data[(ancestry_data['spos'] >= min_position) & (ancestry_data['epos'] <= max_position)])

    ancestry_range = pd.concat(ranges).drop_duplicates()

    return ancestry_range


lai = load_and_restrict_msp(lai, vcf)

vcf = VCF(vcf)
index_to_sample_name = {x: y for (x, y) in enumerate(vcf.samples)}


def retrieve_lai_at(pos, chrom):
    """
    Get the local ancestry of the given range
    """

    lai_pos = lai[(lai.spos <= pos) & (lai.epos >= pos)]
    lai_pos = lai_pos[lai_pos['#chm'].map(str) == str(chrom)]

    if lai_pos.shape[0] < 1:
        raise Exception(f'position {pos} not in local ancestry range')

    if lai_pos.shape[0] > 1:
        print(lai_pos)
        raise Exception(f'position {pos} is in multiple ranges')

    return lai_pos


def sample_ancestry_at(lai_pos, sample_index):
    """
    Get the ancestry of the given sample (index in the vcf)
    at lai_pos
    """
    sample_hap_0 = f'{index_to_sample_name[sample_index]}.0'
    sample_hap_1 = f'{index_to_sample_name[sample_index]}.1'

    anc_hap_0 = lai_pos[sample_hap_0].to_list()[0]
    anc_hap_1 = lai_pos[sample_hap_1].to_list()[0]

    return (anc_hap_0, anc_hap_1)


w = Writer(output, vcf)

i = 0
for variant in vcf:
    # c_: current
    c_pos = variant.POS
    chrom = variant.CHROM
    c_lai = retrieve_lai_at(c_pos, chrom)
    
    if (i % 10000 == 0):
        print(f'masking at position {chrom}-{c_pos} ...')

    i += 1

    for (sample_index, _) in enumerate(variant.genotypes):

        anc_h0, anc_h1 = sample_ancestry_at(c_lai, sample_index)
       # print(f'sample: hap0 = {anc_h0}, hap1 = {anc_h1}')
	# Conver to str to have avoid python reading one as numeric and other as character and showing false when true
        if str(anc_h0) != str(ANC):
          #  print(f'{anc_h0} = {ANC} ?')
            variant.genotypes[sample_index][0] = -1

        if str(anc_h1) != str(ANC):
           # print(f'{anc_h1} = {ANC} ?')
            variant.genotypes[sample_index][1] = -1

    variant.genotypes = variant.genotypes
    w.write_record(variant)

w.close()
vcf.close()

