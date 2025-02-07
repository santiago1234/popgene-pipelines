import click
import numpy as np
import allel
from cyvcf2 import VCF, Writer


def vcf_positions(vcf_file):
    """Extracts variant positions from a VCF file."""
    try:
        callset = allel.read_vcf(vcf_file, fields=['variants/CHROM', 'variants/POS'])
        if not callset or 'variants/POS' not in callset:
            raise ValueError("VCF file is empty or improperly formatted.")

        chrom_set = set(callset['variants/CHROM'])
        if len(chrom_set) != 1:
            raise ValueError(f"VCF contains multiple chromosomes: {chrom_set}. Only one chromosome is supported.")

        return np.sort(callset['variants/POS']), chrom_set.pop()
    except Exception as e:
        raise RuntimeError(f"Error reading VCF positions: {e}")


def split_positions(positions, n_split, chrom):
    """Splits positions into N equal-sized regions."""
    if n_split <= 0:
        raise ValueError("Number of splits must be greater than 0.")
    if len(positions) < n_split:
        raise ValueError(f"VCF file has fewer variants ({len(positions)}) than requested splits ({n_split}).")

    pos_splits = np.array_split(positions, n_split)
    return [f'{chrom}:{ps[0]}-{ps[-1]}' for ps in pos_splits]


def write_to_vcf_files(vcf_file, regions, out_file_prefix):
    """Writes each region to a separate VCF file."""
    vcf = VCF(vcf_file)  # Open VCF only once to avoid repeated I/O operations

    for i, region in enumerate(regions):
        out_vcf = f'{out_file_prefix}_{i+1}.vcf.gz'
        print(f'Writing {region} to {out_vcf}...')

        try:
            writer = Writer(out_vcf, vcf)
            for variant in vcf(region):
                writer.write_record(variant)
            writer.close()
        except Exception as e:
            print(f"Error writing to {out_vcf}: {e}")
        finally:
            writer.close()


@click.command()
@click.option('--vcf_file', required=True, type=click.Path(exists=True), help='Path to the VCF file.')
@click.option('--out_file_prefix', required=True, help='Prefix for the output VCF files.')
@click.option('--n_split', required=True, type=int, help='Number of splits.')
def main(vcf_file, out_file_prefix, n_split):
    """Main function to split a VCF file into multiple smaller files."""
    try:
        positions, chrom = vcf_positions(vcf_file)
        regions = split_positions(positions, n_split, chrom)
        write_to_vcf_files(vcf_file, regions, out_file_prefix)
    except Exception as e:
        print(f"Error: {e}")


if __name__ == '__main__':
    main()
