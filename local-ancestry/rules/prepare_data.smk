rule ref_samples:
    """
    Get a list of the reference samples
    """
    input:
        config['sample_map_file'] 
    output:
        f'{SCRATCH}/ref-samples.txt'
    shell:
        """
        tail -n +2 {input} |\
            cut -f1 >{output}
        """


rule query_vcf:
    input:
        vcf = config['vcf'],
        samples = config['query_samples']
    output:
        f'{SCRATCH}/vcf/query_chr{{chrn}}.vcf.gz'
    shell:
        """
        bcftools view -S {input.samples} {input.vcf} -Oz -o {output}
        """


rule ref_vcf:
    input:
        vcf = config['vcf'],
        samples = f'{SCRATCH}/ref-samples.txt'
    output:
        f'{SCRATCH}/vcf/reference_chr{{chrn}}.vcf.gz'
    shell:
        """
        bcftools view -S {input.samples} {input.vcf} -Oz -o {output}
        """


rule format_genetic_map:
    '''
    Put the genetic map in the format that
    gnomix needs:
        It's a .tsv file with 3 columns; chromosome number, SNP physical
        position and SNP genetic position. There should be no headers unless
        they start with "#".""'
    '''
    input:
        config['genetic_map']
    output:
        f'{SCRATCH}/genetic_maps/chr{{chrn}}.map'
    shell:
        """
        awk '{{ print $1 "\t" $4 "\t" $3 }}' {input} >{output}
        """
