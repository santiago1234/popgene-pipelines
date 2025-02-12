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


rule freeze_samples:
    """
    Extract a vcf containing samples in query and ref
    """
    input:
        config['query_samples'],
        f'{SCRATCH}/ref-samples.txt'
    output:
        f"{SCRATCH}/freeze_samples.txt"
    shell:
        """
        cat {input} | sort | uniq >>{output}
        """

rule freeze_vcf:
    """
    VCF containing both query and ref samples
    """
    input:
        vcf = config['vcf'],
        samples = f"{SCRATCH}/freeze_samples.txt"
    output:
        f"{SCRATCH}/vcf/freeze_chr{{chrn}}.vcf.gz"
    params:
        mac = MAC
    shell:
        """
        if [[ {params.mac} -eq 0 ]]; then
            bcftools view -S {input.samples} {input.vcf} -Oz -o {output}
        else
            bcftools view -c{params.mac} -S {input.samples} {input.vcf} -Oz -o {output}
        fi
        """
        

rule query_vcf:
    input:
        vcf = f"{SCRATCH}/vcf/freeze_chr{{chrn}}.vcf.gz",
        samples = config['query_samples']
    output:
        f"{SCRATCH}/vcf/query_chr{{chrn}}.vcf.gz"
    shell:
        """
        bcftools view -S {input.samples} {input.vcf} -Oz -o {output}
        """


rule ref_vcf:
    input:
        vcf = f"{SCRATCH}/vcf/freeze_chr{{chrn}}.vcf.gz",
        samples = f'{SCRATCH}/ref-samples.txt'
    output:
        f"{SCRATCH}/vcf/reference_chr{{chrn}}.vcf.gz"
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
