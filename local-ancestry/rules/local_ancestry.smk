rule train_model:
    input:
        query_file = f'{SCRATCH}/vcf/query_chr{{chrn}}.vcf.gz',
        genetic_map_file = f'{SCRATCH}/genetic_maps/chr{{chrn}}.map',
        reference_file = f'{SCRATCH}/vcf/reference_chr{{chrn}}.vcf.gz',
        sample_map_file = config['sample_map_file'],
        config_file = config['gnomix_config']
    params:
        chr_nr = "chr{chrn}",
        phase = config['phase'],
        output_folder = f"{SCRATCH}/results/chr{{chrn}}"
    output:
         "logs/train_{chrn}.log"
    threads: 10
    shell:
        """
        python bin/gnomix/gnomix.py \
            {input.query_file} \
            {params.output_folder} \
            {params.chr_nr} \
            {params.phase} \
            {input.genetic_map_file} \
            {input.reference_file} \
            {input.sample_map_file} \
            {input.config_file} &>{output}
        """


rule mv_files:
    """
    Put gnomix output data in a nicer formating
    """
    input:
         "logs/train_{chrn}.log"
    params:
        phased_vcf = f"{SCRATCH}/results/chr{{chrn}}/query_file_phased.vcf"
    output:
        phased_vcf = f"{RESULTS}/query_phased_chr{{chrn}}.vcf.gz",
        index = f"{RESULTS}/query_phased_chr{{chrn}}.vcf.gz.tbi"
    shell:
        """
        bcftools view {params.phased_vcf} -W=tbi -Oz -o {output.phased_vcf}
        """
        
