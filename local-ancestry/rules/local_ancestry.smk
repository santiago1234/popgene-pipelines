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
        output_folder = "results"
    output:
         "train_{chrn}.txt"
    threads: 10
    log: 'logs/train_{chrn}.log'
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
            {input.config_file} &>{log}

        touch {output}
        """
