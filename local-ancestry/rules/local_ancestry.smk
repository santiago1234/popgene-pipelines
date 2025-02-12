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
    log:
         "logs/train_{chrn}.log"
    output:
        f"{SCRATCH}/results/chr{{chrn}}/models/model_chm_chr{{chrn}}/model_chm_chr{{chrn}}.pkl",
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
            {input.config_file} \
            2>&1 | tee {log}
        """


rule mv_files:
    """
    Put gnomix output data in a nicer formating
    """
    input:
        f"{SCRATCH}/results/chr{{chrn}}/models/model_chm_chr{{chrn}}/model_chm_chr{{chrn}}.pkl",
    params:
        phased_vcf = f"{SCRATCH}/results/chr{{chrn}}/query_file_phased.vcf",
        msp = f"{SCRATCH}/results/chr{{chrn}}/query_results.msp",
        fb = f"{SCRATCH}/results/chr{{chrn}}/query_results.fb",
        model = f"{SCRATCH}/results/chr{{chrn}}/models/model_chm_chr{{chrn}}/model_chm_chr{{chrn}}.pkl",
        bed_dir = directory(f'{SCRATCH}/results/chr{{chrn}}/query_results_bed'), # this is a folder
    output:
        phased_vcf = f"{RESULTS}/query_phased_chr{{chrn}}.vcf.gz",
        index = f"{RESULTS}/query_phased_chr{{chrn}}.vcf.gz.tbi",
        msp = f"{RESULTS}/query_results_chr{{chrn}}.msp",
        fb = f"{RESULTS}/query_results_chr{{chrn}}.fb",
        model = f"{RESULTS}/model_chr{{chrn}}.pkl",
        bed_dir = temp(directory(f'{RESULTS}/bed_chr{{chrn}}')), # this is a folder
    shell:
        """
        bcftools view {params.phased_vcf} -W=tbi -Oz -o {output.phased_vcf}
        cp {params.msp} {output.msp}
        cp {params.fb} {output.fb}
        cp {params.model} {output.model}
        cp -r {params.bed_dir} {output.bed_dir}
        """


rule all_predictions:
    input:
        expand(f'{RESULTS}/' + 'bed_chr{chrn}', chrn=CHROMS),
    output:
        directory(f'{RESULTS}/predictions_bed')
    shell:
        """
        python scripts/collect_beds.py {input} {output}
        """
