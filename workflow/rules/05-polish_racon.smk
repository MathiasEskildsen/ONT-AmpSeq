rule polish_racon:
    input:
        combined = os.path.join(config["output_cluster"], "concatenate_otus", "concatenated_otus.fasta"),
        alignment = os.path.join(config["output_mapping"], "mapping", "{sample}_aligned.sam"),
        polish_target = os.path.join(config["output_cluster"], "samples", "{sample}_cluster.fasta")
    output:
        os.path.join(config["output_polish"], "samples", "{sample}_polished.fasta")
    conda:
        "../envs/polish.yml"
    threads:
        config['max_threads']
    resources:
        mem_mb = 1024, # 20GB might need to be adjusted to appropriate value after testing
        runtime = "01:00:00"
    log:
        os.path.join(config["log_dir"], "polish_racon", "{sample}.log")
    shell:
        """
        racon \
        {input.combined} \
        {input.alignment} \
        {input.polish_target} \
        -t {threads} \
        > {output}
        """