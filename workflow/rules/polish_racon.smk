rule polish_racon:
    input:
        combined = os.path.join(config["output_cluster"], "concatenated_otus.fasta"),
        alignment = os.path.join(config["output_mapping"], "{sample}_aligned.sam"),
        polish_target = os.path.join(config["output_cluster"], "{sample}.cluster.fasta")
    output:
        os.path.join(config["output_polish"], "{sample}.polished.fasta")
    conda:
        "../envs/racon.yml"
    threads:
        config['max_threads']
    log:
        os.path.join(config["log_dir"], "polish", "{sample}.polish_racon.log")
    shell:
        """
        racon \
        {input.combined} \
        {input.alignment} \
        {input.polish_target} \
        -t {threads} \
        > {output}
        """