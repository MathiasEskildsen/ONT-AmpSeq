rule polish_racon:
    input:
        combined = os.path.join(config["output_dir"], "vsearch", "samples", "concatenated_otus.fasta"),
        alignment = os.path.join(config["output_dir"], "mapping", "samples", "{sample}_aligned.sam"),
        polish_target = os.path.join(config["output_dir"], "vsearch", "samples", "{sample}_cluster.fasta")
    output:
        os.path.join(config["output_dir"], "polish", "samples", "{sample}_polished.fasta")
    conda:
        "../envs/polish.yml"
    threads:
        config['max_threads']
    resources:
        mem_mb = 10240,
        runtime = 60
    log:
        os.path.join(config["log_dir"], "05-polish_racon", "{sample}.log")
    shell:
        """
        {{
        racon \
        {input.combined} \
        {input.alignment} \
        {input.polish_target} \
        -t {threads} \
        > {output}
        }} > {log} 2>&1
        """