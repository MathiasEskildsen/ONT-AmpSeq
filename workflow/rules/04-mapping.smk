import glob
import os

rule concatenate_otus:
    input:
        expand(os.path.join(config["output_dir"], "vsearch", "samples", "{sample}_cluster.fasta"),sample=get_samples()) # Use aggregate rule to concatenate all files using wildcard.sample
    output:
        os.path.join(config["output_dir"], "vsearch", "samples", "concatenated_otus.fasta")
    conda:
        "../envs/mapping.yml"
    log:
        os.path.join(config["log_dir"], "concat_otus","concatenate_otus.log")
    resources:
        mem_mb = 1024,
        runtime = 60
    shell:
        """
        cat {input} > {output}
        """

rule mapping:
    input:
        combined = os.path.join(config["output_dir"], "vsearch", "samples", "concatenated_otus.fasta"),
        samples = os.path.join(config["output_dir"], "vsearch", "samples", "{sample}_cluster.fasta")
    output:
        os.path.join(config["output_dir"], "mapping", "samples", "{sample}_aligned.sam")
    resources:
        mem_mb = 40960,
        runtime = 2880
    threads:
        config['max_threads']
    conda:
        "../envs/mapping.yml"
    log:
        os.path.join(config["log_dir"], "mapping", "{sample}.log")
    shell:
        """
        minimap2 \
        -ax map-ont \
        -K500M \
        -t {threads} \
        --secondary=no \
        {input.samples} \
        {input.combined} \
        > {output}
    """
