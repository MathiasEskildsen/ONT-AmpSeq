import glob
import os

rule concatenate_otus:
    input:
        expand(os.path.join(config["output_cluster"], "samples", "{sample}_cluster.fasta"),sample=sample_dirs) # Use aggregate rule to concatenate all files using wildcard.sample
    output:
        os.path.join(config["output_cluster"], "concatenate_otus", "concatenated_otus.fasta")
    resources:
        mem_mb = 1024,
        runtime = "01:00:00" # - Adding run time to the rule. Standard on SLURM config is 1 hour. days-hours:minutes:seconds
    shell:
        """
        cat {input} > {output}
        """

rule mapping:
    input:
        combined = os.path.join(config["output_cluster"], "concatenate_otus", "concatenated_otus.fasta"),
        samples = os.path.join(config["output_cluster"], "samples", "{sample}_cluster.fasta")
    output:
        os.path.join(config["output_mapping"], "mapping", "{sample}_aligned.sam")
    resources:
        mem_mb = 20480,
        runtime = "1-00:00:00" # - Adding run time to the rule. Standard on SLURM config is 1 hour. days-hours:minutes:seconds
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
