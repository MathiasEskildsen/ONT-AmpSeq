import glob
import os


rule concatenate_otus:
    input:
        expand(
            os.path.join(
                config["output_dir"], "vsearch", "samples", "{sample}_cluster.fasta"
            ),
            sample=get_samples(),
        ),
    # Use aggregate rule to concatenate all files using wildcard.sample
    output:
        os.path.join(
            config["output_dir"], "vsearch", "samples", "concatenated_otus.fasta"
        ),
    conda:
        "../envs/mapping.yml"
    log:
        os.path.join(
            config["log_dir"], "04-mapping", "concatenate_otus", "concatenate_otus.log"
        ),
    resources:
        mem_mb=1024,
        runtime=60,
    shell:
        """
        {{
        echo "Starting concatenation of OTUs"
        cat {input} > {output}
        echo "Finished concatenation of OTUs"
        }} > {log} 2>&1
        """


rule mapping:
    input:
        combined=os.path.join(
            config["output_dir"], "vsearch", "samples", "concatenated_otus.fasta"
        ),
        samples=os.path.join(
            config["output_dir"], "vsearch", "samples", "{sample}_cluster.fasta"
        ),
    output:
        os.path.join(config["output_dir"], "mapping", "samples", "{sample}_aligned.sam"),
    resources:
        mem_mb = 40960,
        runtime = 2880
    params:
        K = config['K'],
        f = config['f']
    threads:
        config['max_threads']
    conda:
        "../envs/mapping.yml"
    log:
        os.path.join(config["log_dir"], "04-mapping", "mapping", "{sample}.log"),
    shell:
        """
        {{
        minimap2 \
        -ax map-ont \
        -K{params.K} \
        -f {params.f} \
        -t {threads} \
        --secondary=no \
        {input.samples} \
        {input.combined} \
        > {output}
        }} > {log} 2>&1
    """
