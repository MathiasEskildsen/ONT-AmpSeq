import glob
import os

rule mapping:
    input:
        combined=os.path.join(
            config["output_dir"], "vsearch", "denoised.fasta"
        ),
        samples=os.path.join(
            config["tmp_dir"], "samples", "{sample}_filtered.fasta"
        ),
    output:
        touch(os.path.join(config["output_dir"], "mapping", "samples", "{sample}_aligned.sam")),
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
