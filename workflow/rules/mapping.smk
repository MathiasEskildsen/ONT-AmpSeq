import glob
import os
# Create an input function to get all the fasta files in the output_cluster directory
def get_fasta_files():
    fasta_files = glob.glob(os.path.join(config["output_cluster"], "*.cluster.fasta"))
    return fasta_files

rule concatenate_otus:
    input:
        get_fasta_files()
    output:
        os.path.join(config["output_cluster"], "concatenated_otus.fasta")
    shell:
        """
        cat {input} > {output}
        """

rule mapping:
    input:
        combined = os.path.join(config["output_cluster"], "concatenated_otus.fasta"),
        samples = os.path.join(config["output_cluster"], "{sample}.cluster.fasta")
    output:
        os.path.join(config["output_mapping"], "{sample}_aligned.sam")
    resources:
        mem_mb = 10240
    threads:
        config['max_threads']
    conda:
        "../envs/mapping.yml"
    log:
        os.path.join(config["log_dir"], "aligned", "{sample}.log")
    shell:
        """
        minimap2 \
        -ax map-ont \
        -K500M
        -t {threads} \
        --secondary=no \
        {input.samples} \
        {input.combined} \
        > {output}
    """
