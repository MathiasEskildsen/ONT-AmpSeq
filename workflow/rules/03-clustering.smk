configfile: "config/config.yaml"
import os
rule convert_to_fasta:
    input:
        os.path.join(config["tmp_dir"], "samples", "{sample}_filtered.fastq")
    output:
        os.path.join(config["tmp_dir"], "samples", "{sample}_filtered.fasta")
    threads: 1
    resources:
        mem_mb = 512,
        runtime = 60
    conda:
        "../envs/vsearch.yml"
    log:
        os.path.join(config["log_dir"], "convert_to_fasta", "{sample}.log")
    shell:
        """
        sed -n '1~4s/^@/>/p;2~4p' {input} > {output}    
    """
rule vsearch_cluster:
    input:
        os.path.join(config["tmp_dir"], "samples", "{sample}_filtered.fasta")
    output:
        os.path.join(config["output_dir"], "vsearch", "samples", "{sample}_cluster.fasta")
    threads:
        1
    resources:
        mem_mb = 2048,
        runtime = 1440
    conda:
        "../envs/vsearch.yml"
    log:
        os.path.join(config["log_dir"], "vsearch_cluster", "{sample}.log")
    shell:
        """
    vsearch \
        --cluster_unoise {input} \
        --minsize 1 \
        --threads {threads} \
        --centroids {output}
    """


