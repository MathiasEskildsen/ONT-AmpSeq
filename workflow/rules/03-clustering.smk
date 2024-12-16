configfile: "config/config.yaml"


import os


rule convert_to_fasta:
    input:
        os.path.join(config["tmp_dir"], "samples", "{sample}_filtered.fastq"),
    output:
        os.path.join(config["tmp_dir"], "samples", "{sample}_filtered.fasta"),
    threads: 1
    resources:
        mem_mb=512,
        runtime=60,
    conda:
        "../envs/vsearch.yml"
    log:
        os.path.join(
            config["log_dir"], "03-clustering", "convert_to_fasta", "{sample}.log"
        ),
    shell:
        """
        {{
            echo "Starting conversion to FASTA for sample {wildcards.sample}"
            sed -n '1~4s/^@/>/p;2~4p' {input} > {output}
            echo "Finished conversion to FASTA for sample {wildcards.sample}"
        }} > {log} 2>&1  
    """


rule vsearch_cluster:
    input:
        os.path.join(config["tmp_dir"], "samples", "{sample}_filtered.fasta"),
    output:
        os.path.join(
            config["output_dir"], "vsearch", "samples", "{sample}_cluster.fasta"
        ),
    threads: config["max_threads"]
    resources:
        mem_mb=8192,
        runtime=1440,
    conda:
        "../envs/vsearch.yml"
    log:
        os.path.join(
            config["log_dir"], "03-clustering", "vsearch_cluster", "{sample}.log"
        ),
    shell:
        """
        {{
            echo "Starting clustering for sample {wildcards.sample}"
            vsearch \
                --cluster_unoise {input} \
                --minsize 1 \
                --threads {threads} \
                --centroids {output}
            echo "Finished clustering for sample {wildcards.sample}"
        }} > {log} 2>&1
    """
