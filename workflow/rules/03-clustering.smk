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
        runtime=30,
    log:
        os.path.join(
            config["log_dir"], "03-clustering", "convert_to_fasta", "{sample}.log"
        ),
    shell:
        """
        {{
            echo "Starting conversion to FASTA for sample {wildcards.sample}"
            # will this work with interleaved fastq files?
            sed -n '1~4s/^@/>/p;2~4p' {input} > {output}
            echo "Finished conversion to FASTA for sample {wildcards.sample}"
        }} > {log} 2>&1  
    """

rule concatenate_all_samples:
    input:
        expand(
            os.path.join(
                config["tmp_dir"], "samples", "{sample}_filtered.fasta"
            ),
            sample=get_samples(),
        ),
    # Use aggregate rule to concatenate all files using wildcard.sample
    output:
        os.path.join(
            config["tmp_dir"], "samples", "all_samples.fasta"
        ),
    conda:
        "../envs/mapping.yml"
    log:
        os.path.join(
            config["log_dir"], "03-clustering", "concatenate_all_samples.log"
        ),
    resources:
        mem_mb=512,
        runtime=30,
    shell:
        """
        {{
        echo "Starting concatenation of OTUs"
        cat {input} > {output}
        echo "Finished concatenation of OTUs"
        }} > {log} 2>&1
        """

rule concatenated_derep:
    input:
        os.path.join(
            config["tmp_dir"], "samples", "all_samples.fasta"
        ),
    output:
        os.path.join(
            config["output_dir"], "vsearch", "samples", "concatenated_derep.fasta"
        ),
    threads: 1
    resources:
        mem_mb=8192,
        runtime=60,
    conda:
        "../envs/vsearch.yml"
    log:
        os.path.join(
            config["log_dir"], "03-clustering", "concatenated_derep.log"
        ),
    shell:
        """
        {{
            echo "Starting dereplication of concatenated sequences from all samples"
            vsearch --derep_fulllength \
                {input} \
                --sizeout \
                --fasta_width 0 \
                --output {output} \
                --minuniquesize 2
            echo "Finished dereplication of concatenated sequences from all samples"
        }} > {log} 2>&1
    """
    

rule vsearch_denoise:
    input:
        os.path.join(
            config["output_dir"], "vsearch", "samples", "concatenated_derep.fasta"
        ),
    output:
        os.path.join(
            config["output_dir"], "vsearch", "denoised.fasta"
        ),
    threads: config["max_threads"]
    resources:
        mem_mb=8192,
        runtime=1440,
    params:
        unoise_minsize = config['unoise_minsize'],
    conda:
        "../envs/vsearch.yml"
    log:
        os.path.join(
            config["log_dir"], "03-clustering", "vsearch_denoise.log"
        ),
    shell:
        """
        {{
            echo "Starting denoising all reads "
            # unoise is for Illumina reads only
            vsearch \
                --cluster_unoise {input} \
                --minsize {params.unoise_minsize} \
                --threads {threads} \
                --centroids {output} \
                --sizein \
                --sizeout
            echo "Finished denoising all reads"
        }} > {log} 2>&1
    """
