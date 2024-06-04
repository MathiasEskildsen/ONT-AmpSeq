rule relabel:
    input:
        os.path.join(config["output_dir"], "polish", "samples", "{sample}_polished.fasta")
    output:
        os.path.join(config['output_dir'], "polish", "samples_relabeled", "{sample}_relabeled.fasta")
    conda:
        "../envs/vsearch.yml"
    threads:
        1
    resources:
        mem_mb = 2048,
        runtime = 60
    log:
        os.path.join(config['log_dir'], "relabel", "{sample}.log")
    shell:
        """
        vsearch \
            --sortbysize {input} \
            --relabel {wildcards.sample}. \
            --threads {threads} \
            --output {output}
        """

rule relabel_merge:
    input:
        expand(os.path.join(config['output_dir'], "polish", "samples_relabeled", "{sample}_relabeled.fasta"), sample=get_samples())
    output:
        os.path.join(config['output_dir'], "polish", "samples_relabeled", "merged_polished_relabeled.fasta")
    threads:
        1
    resources:
        mem_mb = 512,
        runtime = 60
    conda:
        "../envs/vsearch.yml"
    log:
        os.path.join(config['log_dir'], "relabel", "merged_polished_relabeled.log")
    shell:
        """
        cat {input} > {output}
        """
