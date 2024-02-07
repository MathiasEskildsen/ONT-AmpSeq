rule relabel:
    input:
        os.path.join(config['output_polish'], "{sample}.polished.fasta")
    output:
        os.path.join(config['output_relabeled'], "{sample}.relabeled.fasta")
    conda:
        "../envs/vsearch.yml"
    threads:
        config['max_threads']
    resources:
        mem_mb = 2048
    log:
        os.path.join(config['log_dir'], "relabel" "{sample}.relabeled.log")
    shell:
        """"
        vsearch \
        --sortbysize {input} \
        --relabel {wildcards.sample}. \
        --threads {threads} \
        --output {output}
        """

rule relabel_merge:
    input:
        expand(os.path.join(config['output_relabeled'], "{sample}.relabeled.fasta"), sample=config['input_dir'])
    output:
        os.path.join(config['output_relabeled'], "merged_polished_relabeled.fasta")
    threads:
        config['max_threads']
    resources:
        mem_mb = 512
    log:
        os.path.join(config['log_dir'], "relabel" ,"merged_polished_relabeled.log")
    shell:
        """"
        cat {input} > {output}
        """
