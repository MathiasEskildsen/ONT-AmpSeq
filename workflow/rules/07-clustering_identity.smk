rule cluster_ID:
    input:
        os.path.join(config['output_relabeled'], "merged_polished_relabeled.fasta")
    output:
        otu_table = os.path.join(config['output_OTU'], "{id}", "otu_cluster_{id}.tsv"),
        otu_centroids = os.path.join(config['output_OTU'], "{id}", "otu_{id}.fa")
    threads:
        config['max_threads']
    resources:
        mem_mb = 2048
    conda:
        "../envs/vsearch.yml"
    log:
        os.path.join(config['log_dir'], "{id}", "otu_{id}.log")
    params:
        id = lambda wildcards: float(wildcards.id) / 100
    shell:
        """
        vsearch \
        --cluster_fast {input} \
        --id {params.id} \
        --threads {threads} \
        --relabel OTU_ \
        --sizeout \
        --otutabout {output.otu_table} \
        --centroids {output.otu_centroids}
        """