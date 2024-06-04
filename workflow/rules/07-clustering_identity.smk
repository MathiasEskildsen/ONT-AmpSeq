rule cluster_ID:
    input:
        os.path.join(config['output_dir'], "polish", "samples_relabeled", "merged_polished_relabeled.fasta")
    output:
        otu_table = os.path.join(config['output_dir'], "cluster", "{id}", "otu_cluster_{id}.tsv"),
        otu_centroids = os.path.join(config['output_dir'], "cluster", "{id}", "otu_{id}.fa")
    threads:
        config['max_threads']
    resources:
        mem_mb = 2048,
        runtime = 1440
    conda:
        "../envs/vsearch.yml"
    log:
        os.path.join(config['log_dir'], "cluster_ID", "otu_{id}.log")
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