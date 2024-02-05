rule cluster_97:
    input:
        os.path.join(config['output_relabeled'], "merged_polished_relabeled.fasta")
    output:
        otu_table = os.path.join(config['output_OTU'], "otu_table_97.tsv"),
        otu_centroids = os.path.join(config['output_OTU'], "otu_97.fa")
    threads:
        config['max_threads']
    conda:
        "envs/vsearch.yml"
    log:
        os.path.join(config['log_dir'], "clustering_ID", "otu_97.log")
    shell:
        """
        vsearch \
        --cluster_fast {input} \
        --id 0.97 \
        --threads {threads} \
        --relabel OTU_ \
        --sizeout \
        --otutabout {output.otu_table} \
        --centroids {output.otu_centroids}
        """
rule cluster_99:
    input:
        os.path.join(config['output_relabeled'], "merged_polished_relabeled.fasta")
    output:
        otu_table = os.path.join(config['output_OTU'], "otu_table_99.tsv"),
        otu_centroids = os.path.join(config['output_OTU'], "otu_99.fa")
    threads:
        config['max_threads']
    conda:
        "envs/vsearch.yml"
    log:
        os.path.join(config['log_dir'], "clustering_ID", "otu_99.log")
    shell:
        """
        vsearch \
        --cluster_fast {input} \
        --id 0.99 \
        --threads {threads} \
        --relabel OTU_ \
        --sizeout \
        --otutabout {output.otu_table} \
        --centroids {output.otu_centroids}
        """
