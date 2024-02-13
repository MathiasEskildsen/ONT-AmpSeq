rule taxonomy_sintax:
    input: 
        os.path.join(config['output_OTU'], "{id}", "otu_{id}.fa")
    output:
        os.path.join(config['output_OTU'], "{id}", "otu_taxonomy_{id}.txt")
    threads:
        config['max_threads']
    resources:
        mem_mb = 2048,
        runtime = "02:00:00"
    conda:
        "../envs/vsearch.yml"
    params:
        db_path = config['db_path'],
        id = lambda wildcards: float(wildcards.id) / 100 # <- sintax cutoff (bootstrap value)
    log:
        os.path.join(config['log_dir'], "taxonomy_sintax", "{id}", "otu_taxonomy_{id}.log")
    shell:
        """
        vsearch \
            --sintax {input} \
            --db {params.db_path} \
            --tabbedout {output} \
            --threads {threads} \
            --sintax_cutoff {params.id} \
            --strand both
        """
