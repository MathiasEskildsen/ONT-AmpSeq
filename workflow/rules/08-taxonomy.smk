configfile: "config/config.yaml"
rule taxonomy_sintax:
    input: 
        os.path.join(config['output_dir'], "cluster", "{id}", "otu_{id}.fa")
    output:
        os.path.join(config['output_dir'], "taxonomy", "{id}", "otu_taxonomy_{id}_sintax.txt")
    threads:
        config['max_threads']
    resources:
        mem_mb = 2048,
        runtime = "02:00:00"
    conda:
        "../envs/vsearch.yml"
    params:
        db_path = config['db_path_sintax'],
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
rule taxonomy_blast:
    input: 
        os.path.join(config['output_dir'], "cluster", "{id}", "otu_{id}.fa")
    output:
        os.path.join(config['output_dir'], "taxonomy", "{id}", "otu_taxonomy_{id}_blast.txt")
    threads:
        config['max_threads']
    resources:
        mem_mb = 2048,
        runtime = "1-00:00:00"
    conda:
        "../envs/blast.yml"
    params:
        db_path = config['db_path_blast'],
        e_value = config['evalue']
    log:
        os.path.join(config['log_dir'], "taxonomy_blast", "{id}", "otu_taxonomy_{id}.log")
    shell:
        """
        blastn \
            -evalue {params.e_value} \
            -outfmt 6 \
            -word_size 11 \
            -max_target_seqs 1 \
            -db {params.db_path} \
            -num_threads {threads} \
            -query {input} \
            -out {output}
        """