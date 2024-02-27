configfile: "config/config.yaml"

if config['taxonomy_sintax']:
    rule taxonomy_sintax:
        input: 
            os.path.join(config['output_OTU'], "{id}", "otu_{id}.fa")
        output:
            os.path.join(config['output_OTU'], "{id}", "otu_taxonomy_{id}_sintax.txt")
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
if config['taxonomy_blast']:
    rule taxonomy_blast:
        input: 
            os.path.join(config['output_OTU'], "{id}", "otu_{id}.fa")
        output:
            os.path.join(config['output_OTU'], "{id}", "otu_taxonomy_{id}_blast.txt")
        threads:
            config['max_threads']
        resources:
            mem_mb = 2048,
            runtime = "02:00:00"
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
                -query {input} \
                -word_size 11 \
                -max_taget_seqs 1 \
                -num_threads {threads} \
                -outfmt "6 qseqid sseqid stitle evalue bitscore length pident" \
                -out {output} \
                -db {params.db_path}
            """