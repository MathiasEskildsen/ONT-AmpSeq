rule taxonomy_sintax_97:
    input: 
        os.path.join(config['output_OTU'], 'otu_97.fa')
    output:
        os.path.join(config['output_OTU'], 'otu_taxonomy_97.txt')
    threads:
        config['max_threads']
    conda:
        "envs/vsearch.yml"
    params:
        db_path = config['db_path']
    log:
        os.path.join(config['log_dir'], "vsearch_tax", "otu_taxonomy_97.log")
    shell:
        """"
        vsearch \
        --sintax {input} \
        --db {params.db_path} \
        --tabbedout {output} \
        --threads {threads} \
        --sintax_cutoff 0.97 \
        --strand both
        """
rule taxonomy_sintax_99:
    input: 
        os.path.join(config['output_OTU'], 'otu_99.fa')
    output:
        os.path.join(config['output_OTU'], 'otu_taxonomy_99.txt')
    threads:
        config['max_threads']
    conda:
        "envs/vsearch.yml"
    params:
        db_path = config['db_path']
    log:
        os.path.join(config['log_dir'], "vsearch_tax", "otu_taxonomy_99.log")
    shell:
        """"
        vsearch \
        --sintax {input} \
        --db {params.db_path} \
        --tabbedout {output} \
        --threads {threads} \
        --sintax_cutoff 0.99 \
        --strand both
        """

