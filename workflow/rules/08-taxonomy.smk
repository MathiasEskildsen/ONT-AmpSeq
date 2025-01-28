configfile: "config/config.yaml"


rule taxonomy_sintax:
    input:
        os.path.join(config["output_dir"], "cluster", "{id}", "otu_{id}.fa"),
    output:
        os.path.join(
            config["output_dir"], "taxonomy", "{id}", "otu_taxonomy_{id}_sintax.txt"
        ),
    threads: config["max_threads"]
    resources:
        mem_mb=2048,
        runtime=120,
    conda:
        "../envs/vsearch.yml"
    params:
        db_path=config["db_path_sintax"],
        id=lambda wildcards: float(wildcards.id) / 100,  # <- sintax cutoff (bootstrap value)
    log:
        os.path.join(
            config["log_dir"],
            "08-taxonomy",
            "taxonomy_sintax",
            "otu_taxonomy_{id}.log",
        ),
    shell:
        """
        {{
        vsearch \
            --sintax {input} \
            --db {params.db_path} \
            --tabbedout {output} \
            --threads {threads} \
            --sintax_cutoff {params.id} \
            --strand both
        }} > {log} 2>&1
        """


rule taxonomy_blast:
    input:
        os.path.join(config["output_dir"], "cluster", "{id}", "otu_{id}.fa"),
    output:
        tax_raw=os.path.join(
            config["output_dir"], "taxonomy", "{id}", "otu_taxonomy_{id}_blast.txt"
        ),
        tax_uniq=os.path.join(
            config["output_dir"],
            "taxonomy",
            "{id}",
            "otu_taxonomy_{id}_blast_uniq.txt",
        ),
    threads: config["max_threads"]
    resources:
        mem_mb=40960,
        runtime=10080,
    conda:
        "../envs/blast.yml"
    params:
        db_path=config["db_path_blast"],
        e_value=config["evalue"],
    log:
        os.path.join(
            config["log_dir"], "08-taxonomy", "taxonomy_blast", "otu_taxonomy_{id}.log"
        ),
    shell:
        """
        {{
        blastn \
            -evalue {params.e_value} \
            -outfmt "6 qseqid stitle evalue bitscore length pident" \
            -word_size 11 \
            -max_target_seqs 1 \
            -db {params.db_path} \
            -num_threads {threads} \
            -query {input} \
            -out {output.tax_raw}
        
        awk '!seen[$0]++' {output.tax_raw} > {output.tax_uniq}
        }} > {log} 2>&1
        """
rule without_taxonomy:
    input: 
        table = os.path.join(config['output_dir'], "cluster", "{id}", "otu_cluster_{id}.tsv"),
        fasta = os.path.join(config['output_dir'], "cluster", "{id}", "otu_{id}.fa"),
    output:
        table = os.path.join(config['output_dir'], "final", "{id}", "OTU_table_{id}.tsv"),
        fasta = os.path.join(config['output_dir'], "final", "{id}", "OTUs_{id}.fna"),
    threads:
        1
    resources:
        mem_mb = 512,
        runtime = 60
    conda:
        "../envs/vsearch.yml"
    log:
        os.path.join(config['log_dir'], "08-taxonomy", "without_taxonomy", "otu_taxonomy_{id}.log")
    shell:
        """
        {{
        cp {input.table} {output.table}
        cp {input.fasta} {output.fasta}
        }} > {log} 2>&1
        """
