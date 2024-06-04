include: "common.smk"
rule phyloseq_abund_sintax:
    input:
        os.path.join(config['output_dir'], "cluster", "{id}", "otu_cluster_{id}.tsv")
    output:
        phyloseq_abundance = os.path.join(config['output_dir'], "final", "{id}", 'phyloseq_abundance_{id}_sintax.tsv')
    log: 
        os.path.join(config['log_dir'], "phyloseq_abundance", 'phyloseq_{id}.log')
    threads:
        1
    conda:
        "../envs/generic.yml"
    resources:
        mem_mb = 1024,
        runtime = 60
    shell:
        """
        cp {input} {output}
        sed -i '1 s/#OTU ID/OTU/' {output}
        """

rule phyloseq_sintax:
    input:
        input_OTU = os.path.join(config['output_dir'], "final", "{id}", 'OTUtable_tax_{id}_sintax.tsv')
    output:
        phyloseq_tax = os.path.join(config['output_dir'], "final", "{id}", 'phyloseq_tax_{id}_sintax.tsv'),
        tmp_output = temp(os.path.join(config['tmp_dir'], "phyloseq_tax_{id}_tmp.tsv"))
    log:
        os.path.join(config['log_dir'], "phyloseq_tax", 'phyloseq_{id}.log')
    threads:
        1
    resources:
        mem_mb = 1024,
        runtime = 60
    conda:
        "../envs/generic.yml"
    script:
        "../scripts/phyloseq_tax_sintax.py"