rule individual_outputs_blast:
    input:
        abund = os.path.join(config['output_dir'], "cluster", "{id}", "otu_cluster_{id}.tsv"),
        tax = os.path.join(config["output_dir"], "OTU-tables", "{id}", "otu_taxonomy_{id}_blast_trimmed.txt")
    output:
        abund = os.path.join(config['output_dir'], "final", "{id}", "phyloseq_abundance_{id}_blast.tsv"),
        tax = os.path.join(config['output_dir'], "final", "{id}", "phyloseq_tax_{id}_blast.tsv")
    log: 
        os.path.join(config['log_dir'], "phyloseq_tax", 'phyloseq_{id}.log')
    threads:
        1
    resources:
        mem_mb = 1024,
        runtime = "01:00:00"
    conda:
        "../envs/generic.yml"
    shell:
        """
        cp {input.tax} {output.tax}
        cp {input.abund} {output.abund}
        sed -i '1 s/#OTU ID/OTU/' {output.abund}
        """
rule merge_abund_tax_blast:
    input:
        tax = os.path.join(config["output_dir"], "OTU-tables", "{id}", "otu_taxonomy_{id}_blast_trimmed.txt"),
        abund = os.path.join(config['output_dir'], "cluster", "{id}", "otu_cluster_{id}.tsv")
    output:
        os.path.join(config['output_dir'], "final", "{id}", 'OTUtable_tax_{id}_blast.tsv')
    threads:
        1
    resources:
        mem_mb = 1024,
        runtime = "01:00:00"
    log:
        os.path.join(config['log_dir'], "taxonomy_blast", 'merge_abund_tax_blast_{id}.log')
    run:
        merge_abund_tax_blast(input["tax"], input["abund"], output[0])