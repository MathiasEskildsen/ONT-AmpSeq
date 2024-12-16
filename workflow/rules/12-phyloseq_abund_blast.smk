rule individual_outputs_blast:
    input:
        abund=os.path.join(
            config["output_dir"], "cluster", "{id}", "otu_cluster_{id}.tsv"
        ),
        tax=os.path.join(
            config["output_dir"],
            "OTU-tables",
            "{id}",
            "otu_taxonomy_{id}_blast_trimmed.txt",
        ),
    output:
        abund=os.path.join(
            config["output_dir"], "final", "{id}", "phyloseq_abundance_{id}_blast.tsv"
        ),
        tax=os.path.join(
            config["output_dir"], "final", "{id}", "phyloseq_tax_{id}_blast.tsv"
        ),
    log:
        os.path.join(
            config["log_dir"],
            "12-phyloseq_abund_blast",
            "individual_outputs_blast",
            "phyloseq_{id}.log",
        ),
    threads: 1
    resources:
        mem_mb=1024,
        runtime=60,
    conda:
        "../envs/generic.yml"
    shell:
        """
        {{
        echo "removing "#" from the OTU ID"
        cp {input.tax} {output.tax}
        cp {input.abund} {output.abund}
        sed -i '1 s/#OTU ID/OTU/' {output.abund}
        echo "Finished the process"
        }} > {log} 2>&1
        """


rule merge_abund_tax_blast:
    input:
        tax=os.path.join(
            config["output_dir"],
            "OTU-tables",
            "{id}",
            "otu_taxonomy_{id}_blast_trimmed.txt",
        ),
        abund=os.path.join(
            config["output_dir"], "cluster", "{id}", "otu_cluster_{id}.tsv"
        ),
    output:
        os.path.join(
            config["output_dir"], "final", "{id}", "OTUtable_tax_{id}_blast_full.tsv"
        ),
    threads: 1
    resources:
        mem_mb=1024,
        runtime=60,
    log:
        os.path.join(
            config["log_dir"],
            "12-phyloseq_abund_blast",
            "merge_abund_tax_blast",
            "merge_abund_tax_blast_{id}.log",
        ),
    run:
        merge_abund_tax_blast(input["tax"], input["abund"], output[0])


rule fix_tax_blast:
    input:
        os.path.join(
            config["output_dir"], "final", "{id}", "OTUtable_tax_{id}_blast_full.tsv"
        ),
    output:
        os.path.join(
            config["output_dir"], "final", "{id}", "OTUtable_tax_{id}_blast.tsv"
        ),
    log:
        os.path.join(
            config["log_dir"],
            "12-phyloseq_abund_blast",
            "fix_tax_blast",
            "fix_tax_blast_{id}.log",
        ),
    threads: 1
    resources:
        mem_mb=1024,
        runtime=60,
    conda:
        "../envs/generic.yml"
    script:
        "../scripts/fix_tax_blast.py"
