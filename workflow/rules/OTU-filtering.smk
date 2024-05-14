configfile: "config/config.yml"
rule filter_otu_tables:
    input:
        OTU_sintax_raw  =   os.path.join(config["output_dir"], "final", "{id}", "OTUtable_tax_{id}_sintax.tsv"),
        OTU_blast_raw   =   os.path.join(config["output_dir"], "final", "{id}", "OTUtable_tax_{id}_blast.tsv")
    output:
        OTU_sintax_1    =   os.path.join(config["output_dir"], "final", "{id}", "OTUtable_tax_{id}_sintax_1.tsv"),
        OTU_sintax_5    =   os.path.join(config["output_dir"], "final", "{id}", "OTUtable_tax_{id}_sintax_5.tsv"),
        OTU_sintax_10   =   os.path.join(config["output_dir"], "final", "{id}", "OTUtable_tax_{id}_sintax_10.tsv"),
        OTU_blast_1 =   os.path.join(config["output_dir"], "final", "{id}", "OTUtable_tax_{id}_blast_1.tsv"),
        OTU_blast_5 =   os.path.join(config["output_dir"], "final", "{id}", "OTUtable_tax_{id}_blast_5.tsv"),
        OTU_blast_10    =   os.path.join(config["output_dir"], "final", "{id}", "OTUtable_tax_{id}_blast_10.tsv")
    threads:
        1
    resources:
        mem_mb = 1024
        time = "01:00:00"
    shell:
        """
        bash 