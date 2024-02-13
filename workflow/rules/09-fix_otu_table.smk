rule fix_otu_table:
    input: 
        otu_txt =  os.path.join(config['output_OTU'], "{id}", "otu_taxonomy_{id}.txt"),
        otu_tsv = os.path.join(config['output_OTU'], "{id}", "otu_cluster_{id}.tsv")
    output:
        output_all = os.path.join(config["output_OTU"], "{id}", "otu_table_all_{id}.tsv"),
        output_temp = temp(os.path.join(config["tmp_dir"], "{id}", "otu_taxonomy_{id}_cut_temp.txt"))
    threads:
        1
    resources:
        mem_mb = 1024,
        runtime = "01:00:00"
    log:
        os.path.join(config["log_dir"], "fix_otu_table", "otu_table_all_{id}.log")
    shell:
        """
        awk -F "\t" 'OFS="\t"{{gsub(/;.*/," ",$1);print $1 $4}}' {input.otu_txt} | awk -F ' ' '{{print $1"\t"$2}}' > {output.output_temp}
        awk 'BEGIN {{FS=OFS="\t"}} NR==FNR {{hold[$1]=$2; next}} {{print $0, hold[$1]}}' {output.output_temp} "{input.otu_tsv}" > {output.output_all}
        """