configfile: "config/config.yaml"
if config['taxonomy_sintax']:
    rule fix_otu_table_sintax:
        input: 
            otu_txt =  os.path.join(config['output_OTU'], "{id}", "otu_taxonomy_{id}_sintax.txt"),
            otu_tsv = os.path.join(config['output_OTU'], "{id}", "otu_cluster_{id}.tsv")
        output:
            output_all = os.path.join(config["output_OTU"], "{id}", "otu_table_all_{id}_sintax.tsv"),
            output_temp = temp(os.path.join(config["tmp_dir"], "{id}", "otu_taxonomy_{id}_cut_temp_sintax.txt"))
        threads:
            1
        resources:
            mem_mb = 1024,
            runtime = "01:00:00"
        log:
            os.path.join(config["log_dir"], "fix_otu_table_sintax", "otu_table_all_{id}.log")
        shell:
            """
                awk -F "\t" 'OFS="\t"{{gsub(/;.*/," ",$1);print $1 $4}}' {input.otu_txt} | awk -F ' ' '{{print $1"\t"$2}}' > {output.output_temp}
                awk 'BEGIN {{FS=OFS="\t"}} NR==FNR {{hold[$1]=$2; next}} {{print $0, hold[$1]}}' {output.output_temp} "{input.otu_tsv}" > {output.output_all}
            """
if config['taxonomy_blast']:
    rule fix_otu_table_blast: 
        input: 
            otu_txt =  os.path.join(config['output_OTU'], "{id}", "otu_taxonomy_{id}_blast.txt"),
            otu_tsv = os.path.join(config['output_OTU'], "{id}", "otu_cluster_{id}.tsv")
        output:
            output_all = os.path.join(config["output_OTU"], "{id}", "otu_table_all_{id}_blast.tsv"),
            output_temp = temp(os.path.join(config["tmp_dir"], "{id}", "otu_taxonomy_{id}_cut_temp_blast.txt")),
            output_temp1 = temp(os.path.join(config["tmp_dir"], "{id}", "otu_taxonomy_{id}_cut_temp1_blast.txt"))
        threads:
            1
        resources:
            mem_mb = 1024,
            runtime = "01:00:00"
        log:
            os.path.join(config["log_dir"], "fix_otu_table_blast", "otu_table_all_{id}.log")
        shell:
            """
                awk -F'\t' '{print $1 "\t" $2}' {input.otu_txt} > {output.output_temp}
                sed -i 's/;size=[0-9]\\+\t/\t/' {output.output_temp}
                awk -i inplace -F'\t' '$2 !~ /^[uU]ncultured/' {output.output_temp}
                (echo -e "OTU ID\tgenus\tspecies\tnotes"; awk -F'\t' 'BEGIN {OFS="\t"} {split($2, words, " "); print $1, "g__"words[1], "s__"words[2], $2}' {output.output_temp}) > {output.output_temp1}
                sed -i '1 s/OTU ID/OTU_ID/' {output.output_temp1}
                sed -i '1 s/#OTU ID/OTU_ID/' {input.otu_tsv}
            """