rule fix_otu_table_97:
    input: 
        otu_txt = os.path.join(['output_OTU'], 'otu_taxonomy_97.txt'),
        otu_tsv = os.path.join(['output_OTU'], 'otu_table_97.tsv')
    output:
        output_all = os.path.join(['output_OTU'], 'otu_table_97_fixed_all.tsv'),
        output_temp = temp(os.path.join(['tmp_dir'], "otutable", 'otu_taxonomy_97_cut_temp.txt'))
    threads:
        1
    log:
        os.path.join(["log_dir"], "fix_otu", "otu_table_97_fixed_all.log")
    shell:
        """
        awk -F "\t" 'OFS="\t"{gsub(/;.*/," ",$1);print $1 $4}' {input.otu_txt} | awk -F ' ' '{print $1"\t"$2}' > {output.output_temp}
        awk 'BEGIN {FS=OFS="\t"} NR==FNR {hold[$1]=$2; next} {print $0, hold[$1]}' {output.output_temp} "{input.otu_tsv}" > {output.output_all}
        """
rule fix_otu_table_99:
    input: 
        otu_txt = os.path.join(['output_OTU'], 'otu_taxonomy_99.txt'),
        otu_tsv = os.path.join(['output_OTU'], 'otu_table_99.tsv')
    output:
        output_all = os.path.join(['output_OTU'], 'otu_table_99_fixed_all.tsv'),
        output_temp = temp(os.path.join(['tmp_dir'], "otutable", 'otu_taxonomy_99_cut_temp.txt'))
    threads:
        1
    log:
        os.path.join(["log_dir"], "fix_otu", "otu_table_99_fixed_all.log")
    shell:
        """
        awk -F "\t" 'OFS="\t"{gsub(/;.*/," ",$1);print $1 $4}' {input.otu_txt} | awk -F ' ' '{print $1"\t"$2}' > {output.output_temp}
        awk 'BEGIN {FS=OFS="\t"} NR==FNR {hold[$1]=$2; next} {print $0, hold[$1]}' {output.output_temp} "{input.otu_tsv}" > {output.output_all}
        """
