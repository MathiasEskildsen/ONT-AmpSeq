configfile: "config/config.yaml"
include: "common.smk"
rule prep_for_ampvis2_sintax:
    input:
        input_all = os.path.join(config["output_dir"], "OTU-tables", "{id}", "otu_table_all_{id}_sintax.tsv")
    output:
        output_all = os.path.join(config['output_dir'], "OTU-tables", "{id}", "otu_table_all_fixed_{id}_sintax.tsv")
    log:
        os.path.join(config['log_dir'], "taxonomy_sintax", 'prep_for_ampvis2_{id}.log')
    threads:
        1
    resources:
        mem_mb = 2048,
        runtime = "01:00:00"
    conda:
        "../envs/generic.yml"
    script:
        "../scripts/prep_for_ampvis2_sintax.py"

rule ampvis2_modifications_sintax:
    input:
        os.path.join(config['output_dir'], "OTU-tables", "{id}", "otu_table_all_fixed_{id}_sintax.tsv")
    output:
        final = os.path.join(config['output_dir'], "final", "{id}", 'OTUtable_tax_{id}_sintax.tsv'), # <- this is the output file ready for ampvis2
        tmp = temp(os.path.join(config['tmp_dir'], "{id}", 'OTUtable_tax_{id}_sintax_temp.tsv'))
    threads:
        1
    resources:
        mem_mb = 1024,
        runtime = "01:00:00"
    conda:
        "../envs/generic.yml"
    log:
        os.path.join(config['log_dir'], "taxonomy_sintax", 'ampvis2_modifications_{id}.log')
    shell:
        """
        modifications=(
            's/kingdom\\s*__\\s*\\([^[:space:]]*\\)/k__\\1/g'
            's/phylum\\s*__\\s*\\([^[:space:]]*\\)/p__\\1/g'
            's/class\\s*__\\s*\\([^[:space:]]*\\)/c__\\1/g'
            's/order\\s*__\\s*\\([^[:space:]]*\\)/o__\\1/g'
            's/family\\s*__\\s*\\([^[:space:]]*\\)/f__\\1/g'
            's/genus\\s*__\\s*\\([^[:space:]]*\\)/g__\\1/g'
            's/species\\s*__\\s*\\([^[:space:]]*\\)/s__\\1/g'
        )
        input="{input}"
        output="{output.tmp}"
        cp $input $output
        for modification in "${{modifications[@]}}"; do
                echo "Applying modifications $modification to $output"
                sed -i -e "$modification" $output
        done
        sed -i '1 s/#OTU ID/OTU/' {output.tmp}
        awk -F'\t' 'BEGIN {{OFS="\t"}} {{for (i=1; i<=NF; i++) gsub(/^[ \t]+|[ \t]+$/, "", $i)}} 1' {output.tmp} > {output.final}
    """
