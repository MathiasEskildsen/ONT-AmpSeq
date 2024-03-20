configfile: "config/config.yaml"

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
    run:
        import csv
        input_file = input["input_all"]
        output_file = output["output_all"]
        new_column_headers = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

        # Create a dictionary to map field names to new column values
        field_mapping = {
            "d": "kingdom",
            "k": "kingdom",
            "p": "phylum",
            "c": "class",
            "o": "order",
            "f": "family",
            "g": "genus",
            "s": "species",
            }

        # Function to extract and map the values to the new columns
        def extract_and_map_values(input_row):
            last_column = input_row.pop()  # Remove the last column
            field_values = {field: "" for field in new_column_headers}
            parts = last_column.split(",")
            
            for part in parts:
                field_prefix, value = part.split(":") if ":" in part else (None, part)
                new_column = field_mapping.get(field_prefix)
                if new_column:
                    field_values[new_column] = f"{new_column}__{value}"

            new_row = input_row + [field_values[field] for field in new_column_headers]

            return new_row

        # Open the input and output files
        with open(input_file, mode="r") as infile, open(output_file, mode="w", newline="") as outfile:
            reader = csv.reader(infile, delimiter="\t")
            writer = csv.writer(outfile, delimiter="\t")

            # Write the headers with the new columns
            headers = next(reader)[:-1] + new_column_headers
            writer.writerow(headers)

            # Process and write the data
            for row in reader:
                modified_row = extract_and_map_values(row)
                writer.writerow(modified_row)

rule ampvis2_modifications_sintax:
    input:
        os.path.join(config['output_dir'], "OTU-tables", "{id}", "otu_table_all_fixed_{id}_sintax.tsv")
    output:
        final = os.path.join(config['output_dir'], "final", "{id}", 'OTUtable_tax_{id}_sintax.tsv'), # <- this is the output file ready for ampvis2
        tmp = os.path.join(config['tmp_dir'], "{id}", 'OTUtable_tax_{id}_sintax_temp.tsv')
    threads:
        1
    resources:
        mem_mb = 1024,
        runtime = "01:00:00"
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
        import pandas as pd

        tax = pd.read_csv(input["tax"], sep='\t', header=None, names=['OTU', 'Taxonomy'])
        abund = pd.read_csv(input["abund"], sep='\t')
        abund.rename(columns={'#OTU ID': 'OTU'}, inplace=True)

        merged_data = pd.merge(abund, tax, on='OTU', how='right')
        merged_data.to_csv(output[0], sep='\t', index=False)