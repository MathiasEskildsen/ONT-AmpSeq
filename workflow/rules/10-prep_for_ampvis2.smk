configfile: "config/config.yaml"
if config['taxonomy_sintax']:
    rule prep_for_ampvis2_sintax:
        input:
            input_all = os.path.join(config["output_OTU"], "{id}", "otu_table_all_{id}_sintax.tsv")
        output:
            output_all = os.path.join(config['output_taxonomy'], "{id}", "otu_table_all_fixed_{id}_sintax.tsv")    
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
if config['taxonomy_sintax']:
    rule ampvis2_modifications_sintax:
        input:
            os.path.join(config['output_taxonomy'], "{id}", "otu_table_all_fixed_{id}_sintax.tsv")
        output:
            os.path.join(config['output_taxonomy'], "{id}", 'ampvis2_otu_{id}.tsv')
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
                's/kingdom\s*__\s*\\([^[:space:]]*\\)/k__\\1/g'
                's/phylum\s*__\s*\\([^[:space:]]*\\)/p__\\1/g'
                's/class\s*__\s*\\([^[:space:]]*\\)/c__\\1/g'
                's/order\s*__\s*\\([^[:space:]]*\\)/o__\\1/g'
                's/family\s*__\s*\\([^[:space:]]*\\)/f__\\1/g'
                's/genus\s*__\s*\\([^[:space:]]*\\)/g__\\1/g'
                's/species\s*__\s*\\([^[:space:]]*\\)/s__\\1/g'
            )
            input="{input}"
            output="{output}"
            cp $input $output
            for modification in "${{modifications[@]}}"; do
                    echo "Applying modifications $modification to $output"
                    sed -i -e "$modification" $output
            done
            """
configfile: "config/config.yaml"
if config['taxonomy_blast']:
    rule prep_for_ampvis2_blast:
        input:
            otu_table = os.path.join(config['output_OTU'], "{id}", "otu_cluster_{id}.tsv")
            otu_tax = os.path.join(config["tmp_dir"], "{id}", "otu_taxonomy_{id}_cut_temp1_blast.txt")
        output:
            output_all = os.path.join(config['output_taxonomy'], "{id}", "otu_table_all_fixed_{id}_blast.tsv")    
        log:
            os.path.join(config['log_dir'], "taxonomy_blast", 'prep_for_ampvis2_{id}.log')
        threads:
            1
        resources:
            mem_mb = 2048,
            runtime = "01:00:00"
        run:
            import pandas as pd
            otus_tax = pd.read_csv(input.otu_tax, sep='\t')
            otutable = pd.read_csv(input.otu_table, sep='\t')
            merged_data = pd.merge(otutable, otus_tax, on='OTU_ID', how='right')
            merged_data.to_csv(output[0], sep='\t', index=False)