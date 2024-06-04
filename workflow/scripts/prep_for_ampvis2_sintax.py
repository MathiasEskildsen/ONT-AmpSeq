import csv
import os

# Extract and map values from the input row
def extract_and_map_values(input_row, new_column_headers, field_mapping):
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

def main(input_file, output_file):
    new_column_headers = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
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

    with open(input_file, mode="r") as infile, open(output_file, mode="w", newline="") as outfile:
        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t")

        # Write the headers with the new columns
        headers = next(reader)[:-1] + new_column_headers
        writer.writerow(headers)

        # Process and write the data
        for row in reader:
            modified_row = extract_and_map_values(row, new_column_headers, field_mapping)
            writer.writerow(modified_row)

if __name__ == "__main__":
    input_file = snakemake.input.input_all
    output_file = snakemake.output.output_all
    main(input_file, output_file)
