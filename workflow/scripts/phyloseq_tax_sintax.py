import csv

def process_phyloseq_tax(input_file, tmp_output_file, final_output_file):
    # First part: Process the initial OTU table
    with open(input_file, 'r') as csvinput:
        with open(tmp_output_file, 'w') as csvoutput:
            reader = csv.reader(csvinput, delimiter='\t')
            writer = csv.writer(csvoutput, delimiter='\t')

            for row in reader:
                writer.writerow([row[0]] + row[-7:])

    # Second part: Process the temporary OTU table
    with open(tmp_output_file, 'r') as csvinput:
        with open(final_output_file, 'w') as csvoutput:
            reader = csv.reader(csvinput, delimiter='\t')
            writer = csv.writer(csvoutput, delimiter='\t')

            # Write header
            header = next(reader)
            writer.writerow(header)

            # Process and write each row
            for row in reader:
                modified_fields = [row[0]] + [field[3:] for field in row[1:]]
                writer.writerow(modified_fields)

if __name__ == "__main__":
    # Access the Snakemake variables
    input_file = snakemake.input.input_OTU
    tmp_output_file = snakemake.output.tmp_output
    final_output_file = snakemake.output.phyloseq_tax

    process_phyloseq_tax(input_file, tmp_output_file, final_output_file)
