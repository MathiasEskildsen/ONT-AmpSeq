import csv
import sys
import logging

# Setup logging
log_file = snakemake.log[0]  # The log file path passed from Snakemake
logging.basicConfig(
    filename=log_file,
    level=logging.DEBUG,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

# Redirect stdout and stderr to the log file
sys.stdout = open(log_file, "a")
sys.stderr = open(log_file, "a")

def process_phyloseq_tax(input_file, tmp_output_file, final_output_file):
    logging.info("Starting the processing of the input file.")
    # First part: Process the initial OTU table
    try:
        with open(input_file, 'r') as csvinput:
            with open(tmp_output_file, 'w') as csvoutput:
                reader = csv.reader(csvinput, delimiter='\t')
                writer = csv.writer(csvoutput, delimiter='\t')

                for row in reader:
                    writer.writerow([row[0]] + row[-7:])
        logging.info("Finished processing the initial OTU table.")
    except Exception as e:
        logging.error(f"An error occurred while processing the initial OTU table: {e}")
        raise
    # Second part: Process the temporary OTU table
    try:
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
        logging.info("Finished processing the temporary OTU table.")
    except Exception as e:
        logging.error(f"An error occurred while processing the temporary OTU table: {e}")
        raise

if __name__ == "__main__":
    # Access the Snakemake variables
    input_file = snakemake.input.input_OTU
    tmp_output_file = snakemake.output.tmp_output
    final_output_file = snakemake.output.phyloseq_tax

    process_phyloseq_tax(input_file, tmp_output_file, final_output_file)
