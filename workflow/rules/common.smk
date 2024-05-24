## Wildcards helper function - new version, this is the one located in the exp branch.
import glob
import os

# Wildcards helper function
def listFastq(wildcards):
    sample_path = os.path.join(config['input_dir'])
    fastqs = []
    
    if os.path.isdir(sample_path):
        # Check if there are subdirectories
        subdirs = [d for d in os.listdir(sample_path) if os.path.isdir(os.path.join(sample_path, d))]
        
        if subdirs:
            # If there are subdirectories, use them as wildcards.sample
            for subdir in subdirs:
                if wildcards.sample in subdir:
                    fastqs = glob.glob(os.path.join(sample_path, subdir, "*.fastq*"))
        else:
            # If there are no subdirectories, use filenames as wildcards.sample
            files = [f for f in os.listdir(sample_path) if os.path.isfile(os.path.join(sample_path, f))]
            for file in files:
                file_base = os.path.splitext(file)[0]
                if wildcards.sample in file_base:
                    fastqs.append(os.path.join(sample_path, file))
    return fastqs

# merge_abund_tax_blast - rule "10_prep_for_ampvis2.smk - merge_abund_tax_blast"
import pandas as pd
import sys

def merge_abund_tax_blast(tax_file, abund_file, output_file):
    tax = pd.read_csv(tax_file, sep='\t', header=None, names=['OTU', 'Taxonomy'])
    abund = pd.read_csv(abund_file, sep='\t')
    abund.rename(columns={'#OTU ID': 'OTU'}, inplace=True)

    merged_data = pd.merge(abund, tax, on='OTU', how='right')

    # Remove duplicate rows based on the 'OTU' column
    merged_data = merged_data.drop_duplicates(subset='OTU')

    merged_data.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("Usage: python merge_abund_tax_blast.py <tax_file> <abund_file> <output_file>")
    
    tax_file, abund_file, output_file = sys.argv[1:]

    merge_abund_tax_blast(tax_file, abund_file, output_file)

# process_phyloseq_tax.py for rule "11-prep_for_phyloseq_smk - phyloseq_tax_sintax"
import csv
import sys

def process_phyloseq_tax(input_file, output_file):
    with open(input_file, 'r') as csvinput:
        with open(output_file, 'w') as csvoutput:
            reader = csv.reader(csvinput, delimiter='\t')
            writer = csv.writer(csvoutput, delimiter='\t')

            for row in reader:
                writer.writerow([row[0]] + row[-7:])

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    process_phyloseq_tax(input_file, output_file)