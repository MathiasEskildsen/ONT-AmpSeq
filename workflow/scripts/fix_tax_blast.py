import pandas as pd
import re
import logging
import sys

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

def extract_genus_species(taxonomy):
    logging.debug(f"Extracting genus and species from taxonomy: {taxonomy}")
    genus_match = re.search(r'g:([^,]+)', taxonomy)
    species_match = re.search(r's:([^,;]+)', taxonomy)
    
    if genus_match and species_match:
        genus = f"g__{genus_match.group(1)}"
        species = f"s__{species_match.group(1)}"
    else:
        words = taxonomy.split()
        genus = ""
        species = ""
        
        for i in range(len(words)):
            if words[i] in ["MAG:", "Uncultured", "cf.", "isolate"]:
                continue
            elif (words[i] in ["sp.", "clone"]) and i < len(words) - 1:
                genus = f"g__{words[i-1]}"
                species = f"s__{words[i+1]}"
                break
            else:
                genus = f"g__{words[i]}"
                if i < len(words) - 1 and words[i+1] not in ["sp.", "clone"]:
                    species = f"s__{words[i+1]}"
                break

    logging.debug(f"Extracted genus: {genus}, species: {species}")
    return genus, species

def main():
    logging.info("Starting fix_tax_blast script")
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]

    try:
        logging.info(f"Reading input file: {input_file}")
        data = pd.read_csv(input_file, sep='\t')

        logging.info("Extracting genus and species for each row")
        data[['Genus', 'Species']] = data['Taxonomy'].apply(lambda x: pd.Series(extract_genus_species(x)))
        data.drop(columns=['Taxonomy'], inplace=True)

        logging.info(f"Writing output file: {output_file}")
        data.to_csv(output_file, sep='\t', index=False)

        logging.info("Finished fix_tax_blast script")
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()