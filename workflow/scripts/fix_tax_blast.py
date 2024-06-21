import pandas as pd
import argparse
import re

def extract_genus_species(taxonomy):
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

    return genus, species


def main ():
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]

    data = pd.read_csv(input_file, sep='\t')

    data[['Genus', 'Species']] = data['Taxonomy'].apply(lambda x: pd.Series(extract_genus_species(x)))
    data.drop(columns=['Taxonomy'], inplace=True)
    data.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    main()