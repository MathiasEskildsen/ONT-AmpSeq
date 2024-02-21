import glob
import os
config = {'input_dir': "data/samples"}
wildcards = type('', (), {})()  # Create an empty object
wildcards.sample = '*.fastq'  # Replace with your sample name
def listFastq(wildcards):
    # If it's a file, remove the extension
    sample_basename = os.path.splitext(wildcards.sample)[0]
    search_path = os.path.join(config['input_dir'], f"{sample_basename}*.fastq*")
    print(f"Searching in: {search_path}")
    fastqs = glob.glob(search_path)
    return fastqs

print(listFastq(wildcards))