import glob
import os
config = {'input_dir': "data/samples"}
class Wildcards:
    def __init__(self, sample):
        self.sample = sample

wildcards = Wildcards("fSRR23176512.fastq.gz")  # replace "sample1" with the name of one of your samples
def listFastq(wildcards):   
    fastqs = glob.glob(os.path.join(config['input_dir'], wildcards.sample, "*.fastq*"))
    print(fastqs)  # Debug statement
    return fastqs
print(listFastq(wildcards))