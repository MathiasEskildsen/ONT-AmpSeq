import glob
import os

# helper function to list all fastq files per wildcard (subfolder/sample) 
def listFastq(wildcards):
  fastqs = glob.glob(os.path.join(config['input_dir'], wildcards.sample, "*.fastq.gz"))
  return fastqs

rule concatenate_fastq:
  input:
    listFastq
  output:
    temp(os.path.join(config['tmp_dir'], "samples", "{sample}_concat.fastq.gz"))
  resources:
    mem_mb = 512
  threads: 1
  log:
    os.path.join(config["log_dir"], "concatenate_fastq", "{sample}.log")
  shell:
    "cat {input} > {output}"
