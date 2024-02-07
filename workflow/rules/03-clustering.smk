configfile: "config/config.yaml"

import os
# Function to change fastq headers to fasta headers
# And create a list of all new files with fasta headers
def fastq_to_fasta(input_fastq, output_fasta):
    with open(input_fastq, "r") as f_in, open(output_fasta, "w") as f_out:
        for line in f_in:
            if line.startswith("@"):
                f_out.write(">" + line[1:])
            else:
                f_out.write(line)

# Get a list of all fastq files in the tmp_dir
fastq_files = [f for f in os.listdir(config["tmp_dir"]) if f.endswith("_filtered.fastq")]

rule convert_to_fasta:
    input:
        os.path.join(config["tmp_dir"], "samples", "{sample}_filtered.fastq")
    output:
        os.path.join(config["tmp_dir"], "samples", "{sample}_filtered.fasta")
    threads: 1
    resources:
        mem_mb = 1024
    run:
        fastq_to_fasta({input}, {output})

rule vsearch_cluster:
    input:
        os.path.join(config["tmp_dir"], "samples", "{sample}_filtered.fasta")
    output:
        os.path.join(config["output_cluster"], "{sample}.cluster.fasta")
    threads:
        1
    resources:
        mem_mb = 4096
    conda:
        "../envs/vsearch.yml"
    log:
        os.path.join(config["log_dir"], "{sample}.cluster.log")
    shell:
        """"
    vsearch \
        --cluster_unoise {input} \
        --minsize \
        --threads {threads} \
        --centroids {output}
    """


