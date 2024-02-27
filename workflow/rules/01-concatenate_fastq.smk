import glob
import os
# Helper function to list all fastq files per wildcard (subfolder/sample)
def listFastq(wildcards):
    sample_path = os.path.join(config['input_dir'])
    if os.path.isdir(sample_path) and any(os.path.isdir(os.path.join(sample_path, i)) for i in os.listdir(sample_path)):
        fastqs = glob.glob(os.path.join(config['input_dir'], wildcards.sample, "*.fastq*"))
    else:
        for file_path in sample_dirs:
            filename_without_extension = os.path.splitext(os.path.basename(file_path))[0]
            filename_without_extension = wildcards.sample
            fastqs = glob.glob(os.path.join(config['input_dir'], wildcards.sample))
    return fastqs
rule concatenate_fastq:
    input:
        listFastq
    output:
        temp(os.path.join(config['tmp_dir'], "samples", "{sample}_concat.fastq"))
    resources:
        mem_mb = 512,
        runtime = "01:00:00"
    threads: 1
    log:
        os.path.join(config["log_dir"], "concatenate_fastq", "{sample}.log")
    shell:
        """
        # Check if input files are compressed
        for file in {input}; do
            if gzip -t $file 2>/dev/null; then
                echo "Decompressing $file"
                zcat $file >> {output}
            else
                echo "Copying $file"
                cat $file >> {output}
            fi
        done
        """