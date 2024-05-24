include: "common.smk"

rule concatenate_fastq:
    input:
        lambda wildcards: listFastq(wildcards)
    output:
        concat = temp(os.path.join(config['tmp_dir'], "samples", "{sample}_concat.fastq")),
        total_reads = temp(os.path.join(config['tmp_dir'], "read_count", "{sample}", "{sample}_total_reads_pre_filtering.tsv"))
    resources:
        mem_mb = 512,
        runtime = "01:00:00"
    conda:
        "../envs/generic.yml"
    threads: 1
    log:
        os.path.join(config["log_dir"], "concatenate_fastq", "{sample}.log")
    shell:
        """
        # Check if input files are compressed
        for file in {input}; do
            if gzip -t $file 2>/dev/null; then
                echo "Decompressing $file"
                zcat $file >> {output.concat}
            else
                echo "Copying $file"
                cat $file >> {output.concat}
            fi
        done
        # Count total reads
        num_reads=$(($(wc -l < "{output.concat}") / 4))
        # Put into a temporary file
        echo -e "Sample\tReads_Pre_Filtering\n{wildcards.sample}\t$num_reads" > {output.total_reads}
        """