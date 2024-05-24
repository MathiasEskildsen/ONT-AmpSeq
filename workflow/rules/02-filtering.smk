rule filter_fastq:
    input:
        os.path.join(config['tmp_dir'], "samples", "{sample}_concat.fastq")
    output:
        filtered = temp(os.path.join(config['tmp_dir'], "samples", "{sample}_filtered.fastq")),
        total_reads = temp(os.path.join(config['tmp_dir'], "read_count", "{sample}", "{sample}_total_reads_post_filtering.tsv"))
    threads:
        2
    resources:
        mem_mb = 1024,
        runtime = "01:00:00"
    params:
        length_lower_limit = config['length_lower_limit'],
        length_upper_limit = config['length_upper_limit'],
        quality_cut_off = config['quality_cut_off']
    conda:
        "../envs/filtering.yml"
    log:
        os.path.join(config['log_dir'], "filter_fastq", "{sample}.log")
    shell:
        """
        chopper \
        --minlength {params.length_lower_limit} \
        --maxlength {params.length_upper_limit} \
        --headcrop 22 \
        --tailcrop 22 \
        -q {params.quality_cut_off} \
        -t {threads} \
        < {input} \
        > {output.filtered}

        # Count total reads
        num_reads=$(($(wc -l < "{output.filtered}") / 4))
        # Put into a temporary file
        echo -e "Sample\tReads_Post_Filtering\n{wildcards.sample}\t$num_reads" > {output.total_reads}
        """

rule merge_read_count:
    input:
        pre = expand(os.path.join(config['tmp_dir'], "read_count", "{sample}", "{sample}_total_reads_pre_filtering.tsv"), sample=get_samples()),
        post = expand(os.path.join(config['tmp_dir'], "read_count", "{sample}","{sample}_total_reads_post_filtering.tsv"), sample=get_samples())
    output:
        os.path.join(config['output_dir'], "final", "report", "total_reads.tsv")
    threads:
        2
    resources:
        mem_mb = 1024,
        runtime = "01:00:00"
    conda:
        "../envs/filtering.yml"
    log: 
        os.path.join(config['log_dir'], "merge_read_count", "merge_read_count.log")
    shell:
        """
        echo -e "Sample\tReads_Pre_Filtering\tReads_Post_Filtering" > {output}
        for pre_file in {input.pre}; do
            sample=$(basename $pre_file | cut -d'_' -f1)
            reads_pre=$(sed -n '2p' $pre_file | cut -f2)
            post_file=$(echo $pre_file | sed 's/pre_filtering/post_filtering/')
            reads_post=$(sed -n '2p' $post_file | cut -f2)
            echo -e "$sample\t$reads_pre\t$reads_post" >> {output}
        done
        """