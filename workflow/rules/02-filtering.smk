rule filter_fastq:
    input:
        os.path.join(config['tmp_dir'], "samples", "{sample}_concat.fastq")
    output:
        temp(os.path.join(config['tmp_dir'], "samples", "{sample}_filtered.fastq"))
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
        > {output}
        """