configfile: "config/config.yaml"

rule filter_fastq:
    input: os.path.join(config['tmp_dir'], "samples", "{sample}_concat.fastq.gz")
    output: temp(os.path.join(config['tmp_dir'], "samples", "{sample}_filtered.fastq"))
    threads: config['max_threads']
    params:
        length_lower_limit = config['length_lower_limit'],
        length_upper_limit = config['length_upper_limit'],
        quality_cut_off = config['quality_cut_off']
    conda:
        "../envs/filtering.yml"
    log:
        os.path.join(config['log_dir'], "filtering", "{sample}_filtering.log")
    shell:
        """
    chopper \
        {input} \
        --minlength {params.length_lower_limit} \
        --maxlength {params.length_upper_limit} \
        --headcrop 22 \
        --tailcrop 22 \
        -q {params.quality_cut_off} \
        -t {threads} \
        > {output}
    """
