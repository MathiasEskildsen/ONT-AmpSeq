rule map2db:
  input:
    os.path.join(config['tmp_dir'], "samples", "{sample}_concat.fastq.gz")
  output:
    os.path.join(config['output_dir'], "{sample}.sam")
  resources:
    #depending on the tool memory usage usually scales with threads 
    #and/or input/database file size. Can calculate dynamically
    mem_mb = 10240
  threads: config['max_threads']
  params:
    db_path = config['db_path']
  conda:
    "../envs/map2db.yml"
  log:
    os.path.join(config["log_dir"], "map2db", "{sample}.log")
  shell:
    """
    minimap2 \
      -ax map-ont \
      -K20M \
      -t {threads} \
      --secondary=no \
      {params.db_path} \
      {input} \
      | samtools \
        view \
        -F 4 \
        -F 256 \
        -F 2048 \
        --threads {threads} \
        -o {output}
    """
