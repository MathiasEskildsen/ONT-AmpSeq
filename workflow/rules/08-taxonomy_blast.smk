rule rule_name:
    input:
        input_file = "path/to/input/file",
        # add more input files if needed
    output:
        output_file = "path/to/output/file",
        # add more output files if needed
    params:
        param1 = "value1",
        # add more parameters if needed
    threads: 1
    resources:
        mem_mb = 1000,
        # add more resources if needed
    log:
        "path/to/log/file"
    shell:
        """
        # add your shell command here
        """
