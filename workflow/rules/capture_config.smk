import yaml
rule capture_config:
    output:
        config = os.path.join(config["output_dir"], "config.txt")
    threads:
        1
    resources:
        mem_mb = 512,
        runtime = 30,
    run:
        config_combined = config.copy()
        with open(output["config"], "w") as f:
            yaml.dump(config_combined, f)
