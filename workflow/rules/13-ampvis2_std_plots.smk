include: "common.smk"
rule ampvis2_std_plots_sintax:
    input:
        OTU_sintax = os.path.join(config['output_dir'], "final", "{id}", "OTUtable_tax_{id}_sintax.tsv")
    output:
        Heatmap = os.path.join(config['output_dir'], "final", "figs", "{id}", "Heatmap_{id}_sintax.png"),
        Rarefraction = os.path.join(config['output_dir'], "final", "figs", "{id}", "Rarefaction_{id}_sintax.png"),
        Ordination = os.path.join(config['output_dir'], "final", "figs", "{id}", "Ordination_{id}_sintax.png"),
        Boxplot = os.path.join(config['output_dir'], "final", "figs", "{id}", "Boxplot_{id}_sintax.png"),
        r_env = os.path.join(config['output_dir'], "final", "R-environment", "{id}", "R-environment_{id}_sintax.RData")
    threads:
        4
    resources:
        mem_mb = 4096,
        runtime = 120,
    conda:
        "../envs/ampvis2.yml"
    params:
        metadata = config['metadata']
    log:
        os.path.join(config["log_dir"], "Ampvis2", "{id}", "ampvis2_plots_sintax.log")
    script:
        "../scripts/Std_plots_ampvis.R"
rule ampvis2_std_plots_blast:
    input:
        OTU_blast = os.path.join(config['output_dir'], "final", "{id}", "OTUtable_tax_{id}_blast.tsv")
    output:
        Heatmap = os.path.join(config['output_dir'], "final", "figs", "{id}", "Heatmap_{id}_blast.png"),
        Rarefraction = os.path.join(config['output_dir'], "final", "figs", "{id}", "Rarefaction_{id}_blast.png"),
        Ordination = os.path.join(config['output_dir'], "final", "figs", "{id}", "Ordination_{id}_blast.png"),
        Boxplot = os.path.join(config['output_dir'], "final", "figs", "{id}", "Boxplot_{id}_blast.png"),
        r_env = os.path.join(config['output_dir'], "final", "R-environment", "{id}", "R-environment_{id}_blast.RData")
    threads:
        4
    resources:
        mem_mb = 4096,
        runtime = 120,
    conda:
        "../envs/ampvis2.yml"
    params:
        metadata = config['metadata']
    log:
        os.path.join(config["log_dir"], "Ampvis2", "{id}","ampvis2_plots_blast.log")
    script:
        "../scripts/Std_plots_ampvis.R"
