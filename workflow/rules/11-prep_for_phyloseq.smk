include: "common.smk"
rule phyloseq_abund_sintax:
    input:
        os.path.join(config['output_dir'], "cluster", "{id}", "otu_cluster_{id}.tsv")
    output:
        phyloseq_abundance = os.path.join(config['output_dir'], "final", "{id}", 'phyloseq_abundance_{id}_sintax.tsv')
    log: 
        os.path.join(config['log_dir'], "phyloseq_abundance", 'phyloseq_{id}.log')
    threads:
        1
    conda:
        "../envs/generic.yml"
    resources:
        mem_mb = 1024,
        runtime = "01:00:00"
    shell:
        """
        cp {input} {output}
        sed -i '1 s/#OTU ID/OTU/' {output}
        """

rule phyloseq_tax_sintax:
    input:
        input_OTU = os.path.join(config['output_dir'], "final", "{id}", 'OTUtable_tax_{id}_sintax.tsv')
    output:
        phyloseq_tax = os.path.join(config['tmp_dir'], "{id}", 'phyloseq_tax_tmp_{id}_sintax.tsv')
    log: 
        os.path.join(config['log_dir'], "phyloseq_tax", 'phyloseq_{id}.log')
    threads:
        1
    resources:
        mem_mb = 1024,
        runtime = "01:00:00"
    run:
        import csv
        import sys

        input_file = input['input_OTU']
        output_file = output['phyloseq_tax']
        with open(input_file, 'r') as csvinput:
            with open(output_file, 'w') as csvoutput:
                reader = csv.reader(csvinput, delimiter='\t')
                writer = csv.writer(csvoutput, delimiter='\t')

                for row in reader:
                    writer.writerow([row[0]] + row[-7:])
rule phyloseq_tax1_sintax:
    input:
        os.path.join(config['tmp_dir'], "{id}", 'phyloseq_tax_tmp_{id}_sintax.tsv')
    output:
        os.path.join(config['output_dir'], "final", "{id}", 'phyloseq_tax_{id}_sintax.tsv')
    log: 
        os.path.join(config['log_dir'], "phyloseq_tax", 'phyloseq_{id}.log')
    threads:
        1
    resources:
        mem_mb = 1024,
        runtime = "01:00:00"
    run:
        with open(input[0], 'r') as f_in, open(output[0], 'w') as f_out:
            header = next(f_in)
            f_out.write(header)
            for line in f_in:
                fields = line.split('\t')
                modified_fields = [fields[0]] + [field[3:] for field in fields[1:]]
                f_out.write('\t'.join(modified_fields) + '\n')
