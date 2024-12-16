## Configuration
- `input_dir` : Path to the input folder, containing fastq files in compressed or decompressed format. The pipeline expects the input files to conform to 1 of 2 directory structures.
1: `input_dir` contains subfolders for each sampleID/barcode, if that is the case all fastq files in each subfolder are concatenated and the subfolder name is used as a sample ID downstream. This is usually the "fastq_pass" folder from nanopore sequencing and basecalling output (atleast when using Guppy).
2: `input_dir` contains already concatenated fastq files, directly located in `input_dir`. If that is the case, the pipeline uses the entire filename as a sample ID downstream. This is usually the case output from Dorado re-basecalling with demultiplexing enabled.
- `output_dir`: Output directory with the final results and a few intermediary files, that can be used for other downstream purposes if desired.
- `tmp_dir`: Directory for temporary files.
- `log_dir`: Directory for log files for all invoked rules.
- `db_path_sintax`: Database to infer taxonomy using the SINTAX algorithm. Contains sequenceID, taxonomy string and fasta sequence. 
- `db_path_blast`: Nucleotide blast formatted database to infer taxonomy using BLASTn algorithm.
- `evalue`: E-value cutoff for blast. Default = 1e-10.
- `length_lower_limit`: Argument passed on to `chopper` for filtering reads. Appropriate values depends on amplicon length. This can be checked by running the helper script scripts/nanoplot.sh 
- `length_upper_limit`: Argument passed on to `chopper` for filtering reads. Appropriate values depends on amplicon length. This can be checked by running the helper script scripts/nanoplot.sh 
- `quality_cut_off`: Argument passed on to `chopper` for filtering reads. Appropriate value depends on the quality of your sequencing data. This can be checked by running the helper script scripts/nanoplot.sh. It is recommended to pick a Q-score >20, if your data permits it.
- `max_threads`: Maximum number of threads that can be used for any rule.
- `include_blast_output`: Default = True. If true snakemake will output a final OTU-table with taxonomy infered from a blastn search against a nt blast database.
- `include_sintax_output`: Default = True. If true snakemake will output a final OTU-table with taxonomy infered from a sintax formatted database.
- `ids`: Clustering identity for OTUs. Default is 97% and 99%. Use "." decimal seperator i.e 99.9.
- `primer_f`: Forward primer length to trim off. Default = 22
- `primer_r`: Reverse primer length to trim off. Defeault = 22
- `f`: Minimap2 mapping option. Default = 0.0002 More information at [here](https://lh3.github.io/minimap2/minimap2.html#10).
- `K`: Minimap2 mapping option. Default = 500M More information at [here](https://lh3.github.io/minimap2/minimap2.html#10).