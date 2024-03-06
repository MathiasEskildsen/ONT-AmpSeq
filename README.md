# Snakemake workflow: `smk-ONT_OTU_table`

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.18.2-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/MathiasEskildsen/smk-ONT_OTU_table/workflows/Tests/badge.svg?branch=main)](https://github.com/MathiasEskildsen/smk-ONT_OTU_table/actions?query=branch%3Amain+workflow%3ATests)

## Description
This is a snakemake pipeline, designed to generate OTU-tables from demultiplexed ONT amplicon data. The final outputs are designed to be compatible with R-packages [ampvis2](https://kasperskytte.github.io/ampvis2/index.html) and [phyloSeq](https://github.com/joey711/phyloseq) to visualize the microbial composition of the analyzed samples.
The pipeline expects the input files to be demultiplexed prior to running the pipeline. 
The pipeline filters based on user-input in the config file `config/config.yaml` where it is possible to change filters in regards to amplicon length and quality, using [chopper](https://github.com/wdecoster/chopper). The read characteristics for each sample can be assesed using the shell-script script located at `scripts/nanoplot.sh`. More information regarding usage of the script can be found [here](#usage-of-stats-script).
Biologically meaningful reads from each sample/barcode are clustered into OTU's using [Vsearch](https://github.com/torognes/vsearch) and denoising using [UNOISE3](https://doi.org/10.1093/bioinformatics/btv401) algorithm.
OTU's from every sample/barcode are merged and polished using [Racon](https://github.com/isovic/racon).
Taxonomy is infered to the OTU's by either [Vsearch](https://github.com/torognes/vsearch) using a curated SINTAX database (more information on databases [here](#databases)) or [blastn](https://blast.ncbi.nlm.nih.gov/doc/blast-help/) against a blastn formatted database.

NOTE:
This is a pipeline still being actively developed. 
The workflow was constructed using the following [snakemake template](https://github.com/cmc-aau/snakemake_project_template).

## Requirements
All required tools are automatically installed by Snakemake using conda environments or singularity/apptainer containers, however Snakemake itself needs to be installed first. Load a software module with Snakemake, use a native install, or use the `environment.yml` file to create a conda environment for this particular project using fx `mamba env create -n <snakemake_template> -f environment.yml`.

## Usage through slurm
Adjust the `config.yaml` files under both `config/` and `profiles/` accordingly, then simply run `snakemake --profile profiles/<subfolder>` or submit a SLURM job using the `slurm_submit.sbatch` example script.
The usage of this workflow is also described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<owner>%2F<repo>).

## Usage of stats script
The stats script can be used to visualize read characteristics of the amplicon-sequencing data produced by ONT. The script works on both compressed and decompressed ``` .fastq ``` files. The files can be located in the same directory or individual sub-directores. However, if the files are located in the same directory, files originating from the same barcode have to be merged before running the script. 
Input directory structure examples:
```
../data
└── samples
    ├── barcode01
    │   ├── PAQ88430_pass_barcode01_807aee6b_5f7fc5bf_0.fastq
    │   ├── PAQ88430_pass_barcode01_807aee6b_5f7fc5bf_1.fastq
    │   └── PAQ88430_pass_barcode01_807aee6b_5f7fc5bf_2.fastq
    ├── barcode02
    │   ├── PAQ88430_pass_barcode02_807aee6b_5f7fc5bf_0.fastq
    │   ├── PAQ88430_pass_barcode02_807aee6b_5f7fc5bf_1.fastq
    │   └── PAQ88430_pass_barcode02_807aee6b_5f7fc5bf_2.fastq
    └── barcode03
        ├── PAQ88430_pass_barcode03_807aee6b_5f7fc5bf_0.fastq
        ├── PAQ88430_pass_barcode03_807aee6b_5f7fc5bf_1.fastq
        └── PAQ88430_pass_barcode03_807aee6b_5f7fc5bf_2.fastq
```
```
../data
└── samples
    ├── PAQ88430_pass_barcode01_807aee6b.fastq
    ├── PAQ88430_pass_barcode02_807aee6b.fastq
    └── PAQ88430_pass_barcode03_807aee6b.fastq
```


Usage:
```
-- insert full pipeline name: Nanopore Statistics with NanoPlot
usage: nanoplot [-h] [-o path] [-i path] [-t value] [-j value]

where:
    -h Show this help message.
    -o Path where directories should be created and files should be stored
    -i Full path to .fastq.gz files from Nanopore, example: /Full/Path/to/nanopore_data/ONT_RUN_ID/fastq_pass  
    -j Number of parallel jobs [default = 1]
    -t Number of threads [default = 1]
    Important note:
    Remember to activate your conda environment, containing nanoplot version 1.42.0, before running the script.
    If installed through stats.yml, activate the environment with mamba activate stats.
```
Example command:

```
mamba activate stats
bash ../scripts/nanoplot.sh -o ../out_dir -i ../data/samples -t 1 -j 1 
```
The command will create a directory under your chosen directory (or to a full path) called `out_dir` containing three sub-directories `stats`, `fastqs` and `joblog`. The `stats` sub-directory will contain plots for each sample, in their respective sub-directory, the `LengthvsQualityScatterPlot_dot.png` provides a great overview of the reads. `fastqs` contains unzipped merged fastq files, can be removed. `joblog` contains a text file with the command-line output. Can be useful for debugging.

## Databases
Flesh out description

# TODO
* Create phyloseq-ready output
* Add blast taxonomy option
* Replace `<owner>` and `<repo>` with the correct values in this `README.md` as well as in files under `.github/workflows/`.
* Replace `<snakemake_template>` with the workflow/project name (can be the same as `<repo>`) here as well as in the `environment.yml` and `slurm_submit.sbatch` files.
* Add more requirements to the `environment.yml` file if needed, however tools for each Snakemake rule should **NOT** go here, they should be configured separately for each rule instead in `yaml` files under `envs/`.
* Fill in fields in this `README.md` file, in particular provide a proper description of what the workflow does with any relevant details and configuration.
* The workflow will occur in the public Snakemake workflow catalog once the repository has been made public and the provided GitHub actions finish correctly. Then the link under "Usage" will point to the usage instructions if `<owner>` and `<repo>` were correctly set. If you don't want to publish the workflow just delete the `.github/workflows/` and `.template/` folders and `snakemake-workflow-catalog.yml`.
* Consider the license - [Choose a license](https://choosealicense.com/)
* DELETE this **TODO** section when finished with all of the above, and then start developing your workflow!
