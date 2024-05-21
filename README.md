# Snakemake workflow: `ONT-AmpSeq`

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.18.2-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/MathiasEskildsen/ONT-AmpSeq/workflows/Tests/badge.svg?branch=main)](https://github.com/MathiasEskildsen/ONT-AmpSeq/actions?query=branch%3Amain+workflow%3ATests)

## Description
This is a snakemake workflow, designed to generate OTU-tables from demultiplexed ONT amplicon data. The final outputs are designed to be compatible with R-packages [ampvis2](https://kasperskytte.github.io/ampvis2/index.html) and [phyloSeq](https://github.com/joey711/phyloseq) to visualize the microbial composition of the analyzed samples.
The workflow expects the input files to be demultiplexed and basecalled prior to running the workflow. 
The workflow filters based on user-input in the config file `config/config.yaml` where it is possible to change filters in regards to amplicon length and quality, using [chopper](https://github.com/wdecoster/chopper). The read characteristics for each sample can be assesed using the shell-script script located at `scripts/nanoplot.sh`. More information regarding usage of the script can be found [here](#usage-of-stats-script).
Biologically meaningful reads from each sample/barcode are clustered into OTU's using [Vsearch](https://github.com/torognes/vsearch) and denoising using [UNOISE3](https://doi.org/10.1093/bioinformatics/btv401) algorithm.
OTU's from every sample/barcode are merged and polished using [Racon](https://github.com/isovic/racon).
Taxonomy is infered to the OTU's by either [Vsearch](https://github.com/torognes/vsearch) using a curated SINTAX database (more information on databases [here](#databases)) or [blastn](https://blast.ncbi.nlm.nih.gov/doc/blast-help/) against a blastn formatted database.

NOTE:
This is a workflow still being actively developed.
The workflow was constructed using the following [snakemake template](https://github.com/cmc-aau/snakemake_project_template).

## Table of Contents
- [Requirements](#requirements)
- [Usage of Workflow with Snakedeploy](#usage-of-workflow-snakedeploy)
- [Uage of Workflow AAU Biocloud Users](#usage-of-workflow-aau-biocloud-hpc-users)
- [Outputs](#outputs)
- [Stats script](#usage-of-stats-script)
- [Database Choice](#databases)
## Requirements
Requires a Linux OS or WSL 

This workflow requires conda or mamba to install the required tools for the pipeline. Snakemake will automatically install the correct version of the tools required for the pipeline. However, for the first time use you need to install conda or mamba and create an environment containing snakemake and snakedeploy.

Conda can be installed by following this [guide](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and Mamba can be installed by following this [guide]().
It is recommended to follow the original documentation, however below is the commands used to freshly install the software on a Linux machine as per their documentation (14-05-2024).
Miniconda:
```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
```  
Initialize Miniconda:
```
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```
Add channels (for fresh install):
```
conda config --add channels conda-forge
conda config --add channels bioconda
```

Now your conda install is set-up and you're ready to set-up your environment.

## Usage of workflow (Snakedeploy)
The usage of this workflow is also described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=MathiasEskildsen/ONT-AmpSeq).
## Usage with snakedeploy
### Step 1: Install Snakemake and Snakedeploy
Snakemake and Snakedeploy are best installed via the [Mamba package manager](https://github.com/mamba-org/mamba) (a drop-in replacement for conda). If you have neither Conda nor Mamba, it can be installed via [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge). For other options see [here](https://github.com/mamba-org/mamba). Or refer to the Miniconda installation and quick set-up in [Requirements](#requirements).

Given that Mamba (If Miniconda is installed, the mamba can be changed for conda) is installed, set-up your enviroment by running:
`mamba create -c conda-forge -c bioconda --name snakemake snakemake=7.18.2 snakedeploy`
This installs both Snakemake and Snakedeploy in their own isolated environment. For all the following commands ensure that this environment is activated with the following command:
`mamba activate snakemake`
### Step 2: Deploy workflow
Given that Snakemake and Snakedeploy are installed and activated (see step 1), the workflow can be deployed as follows:
First, create an appropriate project working directory on your system and choose it:
```
mkdir -p path/to/project-workdir
cd path/to/project-workdir
```
In the following steps, we will assume that you have chosen and are inside of that directory.
Second, run:
```
snakedeploy deploy-workflow https://github.com/MathiasEskildsen/ONT-AmpSeq . --branch main
```

Snakedeploy will create two folders `workflow` and `config`. The former contains the deployment of the chosen workflow as a [Snakemake module](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-and-combining-pre-exising-workflows), the latter contains configuration files which will be modified in the next step in order to configure the workflow to your needs. Later, when executing the workflow, Snakemake will automatically find the main `Snakefile` in the `workflow` subfolder. 

Third, consider to put this directory under version control, e.g. by managing it [via a (private) Github repository](https://docs.github.com/en/migrations/importing-source-code/using-the-command-line-to-import-source-code/adding-locally-hosted-code-to-github)

### Step 3: Configure the workflow
The workflow needs to be configured according to your directory structure and according to your needs. Below is an explaination of the settings that can be configured. These can be changed directly in the config file located at `path/to/project-workdir/config/config.yaml` or changed by command line arguments, explained later.
* `input_dir`: Path to the input folder, containing fastq files in compressed or decompressed format. The pipeline expects the input files to conform to 1 of 2 directory structures, see [here](#usage-of-stats-script) for more information of directory structures.
* `output_dir`: Path to output directory with the final results of the pipeline and a few intermediary files, thay might prove useful for other purposes.
* `tmp_dir`: Directory for temporary files, temporary files will be removed after a succesful run.
* `log_dir`: Directory for log files for all invoked rules.
* `db_path_sintax`: Database to infer taxonomy using the SINTAX algorithm. Should contain SequenceID, taxonomy string and a fasta sequence. More information at [databases](#databases).
* `db_path_blast`: Nucleotide blast formatted database to infer taxonomy using BLASTn algorithm. More information at [databases](#databases).
* `evalue`: E-value cutoff for blast. Default = 1e-10.
* `length_lower_limit`: Default = 1200. Argument passed on to `chopper` for filtering reads. Appropriate values depends on amplicon length. This can be checked by running the helper [stats script](#usage-of-stats-script) scripts/nanoplot.sh.
* `length_upper_limit`: Default = 1600. Argument passed on to `chopper` for filtering reads. Appropriate values depends on amplicon length. This can be checked by running the helper [stats script](#usage-of-stats-script) scripts/nanoplot.sh.
* `quality_cut_off`: Default = 23. Argument passed on to `chopper` for filtering reads. Appropriate value depends on the quality of your sequencing data. This can be checked by running the helper [stats script](#usage-of-stats-script) scripts/nanoplot.sh. It is recommended to pick a Q-score >20, if your data permits it. 
* `max_threads`: Maximum number of threads that can be used for any given rule.
* `include_blast_output`: Default = true. If true snakemake will output a final OTU-table with taxonomy infered from a blastn search against a nt blast database.
* `include_sintax_output`: Default = true. If true snakemake will output a final OTU-table with taxonomy infered from a sintax formatted database.

As previously mentioned, the workflow configurations can also be changed directly in the command line. Every configuration can be changed, keep in mind, that configurations without specified changes in the command line will use the value specified in the configuration file (`path/to/project-workdir/config/config.yaml`).
Example of changing a few configurations from their default:
```
cd path/to/project-workdir
mamba activate snakemake
snakemake --cores all --use-conda --config include_blast_output=False db_path_sintax=/path/to/SINTAX_DATABASE.fa length_lower_limit=400 length_upper_limit=800 quality_cut_off=20
```
The code snippet above will; choose your project working directory, enabling snakemake to locate the snakefile and configuration file. Activate your environment containing snakemake, as described in [step1](#step-1-install-snakemake-and-snakedeploy). Finally, it will run the snakemake workflow, filtering out reads shorter than 400bp, longer than 800bp and with a Q-score <20. It will output OTU-tables with taxonomy annotated by a sintax database. 
### Step 4: Run the workflow
Given that the workflow has been properly deployed and configured, it can be executed as follows.
For running the workflow while deploying any necessary software via conda using the [Mamba package manager](https://github.com/mamba-org/mamba), run Snakemake with:
```
snakemake --cores all --use-conda
```
Given that you have chosen your project working-directory as previously stated. Snakemake will automatically detect the main `Snakefile` in the `workflow` subfolder and execute the workflow module that has been defined by the deployment in step 2.

For further options, fx. for cluster and cloud execution, see [the docs](https://snakemake.readthedocs.io/en/stable/). If you are an AAU user, see [this](#usage-of-workflow-through-slurm-aau-biocloud-users) section.
### Step 5: Generate report
After finalizing your data analysis, you can automatically generate an interactive visual HTML report for inspection of results together with parameters and code inside of the browser using:
```
snakemake --report report.zip
```
The resulting `report.zip` file can be passed on to collaborators, provided as a supplementary file in publications, or uploaded to a service like [Zenodo](https://zenodo.org/) in order to obtain a citable [DOI](https://en.wikipedia.org/wiki/Digital_object_identifier).

## Usage of workflow (AAU BioCloud HPC users)
AAU BioCloud HPC users can also use the snakedeploy [step-by-step](#usage-with-snakedeploy), however it is recommended to follow the guide below as this will include scripts to help you submit jobs via. SLURM. If you want further guidance on snakemake and BioCloud usage refer to the [user guide](https://cmc-aau.github.io/biocloud-docs/guides/snakemake/intro/).

Change the path to your project-directory
```
cd /path/to/project-dir
wget -O https://github.com/MathiasEskildsen/ONT-AmpSeq/archive/refs/heads/main.tar.gz | tar -xz
```
Install dependencies (BioCloud users already have mamba installed natively)
```
cd ONT-AmpSeq-main
mamba env create -f environment.yml
```
Configure `config/config.yaml` as described [previously](#step-3-configure-the-workflow). Then simply run `snakemake --profile profiles/biocloud` or submit a SLURM job using the `slurm_submit.sbatch` example script. If running the pipeline through `slurm_submit.sbatch`, remember to change the `#SBATCH` arguments at the top of the script, to fit your run (--job-name, --mail and --time). Snakemake will automatically queue jobs with the necesarry ressources so you do not need to change ressources specified in `slurm_submit.sbatch`.

## Outputs
NOTE: `{id}` refers to the percentage identity which reads should be similar in order to be clustered into an OTU. Defaults = [97%, 99%]. Outputs can be annotated with either a blastn or SINTAX formatted database or both.
* `OTUtable_tax_{id}_sintax.tsv`: Matrix containing number of reads per sample per OTU. Taxonomy of each OTU is annotated by a SINTAX formatted [database](#databases). The OTU table is formatted to be ready for data analysis using [Ampvis2](https://kasperskytte.github.io/ampvis2/index.html).
* `phyloseq_tax_{id}_sintax.tsv`: Matrix containing taxonomy for each OTU annotated by a SINTAX formatted [database](#databases).
* `phyloseq_abundance_{id}_sintax.tsv`: Matrix containing OTU abundance information for each sample.
* `OTUtable_tax_{id}_blast.tsv`: Matrix containing number of reads per sample per OTU. Taxonomy of each OTU is annotated by a blastn formatted [database](#databases). The OTU table is formatted to be ready for data analysis using [Ampvis2](https://kasperskytte.github.io/ampvis2/index.html). NOTE: Blastn results can give a lot of edge-cases in relation to the formatting of output taxonomy. Be sure to double check annotated taxonomy when using this approach. 
* `phyloseq_tax_{id}_blast.tsv`: Matrix containing taxonomy for each OTU annotated by a blast formatted [database](#databases).
* `phyloseq_abundance_{id}_blast.tsv`: Matrix containing OTU abundance information for each sample.
* `total_reads.tsv`: Number of reads in each sample pre- and post-filtering.
## Usage of stats script
The stats script can be used to visualize read characteristics of the amplicon-sequencing data produced by ONT. The script works on both compressed and decompressed `.fastq` files. The files can be located in the same directory or individual sub-directores. However, if the files are located in the same directory, files originating from the same barcode have to be merged before running the script. 
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
MiDAS SINTAX database can be downloaded [here](https://www.midasfieldguide.org/guide/downloads)

For more information regarding blastn databases look [here](https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.makeblastdb_application_opt/)

# TODO
* Add release version
* Add description of final outputs
* Flesh out readme with database choice
* Add split_taxonomy.py script to better handle edge case annotation from blast
* Add .test to pass linting test