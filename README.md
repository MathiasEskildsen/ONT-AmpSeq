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

## Requirements
All required tools are automatically installed by Snakemake using conda environments, however Snakemake itself needs to be installed first. Load a software module with Snakemake, use a native install, or use the `environment.yml` file to create a conda environment for this particular project using fx `mamba env create -n <snakemake_template> -f environment.yml`.

## Usage of workflow
The usage of this workflow is also described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=MathiasEskildsen/ONT-AmpSeq).

## Usage
### Step 1: Install Snakemake and Snakedeploy
Snakemake and Snakedeploy are best installed via the [Mamba package manager](https://github.com/mamba-org/mamba) (a drop-in replacement for conda). If you have neither Conda nor Mamba, it can be installed via [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge). For other options see [here](https://github.com/mamba-org/mamba).
Given that Mamba is installed, run:
`mamba create -c conda-forge -c bioconda --name snakemake snakemake snakedeploy`
This install both Snakemake and Snakedeploy in their own isolated environment. For all the following commands ensure that this environment is activated with the following command:
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

Snakedeploy will create two folders `workflow` and `config`. The former contains the deployment of the chosen workflow as a [Snakemake module](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-and-combining-pre-exising-workflows), the latter contains configuration files which will be modified in the next step in order to configure the workflow to your needs. Later, when executing the workflow, SNakemake will automatically find the main `Snakefile` in the `workflow` subfolder. 

Third, consider to put this directory under version contraol, e.g. by managing it [via a (private) Github repository](https://docs.github.com/en/migrations/importing-source-code/using-the-command-line-to-import-source-code/adding-locally-hosted-code-to-github)
### Step 3: Configure the workflow
* `input_dir`: Path to the input folder, containing fastq files in compressed or decompressed format. The pipeline expects the input files to conform to 1 of 2 directory structures, see [here](## Usage of stats script) for more information of directory structures.
* `output_dir`: Path to output directory with the final results of the pipeline and a few intermediary files, thay might prove useful for other purposes.
* `tmp_dir`: Directory for temporary files.
* `log_dir`: Directory for log files for all invoked rules.
* `db_path_sintax`: Database to infer taxonomy using the SINTAX algorithm. Should contain SequenceID, taxonomy string and a fasta sequence.
* `db_path_blast`: Nucleotide blast formatted database to infer taxonomy using BLASTn algorithm.
* `evalue`: E-value cutoff for blast. Default = 1e-10.
* `length_lower_limit`: Argument passed on to `chopper` for filtering reads.
* `length_upper_limit`: Argument passed on to `chopper` for filtering reads.
* `quality_cut_off`: Argument passed on to `chopper` for filtering reads.
* `max_threads`: Maximum number of threads that can be used for any given rule

### Step 4: Run the workflow
Given that the workflow has been properly deployed and configured, it can be executed as follows.
For running the workflow while deploying any necessary software via conda using the [Mamba package manager](https://github.com/mamba-org/mamba), run Snakemake with:
```
snakemake --cores all --use-conda
```
Given that you have chosen your project working-directory as previously stated. Snakemake will automatically detect the main `Snakefile` in the `workflow` subfolder and execute the workflow module that has been defined by the deployment in step 2.

For further options, fx. for cluster and cloud execution, see [the docs](https://snakemake.readthedocs.io/en/stable/). If you are an AAU user, see [this](## Usage of workflow through SLURM (AAU Biocloud users)) section.
### Step 5: Generate report
After finalizing your data analysis, you can automatically generate an interactive visual HTML report for inspection of results together with parameters and code inside of the browser using:
```
snakemake --report report.zip
```
The resulting `report.zip` file can be passed on to collaborators, provided as a supplementary file in publications, or uploaded to a service like [Zenodo](https://zenodo.org/) in order to obtain a citable [DOI](https://en.wikipedia.org/wiki/Digital_object_identifier).

## Usage of workflow through SLURM (AAU Biocloud users)
Adjust the `config.yaml` files under both `config/` and `profiles/` accordingly, then simply run `snakemake --profile profiles/biocloud` or submit a SLURM job using the `slurm_submit.sbatch` example script.


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
* Replace `<owner>` and `<repo>` with the correct values in this `README.md` as well as in files under `.github/workflows/`.
* Replace `<snakemake_template>` with the workflow/project name (can be the same as `<repo>`) here as well as in the `environment.yml` and `slurm_submit.sbatch` files.
* The workflow will occur in the public Snakemake workflow catalog once the repository has been made public and the provided GitHub actions finish correctly. Then the link under "Usage" will point to the usage instructions if `<owner>` and `<repo>` were correctly set. If you don't want to publish the workflow just delete the `.github/workflows/` and `.template/` folders and `snakemake-workflow-catalog.yml`.
* Consider the license - [Choose a license](https://choosealicense.com/)
* DELETE this **TODO** section when finished with all of the above, and then start developing your workflow!
