#!/usr/bin/bash -l
#SBATCH --job-name=<snakemake_template>
#SBATCH --output=/home/bio.aau.dk/mk20aj/joblog/job_%j_%x.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=general
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mhes@bio.aau.dk

# Exit on first error and if any variables are unset
set -eu

# Activate conda environment with only snakemake
conda activate snakemake_test

# Start workflow using ressources defined in the profile. Snakemake itself 
# requires nothing, 1CPU + 1G mem is enough

# Render a DAG to visualize the workflow (optional)
#snakemake --dag | dot -Tsvg > results/dag.svg

# Main workflow
snakemake --profile profiles/biocloud --config taxonomy_blast=true

# Generate a report once finished (optional)
#snakemake --report results/report.html
