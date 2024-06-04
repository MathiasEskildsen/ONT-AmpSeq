#!/bin/bash -l
# DESCRIPTION
#   Script used for the generation of statistics and file management of amplicon-sequencing
# Authors
#   Mathias Eskildsen (mk20aj@bio.aau.dk)
#   Patrick Skov Schaksen (insert mail)
#
#   license GNU General Public License
#
# TO DO
#
### DESCRIPTION -------------------------------------------------------------------------
# Last modified 15-01-2024

USAGE="
-- insert full pipeline name: Nanopore Statistics with NanoPlot
usage: $(basename "$0" .sh) [-h] [-o path] [-i path] [-t value] [-j value]

where:
    -h Show this help message.
    -o Output directory path for the processed files. The script automatically creates a directory based on the specified user-defined output name.
    -i Full path to the input files. If using the default fastq folders created by MinKNOW, a example path could be: /Full/Path/to/nanopore_data/ONT_RUN_ID/fastq_pass.  
    -j Number of parallel jobs [default = 1]
    -t Number of threads [default = 1]
    Important note:
    Remember to activate your conda environment, containing nanoplot version 1.42.0, before running the script.
    If installed through stats.yml, activate the environment with mamba activate stats.
    The total number of threads used by the script is number of parallel jobs x number of threads, e.g. -j 2 -t 10 will use 20 threads.
"
# Process command-line options
while getopts 'o:i:t:j:h' OPTION; do
    case $OPTION in
        h) echo "$USAGE"; exit 1;;
        o) project_dir=$OPTARG;;
        i) input_fastq=$OPTARG;;
        j) JobNr=$OPTARG;;
        t) threads=$OPTARG;;
        :) printf "Missing argument for option -$OPTARG\n" >&2; exit 1;;
        \?) printf "Invalid option -$OPTARG\n" >&2; exit 1;;
    esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${project_dir+x} ]; then echo "-o $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${threads+x} ]; then echo "-t missing. Defaulting to 10 threads"; threads=10; fi;
if [ -z ${JobNr+x} ]; then echo "-j missing. Defaulting to 1 job"; JobNr=1; fi;
if [ -z ${input_fastq+x} ]; then echo "-i $MISSING"; echo "$USAGE"; exit 1; fi; 

#Create directories 
mkdir -p $project_dir
mkdir -p $project_dir/stats
mkdir -p $project_dir/fastqs
mkdir -p $project_dir/joblog 
exec > >(tee -a $project_dir/joblog/joblog_stat.txt) 2>&1

echo "project_dir: $project_dir"
echo "input_fastq: $input_fastq"
echo "JobNr": $JobNr
echo "Threads: $threads"

#¤ Move unzip and concatenate passed .fastq reads
## Move unzip
### Input files. Directory path leading to passed fastq.gz files
input=$input_fastq
output="$project_dir/fastqs"
files=($(find "$input" -type f \( -name "*.gz" -o -name "*.fastq" -o -name "*.fasta" \)))
for file in "${files[@]}"; do
    parent_dir=$(dirname "$file")
    # Check if the file is in a subdirectory
    if [[ "$parent_dir" != "$input" ]]; then
        # If the file is in a subdirectory, create the subdirectory in the output directory
        subdirectory_name=$(basename "$parent_dir")
        output_dir="${output}/${subdirectory_name}"
        mkdir -p "$output_dir"
        if gzip -t $file 2>/dev/null; then
            echo "Decompressing $file"
            gunzip -c "$file" > "$output_dir/$(basename "$file" ".gz")"
        else
            echo "Moving $file"
            cp "$file" "$output_dir/$(basename "$file")"
        fi
    else
        filename=$(basename -- "$file")
        filename="${filename%.*}"
        filename="${filename%.*}"
        # Create the output directory based on the filename name
        output_dir="${output}/${filename}"
        mkdir -p "$output_dir"
        if gzip -t $file 2>/dev/null; then
            echo "Decompressing $file"
            gunzip -c "$file" > "$output_dir/$(basename "$file" ".gz")"
        else
            echo "Moving $filename"
            cp "$file" "$output_dir/$(basename "$file")"
        fi
    fi
done

echo "Finished moving and unzipping"

## Concatenate .fastq files
base_dir="$project_dir/fastqs"
# Iterate over each subdirectory
for sub_dir in "$base_dir"/*; do
    if [ -d "$sub_dir" ]; then
        sub_dir_name=$(basename "$sub_dir")
        echo "Processing files in: $sub_dir_name"
        # Concatenate all files in the subdirectory
        cat "${sub_dir}"/* > "${sub_dir}/${sub_dir_name}_merged.fastq"
        echo "Merged all files in $sub_dir_name"
    fi
done
# Finished unzipping, moving and concatenating passed .fastq files


# Start of statistics workflow 
input="$project_dir/fastqs"
output="$project_dir/stats"
#Create output directories for each input file in the output directory
# Use find to locate files with the pattern "_concatenated.fastq"
files=$(find "$input"/* -type f -name "*_merged.fastq")
for file in $files; do
    # Extract the subdirectory name from the file path
    subdirectory_name=$(basename "$file" _merged.fastq)
    # Construct the output directory path based on the subdirectory name
    output_dir="${output}/${subdirectory_name}"
    # Create the output directory if it doesn't exist
    echo "created directory for $subdirectory_name in $output"
    mkdir -p "$output_dir"
done

input_dir="$project_dir/fastqs"
output_dir="$project_dir/stats"

files=( $(ls ${input_dir}/*/*_merged.fastq) )
# Create function for NanoPlot
statistics() {
    local file="$1"
    local threads="$2"
    local output_dir="$3"
    
    # Extract subdirectory name
    local subdirectory_name=$(basename "$(dirname "$file")")
    
    # Create output directory within the specified output_dir
    local output="$output_dir/$subdirectory_name"
    
    echo "Creating plots and stats for $subdirectory_name in $output"
    # Debug statement to print the actual file being processed
    echo "Processing file: $file"

    NanoPlot --fastq "$file" --plots dot -t "$threads" -o "$output"
    
    echo "Finished plots and stats for $subdirectory_name"
}

export -f statistics

parallel -j "$JobNr" statistics ::: "${files[@]}" ::: $((threads / JobNr)) ::: "$output_dir" 
module purge


echo "Check your amplicon size and quality before setting parameters for chopper filtering, used in the snakemake workflow. Filtering settings can be changed in the config file. \"config/config.yaml\""