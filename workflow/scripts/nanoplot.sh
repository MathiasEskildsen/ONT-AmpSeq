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
input="$input_fastq"
output="$project_dir/fastqs"
mkdir -p "$output"

# Function to process a barcode directory or individual file
process_files() {
    local path="$1"
    local output_dir="$2"

    if [ -d "$path" ]; then
        # Process barcode directory
        local barcode_name=$(basename "$path")

        # Skip the "unclassified" directory
        if [[ "$barcode_name" == "unclassified" ]]; then
            return
        fi

        local output_file="${output_dir}/${barcode_name}.fastq"
        local files=($(find "$path" -type f \( -name "*.gz" -o -name "*.fastq" -o -name "*.fasta" \)))

        for file in "${files[@]}"; do
            if gzip -t "$file" 2>/dev/null; then
                gunzip -c "$file" >> "$output_file"
            else
                cat "$file" >> "$output_file"
            fi
        done
    else
        # Process individual file
        local filename=$(basename "$path")
        if gzip -t "$path" 2>/dev/null; then
            gunzip -c "$path" > "${output_dir}/${filename%.gz}"
        else
            cp "$path" "$output_dir/$filename"
        fi
    fi
}

export -f process_files

if compgen -G "$input/*/" > /dev/null; then
    # Get all directories and files in the input
    dirs_or_files=($(find "$input" -mindepth 1 -maxdepth 1))
else
    # Get only files if no subdirectories exist
    dirs_or_files=($(find "$input" -type f \( -name "*.gz" -o -name "*.fastq" -o -name "*.fasta" \)))
fi

# Run the processing in parallel
parallel -j "$JobNr" process_files ::: "${dirs_or_files[@]}" ::: "$output"

echo "Finished moving and unzipping"

combined_file="${output}/barcodes_combined.fastq"
cat "${output}"/*.fastq > "$combined_file"
echo "Concatenation complete."

# Finished unzipping, moving and concatenating passed .fastq files

project_dir="${project_dir%/}"

input_dir="$project_dir/fastqs"
output_dir="$project_dir/stats"

# Correctly list .fastq files in the input directory
files=( $(ls "${input_dir}"/*.fastq) )

# Create function for NanoPlot
statistics() {
    local file="$1"
    local threads="$2"
    local output_base_dir="$3"
    
    # Extract file name without extension
    local filename_no_ext=$(basename "$file" .fastq)

    # Define output directory using output_base_dir and filename
    local output="$output_base_dir/$filename_no_ext"

    NanoPlot --fastq "$file" --plots dot -t "$threads" -o "$output"
    
    echo "Finished plots and stats for $filename_no_ext"
}

export -f statistics

parallel -j "$JobNr" statistics ::: "${files[@]}" ::: $((threads)) ::: "$output_dir"


# remove the combined file as it should not be run in the ONT-AmpSeq pipeline
rm -f "$project_dir/fastqs/barcodes_combined.fastq"

echo "Check your amplicon size and quality before setting parameters for chopper filtering, used in the snakemake workflow. Filtering settings can be changed in the config file. \"config/config.yaml\""
