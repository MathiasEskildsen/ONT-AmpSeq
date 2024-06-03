#!/bin/bash -l
# DESCRIPTION
#   Script used for the generation of OTU-table from amplicon-sequencing data obtained from ONT
#   This script should be run after "Amp-Seq_Statistics.sh"
# Authors
#   Mathias Eskildsen (mk20aj@bio.aau.dk)   
#   Patrick Skov Schaksen (insert mail)
#   
#
#   License    
# 
# 
# 
### DESCRIPTION -------------------------------------------------------------------------
# Usage message
USAGE="
-- ONT-AmpSeq: Workflow for generation of OTU table. Z-clustered OTU consensus polish with Medaka
usage: $(basename "$0" .sh) [-h] [-i path] [-o path] [-t value] [-j value] [-l value] [-u value] [-q value] [-m string] [-M] [-r string] [-db path]


where:
    -h Show this help message
    -i Directory path of unzipped raw ".fastq" files, from statistics workflow, /path/to/unzipped/raw/fastq/1_raw
    -o Path where directories and files should be stored
    -t Number of threads [default = 10]
    -j Number of parallel jobs [default = 1]
    -l Minimum length of reads (Check size distribution from statistics)
    -u Maximum length of reads (Check size distribution from statistics)
    -q Minimum q-score of reads (Check quality distribution of reads)
    -m ONT Device and basecalling model [default = r1041_e82_400bps_sup_v4.2.0]
    -M List available models for medaka polishing
    -r Method for taxonomic classification, SINTAX or blastn
    -d Full path to database for taxonomic classication, examples: /space/databases/midas/MiDAS4.8.1_20210702/output/FLASVs_w_sintax.fa or /space/databases/blast/nt_2022_07_28/nt
"

MODELS="r103_fast_g507, r103_fast_snp_g507, r103_fast_variant_g507, r103_hac_g507, r103_hac_snp_g507, r103_hac_variant_g507, r103_min_high_g345, r103_min_high_g360, r103_prom_high_g360, r103_prom_snp_g3210, r103_prom_variant_g3210, r103_sup_g507, r103_sup_snp_g507, r103_sup_variant_g507, r1041_e82_260bps_fast_g632, r1041_e82_260bps_fast_variant_g632, r1041_e82_260bps_hac_g632, r1041_e82_260bps_hac_v4.0.0, r1041_e82_260bps_hac_v4.1.0, r1041_e82_260bps_hac_variant_g632, r1041_e82_260bps_hac_variant_v4.1.0, r1041_e82_260bps_sup_g632, r1041_e82_260bps_sup_v4.0.0, r1041_e82_260bps_sup_v4.1.0, r1041_e82_260bps_sup_variant_g632, r1041_e82_260bps_sup_variant_v4.1.0, r1041_e82_400bps_fast_g615, r1041_e82_400bps_fast_g632, r1041_e82_400bps_fast_variant_g615, r1041_e82_400bps_fast_variant_g632, r1041_e82_400bps_hac_g615, r1041_e82_400bps_hac_g632, r1041_e82_400bps_hac_v4.0.0, r1041_e82_400bps_hac_v4.1.0, r1041_e82_400bps_hac_v4.2.0, r1041_e82_400bps_hac_variant_g615, r1041_e82_400bps_hac_variant_g632, r1041_e82_400bps_hac_variant_v4.1.0, r1041_e82_400bps_hac_variant_v4.2.0, r1041_e82_400bps_sup_g615, r1041_e82_400bps_sup_v4.0.0, r1041_e82_400bps_sup_v4.1.0, r1041_e82_400bps_sup_v4.2.0, r1041_e82_400bps_sup_variant_g615, r1041_e82_400bps_sup_variant_v4.1.0, r1041_e82_400bps_sup_variant_v4.2.0, r104_e81_fast_g5015, r104_e81_fast_variant_g5015, r104_e81_hac_g5015, r104_e81_hac_variant_g5015, r104_e81_sup_g5015, r104_e81_sup_g610, r104_e81_sup_variant_g610, r10_min_high_g303, r10_min_high_g340, r941_e81_fast_g514, r941_e81_fast_variant_g514, r941_e81_hac_g514, r941_e81_hac_variant_g514, r941_e81_sup_g514, r941_e81_sup_variant_g514, r941_min_fast_g303, r941_min_fast_g507, r941_min_fast_snp_g507, r941_min_fast_variant_g507, r941_min_hac_g507, r941_min_hac_snp_g507, r941_min_hac_variant_g507, r941_min_high_g303, r941_min_high_g330, r941_min_high_g340_rle, r941_min_high_g344, r941_min_high_g351, r941_min_high_g360, r941_min_sup_g507, r941_min_sup_snp_g507, r941_min_sup_variant_g507, r941_prom_fast_g303, r941_prom_fast_g507, r941_prom_fast_snp_g507, r941_prom_fast_variant_g507, r941_prom_hac_g507, r941_prom_hac_snp_g507, r941_prom_hac_variant_g507, r941_prom_high_g303, r941_prom_high_g330, r941_prom_high_g344, r941_prom_high_g360, r941_prom_high_g4011, r941_prom_snp_g303, r941_prom_snp_g322, r941_prom_snp_g360, r941_prom_sup_g507, r941_prom_sup_snp_g507, r941_prom_sup_variant_g507, r941_prom_variant_g303, r941_prom_variant_g322, r941_prom_variant_g360, r941_sup_plant_g610, r941_sup_plant_variant_g610
Default consensus:  r1041_e82_400bps_sup_v4.2.0
Default variant:  r1041_e82_400bps_sup_variant_v4.2.0"

# Process command-line options
while getopts 'i:o:t:j:l:u:q:m:r:d:Mh' OPTION; do
    case $OPTION in
        h) echo "$USAGE"; exit 1;;
        i) input_fastq=$OPTARG;;
        o) project_dir=$OPTARG;;
        j) JobNr=$OPTARG;;
        t) threads=$OPTARG;;
        m) model=$OPTARG;;
        l) read_cut_off_low=$OPTARG;;
        u) read_cut_off_high=$OPTARG;;
        q) q_score_cutoff=$OPTARG;;
        M)
            echo "$MODELS"
            exit 0;;
        r) method=$OPTARG;;
        d) database=$OPTARG;;
        :) printf "Missing argument for option -$OPTARG\n" >&2; exit 1;;
        \?) printf "Invalid option -$OPTARG\n" >&2; exit 1;;
    esac
done

# Check if the chosen method is blastn and prompt the user
if [ "$method" == "blastn" ]; then
    read -p "HEADS UP: You have chosen the method blastn for taxonomic classification. Make sure you have the blastn tool installed and the necessary database configured. Furthermore; blastn will only give you the best hit for the given OTU, however this is not always the correct classification. Do you want to proceed? (y/n): " choice
    case "$choice" in
        y|Y ) ;;
        n|N ) echo "Exiting."; exit 1;;
        * ) echo "Invalid choice. Exiting."; exit 1;;
    esac
fi

# Debugging



# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${input_fastq+x} ]; then echo "-i $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${project_dir+x} ]; then echo "-o $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${threads+x} ]; then echo "-t missing. Defaulting to 10 threads"; threads=10; fi;
if [ -z ${JobNr+x} ]; then echo "-j missing. Defaulting to 1 job"; JobNr=1; fi;
if [ -z ${model+x} ]; then echo "-m missing. Defaulting to model r1041_e82_400bps_sup_v4.2.0"; model=r1041_e82_400bps_sup_v4.2.0; fi;
if [ -z ${read_cut_off_low+x} ]; then echo "-l $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${read_cut_off_high+x} ]; then echo "-u $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${q_score_cutoff+x} ]; then echo "-q $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${database+x} ]; then echo "-d $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${method+x} ]; then echo "-r $MISSING"; echo "$USAGE"; exit 1; fi;

-set eu 

exec > >(tee -a $project_dir/joblog/joblog_main.txt) 2>&1

#Create directories 
mkdir -p $project_dir
mkdir -p $project_dir/joblog
mkdir -p $project_dir/2_filtering
mkdir -p $project_dir/3_fastq-fa
mkdir -p $project_dir/4_zotus
mkdir -p $project_dir/5_zotus_cat
mkdir -p $project_dir/6_medaka
mkdir -p $project_dir/7_medaka_relabel
mkdir -p $project_dir/7_medaka_relabel/combined_polished
mkdir -p $project_dir/8_OTUtable
mkdir -p $project_dir/8_OTUtable/97
mkdir -p $project_dir/8_OTUtable/99


#load mamba environment
mamba activate OTUtable

#Filter based upon information gained from "Amp-Seq_Statistics" workflow
input="$input_fastq" ##Input path to folder containing merged .fastq files
output="$project_dir/2_filtering" ## Output path

# Loop to create output directories for each input file in the output directory
# Use find to locate files with the pattern "_merged.fastq"
files=$(find "$input"/* -type f -name "*_merged.fastq")
for file in $files; do
    subdirectory_name=$(basename "$file" _merged.fastq)
    output_dir="${output}/${subdirectory_name}"
    echo "created directory for $subdirectory_name in $output"
    mkdir -p "$output_dir"
done

# Filter reads for q-score and length and trim adapters (22 nt's)
input="$input_fastq"
output="$project_dir/2_filtering"

files=($(find "$input"/* -type f -name "*_merged.fastq"))
run_filtering () {
    local file="$1"
    local output="$2" 
    local q_score="$3"
    local threads="$4"
    local min_length="$5"
    local max_length="$6"
    local base_name=$(basename "$file" _merged.fastq)
    local output_file="$output/${base_name}/${base_name}.fastq"

    echo "filtering from $min_length to $max_length and q $q_score in $base_name"
    # Fake cat as chopper can't work directly on previously concatenated file... whack
    cat $file | chopper \
        --minlength $min_length \
        --maxlength $max_length \
        --headcrop 22 \
        --tailcrop 22 \
        -q $q_score \
        -t $threads \
        > "$output_file"
    echo "done filtering from $min_length to $max_length in $base_name"
}

export -f run_filtering

parallel -j $JobNr run_filtering ::: "${files[@]}" ::: "$output" ::: "$q_score_cutoff" ::: $((threads / JobNr)) ::: "$read_cut_off_low" ::: "$read_cut_off_high"

echo "Finished filtering"

# Change headers '@' -> '>' i.e fastq -> fasta
input="$project_dir/2_filtering"
output="$project_dir/3_fastq-fa"

files=($(find "$input"/* -type f -name "*.fastq"))
for file in "${files[@]}"; do
    # Extract the base file name without extension
     base_name=$(basename "$file" .fastq)
     # Define the output file name
     output_file="$output/${base_name}_newheaders.fa"
    # Process the file and save the output in the output directory
    echo "changing '@' headers for $base_name to '>'"
    sed -n '1~4s/^@/>/p;2~4p' "$file" > "${output_file}"
done

echo "Finished changing headers"

# Z-Cluster OTU
input="$project_dir/3_fastq-fa"
output="$project_dir/4_zotus"

echo "Generating Z-clustered OTU's"
find "$input" -type f -iname '*.fa' | parallel -j $JobNr --eta vsearch --cluster_unoise {} --minsize 1 --threads $((threads / JobNr)) --centroids "${output}/{/.}_zotus.fa"
echo "Finished generating Z-clustered OTU's"


echo "Concatenating z-clustered OTU's" 
cat $project_dir/4_zotus/*.fa > $project_dir/5_zotus_cat/zotus_concatenated.fa
echo "Finished Concatenatig z-clustered OTU's"
# Z-Cluster OTU section done

# Medaka Polish section
echo "Starting polishing outputs saved are in $project_dir/6_medaka"

# Alignment part
input="$project_dir/4_zotus"
combined="$project_dir/5_zotus_cat/zotus_concatenated.fa"
output="$project_dir/6_medaka"

files=( $(ls ${input}/*.fa) ) 
run_mini_align() {
    local file="$1"
    local combined="$2"
    local output="$3"
    local threads="$4"

    file_prefix=$(basename "$file" | cut -d '_' -f 1)
    sub_dir="${output}/${file_prefix}"
    mkdir -p "$sub_dir"
    output_file="${sub_dir}/${file_prefix}_calls_to_draft"

    mini_align  -r "$file" -i "$combined" -m \
        -p "$output_file" \
        -t $threads
}
export -f run_mini_align

parallel -j $JobNr run_mini_align ::: "${files[@]}" ::: "$combined" ::: "$output" ::: $((threads / JobNr))

echo "run_mini_align finished"

# Medaka consensus algorithm (2nd step)
# Run medaka consensus algorithm 
input="$project_dir/6_medaka"
output="$project_dir/6_medaka"

files=( $(ls ${input}/*/*.bam) )
run_medaka_consensus() {
    local file="$1"
    local output="$2"
    local model="$3"
    local threads="$4"
    local filepath=$(dirname "$file")
    local sample=$(basename "$filepath")
    local output="${filepath}/${sample}_consensus.hdf"

    medaka consensus "$file" \
    "$output" \
    --model $model \
    --threads $threads
}
export -f run_medaka_consensus

parallel -j "$JobNr" run_medaka_consensus ::: "${files[@]}" ::: "$output" ::: "$model" ::: $((threads / JobNr))

# Medaka stich - creating consensus sequence of polished reads (3rd step)
input="$project_dir/6_medaka"
output="$project_dir/6_medaka"
zotus_dir="$project_dir/4_zotus"
draft_base="barcode"

files=( "${input}"/*/*.hdf )
for file in "${files[@]}"; do
    filepath=$(dirname "${file}")
    barcode=$(basename "${filepath}")
    draft="${zotus_dir}/${barcode}_newheaders_zotus.fa"

    medaka stitch "${file}" \
        "${draft}" \
        "${filepath}/polished.assembly.fasta"
done
# Medaka stich finished
# Renaming files to contain barcode#/BC# for later relabeling of reads
base_directory="$project_dir/6_medaka"
for dir in "$base_directory"/*; do 
    if [ -d "$dir" ]; then
    # Check if the directory contains a file named "polished.assembly.fasta"
    if [ -e "$dir/polished.assembly.fasta" ]; then
            #Construct new dir name 
            new_dir_name=$(basename "$dir")
            # Construct the new filename
            new_filename="${new_dir_name}_polished.assembly.fasta"
            # Construct the full path of the new file
            new_file_path="$dir/$new_filename"

            #rename the file
            mv "$dir/polished.assembly.fasta" "$new_file_path"

            echo echo "Renamed: $dir/polished.assembly.fasta -> $new_file_path"
        fi 
    fi
done

# Relabel reads using vsearch
input="$project_dir/6_medaka"
output="$project_dir/7_medaka_relabel"
files=$(find "$input"/* -type f -iname '*.fasta')
for file in $files; do
    sample=$(basename "$file")
 vsearch \
  --sortbysize $file \
  --relabel "$sample." \
  --threads $threads \
  --output "${output}/$(basename "${file%.*}")_relabel.fa"
done

cat $project_dir/7_medaka_relabel/*.fa > $project_dir/7_medaka_relabel/combined_polished/polished_comb_reads.fa

# 97% clustering using V-search
input="$project_dir/7_medaka_relabel/combined_polished/polished_comb_reads.fa"
output="$project_dir/8_OTUtable/97"

vsearch --cluster_fast $input -id 0.97 --threads $threads --relabel OTU_ --sizeout --otutabout $output/otutabe_0.97.tsv --centroids $output/otu_97.fa

# 99% clustering using V-search
input="$project_dir/7_medaka_relabel/combined_polished/polished_comb_reads.fa"
output="$project_dir/8_OTUtable/99"

vsearch --cluster_fast $input -id 0.99 --threads $threads --relabel OTU_ --sizeout --otutabout $output/otutabe_0.99.tsv --centroids $output/otu_99.fa

# Defining functions for Sintax or blastn database search:
taxonomy_sintax() {
    local input="$1"
    local output="$2"
    local database="$3"
    local threads="$4"
    local sintax_cutoff="$5"

    vsearch --sintax "$input" -db "$database" --tabbedout "$output/otu_cut.txt" \
    --threads "$threads" --quiet --sintax_cutoff "$sintax_cutoff" --strand both

    awk -F "\t" 'OFS="\t"{gsub(/;.*/," ",$1);print $1 $4}' "$output/otu_cut.txt" | awk -F ' ' '{print $1"\t"$2}' > "$output/tax2_cut.txt" 

    input2="$(dirname "$input")/otutabe_${sintax_cutoff}.tsv"
    awk 'BEGIN {FS=OFS="\t"} NR==FNR {hold[$1]=$2; next} {print $0, hold[$1]}' "$output/tax2_cut.txt" "$input2" > "$output/otutable_cut.tsv"
    
}

taxonomy_blastn() {
    local input="$1"
    local output="$2"
    local database="$3"
    local threads="$4"
    local clustering_percentage="$5"

    blastn -query "$input" -word_size 11 -max_target_seqs 1 -num_threads "$threads" \
    -evalue 1e-10 -outfmt "6 qseqid sseqid stitle evalue bitscore length pident" \
    -out "$output/otus_tax.txt" -db "$database"

    #Remove rows containing information regarding statistics
    awk -F'\t' '{print $1 "\t" $3}' $output/otus_tax.txt > $output/otus_tax_mod.txt
    #Remove size
    sed -i 's/;size=[0-9]\+\t/\t/' $output/otus_tax_mod.txt
    #Remove uncultured rows
    awk -i inplace -F'\t' '$2 !~ /^[uU]ncultured/' $output/otus_tax_mod.txt
    #Split words and add suffixes
    (echo -e "OTU ID\tgenus\tspecies\tnotes"; awk -F'\t' 'BEGIN {OFS="\t"} {split($2, words, " "); print $1, "g__"words[1], "s__"words[2], $2}' $output/otus_tax_mod.txt) > $output/otus_tax_mod1.txt
    #Rename headers
    sed -i '1 s/OTU ID/OTU_ID/' "$output/otus_tax_mod1.txt"
    sed -i '1 s/#OTU ID/OTU_ID/' "$output/otutabe_$clustering_percentage.tsv"
    #Merge datasets to final OTUtable
    python3 <<EOF
    import pandas as pd
    import sys       

    otus_tax = pd.read_csv("$output/otus_tax_mod1.txt", sep='\t')
    otutabe = pd.read_csv("$output/otutabe_$clustering_percentage.tsv", sep='\t')
    merged_data = pd.merge(otutabe, otus_tax, on='OTU_ID', how='right')
    merged_data.to_csv("$output/OTUtable_$clustering_percentage.tsv", sep='\t', index=False)
EOF
}


input_97="$project_dir/8_OTUtable/97/otu_97.fa"
output_97="$project_dir/8_OTUtable/97"
input_99="$project_dir/8_OTUtable/99/otu_99.fa"
output_99="$project_dir/8_OTUtable/99"

if [ "$method" = "SINTAX" ]; then
    taxonomy_sintax "$input_97" "$output_97" "$database" "$threads" 0.97
    taxonomy_sintax "$input_99" "$output_99" "$database" "$threads" 0.99
else
    if [ "$method" = "blastn" ]; then
    taxonomy_blastn "$input_97" "$output_97" "$database" "$threads" 0.97
    taxonomy_blastn "$input_99" "$output_99" "$database" "$threads" 0.99
else
    echo "Invalid method. Please enter SINTAX or blastn"
    exit 1
    fi
    exit 0
fi

# 2nd tax OTU table with reads more than 1 time and more than 10 times

input="$project_dir/8_OTUtable/97/otutable_cut.tsv"
output="$project_dir/8_OTUtable/97"
awk 'FNR == 1{print;next}{sum=0; for(i=3; i<=NF; i++) {sum += $i} if (sum > 1) print}' $input > $output/otutable_tax_filtered_1read_97.tsv

awk 'FNR == 1{print;next}{sum=0; for(i=3; i<=NF; i++) {sum += $i} if (sum > 10) print}' $input > $output/otutable_tax_filtered_10read_97.tsv


# 2nd tax OTU table with reads more than 1 time and more than 10 times

input="$project_dir/8_OTUtable/99/otutable_cut.tsv"
output="$project_dir/8_OTUtable/99"
awk 'FNR == 1{print;next}{sum=0; for(i=3; i<=NF; i++) {sum += $i} if (sum > 1) print}' $input > $output/otutable_tax_filtered_1read_99.tsv

awk 'FNR == 1{print;next}{sum=0; for(i=3; i<=NF; i++) {sum += $i} if (sum > 10) print}' $input > $output/otutable_tax_filtered_10read_99.tsv

##### V-search tax and clustering finished 

# Python script content
python_script=$(cat <<END
import csv
import sys

# Rest of the Python script
input_file = sys.argv[1]
output_file = sys.argv[2]
new_column_headers = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

# Create a dictionary to map field names to new column values
field_mapping = {
    "d": "kingdom",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}

# Function to extract and map the values to the new columns
def extract_and_map_values(input_row):
    last_column = input_row.pop()  # Remove the last column
    field_values = {field: "" for field in new_column_headers}

    parts = last_column.split(",")
    
    for part in parts:
        field_prefix, value = part.split(":") if ":" in part else (None, part)
        new_column = field_mapping.get(field_prefix)
        if new_column:
            field_values[new_column] = f"{new_column}__{value}"

    new_row = input_row + [field_values[field] for field in new_column_headers]

    return new_row

# Open the input and output files
with open(input_file, mode="r") as infile, open(output_file, mode="w", newline="") as outfile:
    reader = csv.reader(infile, delimiter="\t")
    writer = csv.writer(outfile, delimiter="\t")

    # Write the headers with the new columns
    headers = next(reader)[:-1] + new_column_headers
    writer.writerow(headers)

    # Process and write the data
    for row in reader:
        modified_row = extract_and_map_values(row)
        writer.writerow(modified_row)


END
)

# End of python script content

# Execute the Python script for 97 1 read
input_file="$project_dir/8_OTUtable/97/otutable_tax_filtered_1read_97.tsv"
output_file="$project_dir/8_OTUtable/97/otutable_tax_fixed_filtered_1read_97.tsv"
python3 -c "$python_script" "$input_file" "$output_file"

# Execute the Python script for 97 10 read
input_file="$project_dir/8_OTUtable/97/otutable_tax_filtered_10read_97.tsv"
output_file="$project_dir/8_OTUtable/97/otutable_tax_fixed_filtered_10read_97.tsv"
python3 -c "$python_script" "$input_file" "$output_file"

# Execute the Python script for 99 1 read
input_file="$project_dir/8_OTUtable/99/otutable_tax_filtered_1read_99.tsv"
output_file="$project_dir/8_OTUtable/99/otutable_tax_fixed_filtered_1read_99.tsv"
python3 -c "$python_script" "$input_file" "$output_file"

# Execute the Python script for 99 10 read
input_file="$project_dir/8_OTUtable/99/otutable_tax_filtered_10read_99.tsv"
output_file="$project_dir/8_OTUtable/99/otutable_tax_fixed_filtered_10read_99.tsv"
python3 -c "$python_script" "$input_file" "$output_file"

modifications=(
  's/kingdom__\([^[:space:]]*\)/k__\1/g'
  's/phylum__\([^[:space:]]*\)/p__\1/g'
  's/class__\([^[:space:]]*\)/c__\1/g'
  's/order__\([^[:space:]]*\)/o__\1/g'
  's/family__\([^[:space:]]*\)/f__\1/g'
  's/genus__\([^[:space:]]*\)/g__\1/g'
  's/species__\([^[:space:]]*\)/s__\1/g'
)
for subdir in "$project_dir"/8_OTUtable/*;do
    if [ -d "$subdir" ]; then
        for input in "$subdir"/*_fixed_*; do
            if [ -f "$input" ]; then
                for modification in "${modifications[@]}"; do
                echo "Applying modifications $modification to $input"
                sed -i -e "$modification" "$input"
                done
            fi
        done
    fi
done



## Prepare tax table for phyloseq
## Raw vsearch output is 1 of the inputs for phyloseq ex. ($project_dir/8_OTUtable/97/otutabe_0.97.tsv) for each clustering percentage
input="$project_dir/8_OTUtable"
output="$project_dir/8_OTUtable"
# Use a loop to directly iterate over the results of find
find "$input" -type f -name "otutable_tax_fixed_filtered_*read_*.tsv" -print0 | while IFS= read -r -d '' file; do
    # Check if the file is a regular file
    if [ -f "$file" ]; then
        # Get the base file name without extension
        base_name=$(basename "$file" .tsv)
        # Output directory
        output_dir=$(dirname "$file")
        # Define the output file name
        output_file="$output_dir/${base_name}_phyloseq_tax.tsv"
        # Create the output file
        # Run the Python script
        python3 <<EOF
import csv
import sys

input_file = "$file"
output_file = "$output_file"
with open(input_file, 'r') as csvinput:
    with open(output_file, 'w') as csvoutput:
        reader = csv.reader(csvinput, delimiter='\t')
        writer = csv.writer(csvoutput, delimiter='\t')

        for row in reader:
            if 'kingdom' in row:
                kingdom_index = row.index('kingdom')
            writer.writerow(row[:1] + row[kingdom_index:])

print("Processing completed. Output written to", output_file)
EOF
    fi
done


##Remove suffixes from phyloseq OTU-table
project_dir=/home/bio.aau.dk/mk20aj/Projects/zymo_test_data/test_outdir
input=$project_dir/8_OTUtable
output=$project_dir/8_OTUtable

# Use a loop to directly iterate over the results of find
find "$input" -type f -name "*phyloseq_tax.tsv" -print0 | while IFS= read -r -d '' file; do
    # Get the base file name without extension
    base_name=$(basename "$file" .tsv)
    # Output directory
    output_dir=$(dirname "$file")
    # Define the output file name
    output_file="$output_dir/${base_name}_modified.tsv"
    # Get the header line
    header=$(head -n 1 "$file")
    # Write the header line to the output file
    echo "$header" > "$output_file"
    # Process the rest of the file
    tail -n +2 "$file" | awk -F'\t' 'BEGIN {OFS="\t"} {printf $1; for (i=2; i<=NF; i++) printf "\t%s", substr($i, 4); printf "\n"}' >> "$output_file"

    echo "Processing completed for $file. Output written to $output_file"
done
###### TEST SHIT AT BOTTOM 


