#!/bin/bash

# Apply filter to extract SNPs for further investigation after GWAS

# I/O parameters
in_dir=$HOME/saige_pipe_test/
out_dir=$HOME/bugbank_data/postgwas_regenie/sig_snp

# check if parameter is provided
if [ -z "$1" ]
then
    echo "please provide the stem of the summary files"
    exit
fi

# relevant column positions (0-based index)
col_n=7
col_a1freq=5
col_info=6
col_pval=12

# Filtering thresholds
min_mac=300
min_info=0.3
min_log10p=7.301

# define a function to read file line by line and apply filter
# Define a function for processing a file
apply_filter() {
    local file=$1
    local out_file=$2
    # skip the first line
    read -r header < "$file"
    echo "$header" > "$output_file"
    # read the file line by line
    while IFS=$'\t' read -r -a line; do
        local n=${line[$col_n]}
        local a1freq=${line[$col_a1freq]}
        local info=${line[$col_info]}
        local pval=${line[$col_pval]}

        # Check if A1FREQ represents the minor allele frequency
        if (( $(echo "$a1freq > 0.5" | bc -l) )); then
            a1freq=$(echo "1 - $a1freq" | bc)
        fi
        # Calculate MAC using bc
        local mac=$(echo "2 * $n * $a1freq" | bc)

        # Apply filtering criteria
        if (( $(echo "$mac >= $min_mac" | bc -l) )) && \
           (( $(echo "$info >= $min_info" | bc -l) )) && \
           (( $(echo "$pval >= $min_log10p" | bc -l) )); then
            echo "${line[*]}" >> "$out_file"
        fi
    done < "$file"
    echo "$out_file done"
}


# get the list of summary files
summary_files=$(ls $in_dir/summary.${STEM}.*.txt.gz)

# loop through the summary files
for summary_file in $summary_files
do
    # make sure there is no more than 40 backgound jobs running before running a new command
    while [ $(ps | grep sh | wc -l) -gt 40 ]
    do
        echo "too many jobs running, sleeping for 5 seconds..."
        sleep 10
    done
    # define the output file
    f_out=$(echo $summary_file | sed 's/summary/filtered/g' | 
        sed "s/txt.gz/mac${min_mac}info${min_info}p${min_log10p}.txt/g")
    # apply snp filters to the summary file
    echo "applying filter to $summary_file..."
    apply_filter $summary_file $f_out &
done


