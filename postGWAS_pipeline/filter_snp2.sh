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
STEM=$1
# relevant column positions (0-based index)
col_n=8
col_a1freq=6
col_info=7
col_pval=13

# Filtering thresholds
min_maf=0.001
max_maf=0.999
min_info=0.3
min_log10p=7.301

# get the list of summary files
summary_files=$(ls $in_dir/summary.${STEM}.*.txt.gz)
#summary_files=("/users/bag/hlq763/saige_pipe_test/summary.05062023_sgss_species.Escherichia_coli.species.all.txt.gz")
# loop through the summary files
for summary_file in $summary_files
do
    # check if there are more than 40 background jobs running
    # if so wait for 10 seconds and recheck again
    while [ $(ps | grep sh | wc -l) -gt 40 ]
    do
        echo "too many jobs running, sleeping for 10 seconds..."
        sleep 10
    done

    echo "filtering $summary_file..."
    # outfile name
    outname=$(echo $summary_file | sed 's/^.*\///g' | sed 's/txt.gz/snp/g' | sed 's/summary/filtered/g')
    out_file=${out_dir}/${outname}

    # print the filtering criteria into header
    echo "# min_maf=${min_maf}, max_maf=${max_maf}, min_info=${min_info}, min_log10p=${min_log10p}" > $out_file
    # print the header from the summary file
    zcat $summary_file | head -n 1 >> $out_file
    # filter on each relevant column
    (
        zcat "$summary_file" |
        awk -v col_a1freq="$col_a1freq" \
            -v min_maf="$min_maf" \
            -v max_maf="$max_maf" \
            -v col_info="$col_info" \
            -v min_info="$min_info" \
            -v col_pval="$col_pval" \
            -v min_log10p="$min_log10p" \
            '{
                if (($(col_info) >= min_info && 
                    $(col_a1freq) >= min_maf && 
                    $(col_a1freq) <= max_maf &&  
                    $(col_pval) >= min_log10p)) print $0
            }' >> "$out_file"
        echo "$out_file done"
    )&
done