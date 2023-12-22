# get the SNPs that are at least nominally significant from the summary files

# set the parameters
in_dir=$HOME/saige_pipe_test/
out_dir=$HOME/bugbank_data/postgwas_regenie
# check if parameter is provided
if [ -z "$1" ]
then
    echo "please provide the stem of the summary files"
    exit
fi
STEM=$1

# get the list of summary files
summary_files=$(ls $in_dir/summary.${STEM}.*.txt.gz)

# loop through the summary files and get the SNPs that are at least nominally significant
for summary_file in $summary_files
do
    f_1e6=$(echo $summary_file | sed 's/^.*\///g' | sed 's/txt.gz/1e6_snp.txt/g')
    # make sure there is no more than 40 backgound jobs running before running a new command
    while [ $(ps | grep sh | wc -l) -gt 40 ]
    do
        echo "too many jobs running, sleeping for 5 seconds..."
        sleep 10
    done
    # run commanp to keep the rows where -log10(p) >= 6, run it in the backgournd
    (
        echo "retrieving 1e6 snps for $f_1e6...."
        zcat $summary_file | awk  '{if ($13 >= 6) print $0}' > ${out_dir}/${f_1e6}
        echo "$f_1e6 done"
    )&
done

echo "all files done"
