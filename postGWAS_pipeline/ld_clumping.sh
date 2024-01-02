# perform LD clumping on the GWAS results

# set the parameters
in_dir=$HOME/saige_pipe_test/
out_dir=$HOME/bugbank_data/postgwas_regenie/

# check if parameter is provided
if [ -z "$1" ]
then
    echo "please provide the stem of the summary files"
    exit
fi
stem=$1

# list the GWAS result files
filtered_snp_fs=$(ls $in_dir/sumary.${stem}.*.txt.gz)
# iterate through the files and perform clumping
for filtered_snp_f in $filtered_snp_fs
do
    plink2 
done

