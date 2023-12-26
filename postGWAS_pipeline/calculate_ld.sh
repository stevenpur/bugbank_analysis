# To determine indepedent GWAS hits, we need to do LD clumping.
# This is for calculating LD for the significant snps from the post-GWAS analysis

# set the parameters
in_dir=$HOME/bugbank_data/postgwas_regenie/sig_snp
out_dir=$HOME/bugbank_data/postgwas_regenie/sig_snp_ld
bgen_dir=/well/ukbb-wtchg/v3/imputation/
# check if parameter is provided
if [ -z "$1" ]
then
    echo "please provide the stem of the summary files"
    exit
fi
STEM=$1

# list the files of nominally significant snps
filtered_snp_fs=$(ls $in_dir/filtered.${STEM}.*)

# calculating LD requires us to read the original .bgen files
# this process is very slow, so it is more idea if we pull all snps we need in one go
# Therefore, we pool the snps from the filtered snp files together
for filtered_snp_f in $filtered_snp_fs
do
    awk '{print $1" "$2" "$3}' $filtered_snp_f >> ${STEM}_all_filtered_snps.lst
done
