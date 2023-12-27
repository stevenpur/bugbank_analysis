# To determine indepedent GWAS hits, we need to do LD clumping.
# This is for calculating LD for the significant snps from the post-GWAS analysis

# set the parameters
in_dir=$HOME/bugbank_data/postgwas_regenie/sig_snp
out_dir=$HOME/bugbank_data/postgwas_regenie/sig_snp_ld
bgen_dir=/well/ukbb-wtchg/v3/imputation/
snpid_col=3
# temporary output file for pooling all filtered snps together
pool_snp_prefix=$out_dir/tmp_poolsnps_lst

# check if parameter is provided
if [ -z "$1" ]
then
    echo "please provide the stem of the summary files"
    exit
fi
STEM=$1

# list the files of nominally significant snps
filtered_snp_fs=$(ls $in_dir/filtered.${STEM}.*)
echo $filtered_snp_fs

# calculating LD requires us to read the original .bgen files
# this process is very slow, so it is more idea if we pull all snps we need in one go
# Therefore, we pool the snps from the filtered snp files together

for filtered_snp_f in $filtered_snp_fs
do
    # skip the commenting lines start with "#"
    # after that, skip the header line
    # then, print the snpid and split the result by chromosome
    cat $filtered_snp_f | 
        awk '{if ($0 !~ /^#/) print $0}' | 
        tail -n +2 |
        awk -v pool_snp_prefix=${pool_snp_prefix} -v snpid_col=${snpid_col} '{print $0 >> pool_snp_prefix".chr"$1".lst"}'
done
