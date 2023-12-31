# To determine indepedent GWAS hits, we need to do LD clumping.
# This is for calculating LD for the significant snps from the post-GWAS analysis

# set the parameters
in_dir=$HOME/bugbank_data/postgwas_regenie/sig_snp
out_dir=$HOME/bugbank_data/postgwas_regenie/sig_snp_ld
bgen_dir=/well/ukbb-wtchg/v3/imputation
samp_dir=/well/bag/wilson/ukb
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


# loop through the autosomes and calculate LD
for chrs in {1..22}
do
    # check if the file exist, create empty file and skip this chromosome
    if [ ! -f $pool_snp_prefix.chr${chrs}.lst ]
    then
        touch $pool_snp_prefix.chr${chrs}.lst
        continue
    fi
    # get the corresponding genotype file
    bgen_f=$bgen_dir/ukb_imp_chr${chrs}_v3.bgen
    # to avoid running out of disk space, run at most 5 chromosomes at the same time
    while [ $(ps | grep plink | wc -l) -gt 5 ]
    do
        echo "too many jobs running, sleeping for 300 seconds..."
        sleep 5
    done
    # use plink to calculate the pairwise LD, run in background for parallelization
    (
        plink2 --bgen $bgen_f ref-first\
            --sample ${samp_dir}/ukb53100_imp_chr1_v3_s487296.sample \
            --extract $pool_snp_prefix.chr${chrs}.lst \
            --r2-phased \
            --ld-window-kb 1000 \
            --out $out_dir/ld.${stem}.filtered_snps.chr${chrs}
    )&
done

# chromosome X is a bit different
# we need to extract the snps from both chrX and chrXY from the .bgen files
bgen_x=$bgen_dir/ukb_imp_chrX_v3.bgen
bgen_xy=$bgen_dir/ukb_imp_chrXY_v3.bgen
# extract from bgen_x
(
    plink2 --bgen $bgen_x ref-first\
        --sample ${samp_dir}/ukb53100_imp_chrX_v3_s486645.sample  \
        --extract $pool_snp_prefix.chrX.lst \
        --make-bed \
        --out $out_dir/filtered_snps.${stem}.chrX
)
# extract from bgen_xy
(
    plink2 --bgen $bgen_xy ref-first\
        --sample ${samp_dir}/ukb53100_imp_chrXY_v3_s486331.sample \
        --extract $pool_snp_prefix.chrX.lst \
        --make-bed \
        --out $out_dir/filtered_snps.${stem}.chrXY
)
# merge the two files
(
    plink2 --bfile $out_dir/filtered_snps.chrX \
        --bmerge $out_dir/filtered_snps.chrXY.bed $out_dir/filtered_snps.chrXY.bim $out_dir/filtered_snps.chrXY.fam \
        --make-bed \
        --out $out_dir/filtered_snps.${stem}.chrX_merged
)
# calculate LD
(
    plink2 --bfile $out_dir/filtered_snps.chrX \
        --r2-phased \
        --ld-window-kb 1000 \
        --out $out_dir/ld.${stem}.filtered_snps.chrX
)
