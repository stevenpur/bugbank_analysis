#!/bin/bash
#SBATCH -D ./
#SBATCH -A bag.prj
#SBATCH -p short
#SBATCH --cpus-per-task=15

#TMPDIR=/well/bag/clme1992/tmp
#export TMPDIR

echo "job array ID: "${SLURM_ARRAY_TASK_ID}
CHR=23
#STEM=$1
#STRAT=$2
#UKBDEDIR=$3
#UKBDIR=$4
#WKDIR=$5

SAMP=$UKBDEDIR/ukb53100_imp_any_chr_v3_s487296.sample

cd /well/bag/clme1992/saige_pipe_test
#cd $WKDIR
#if [ "$CHR" = "23" ]; then
#	CHR="X"
#	SAMP=$UKBDEDIR/ukb53100_imp_chrX_v3_s486645.sample
#fi

#if [ "$CHR" = "24" ]; then
#	CHR="XY"
#	SAMP=${UKBDEDIR}/ukb53100_imp_chrXY_v3_s486331.sample   
#fi

/well/bag/clme1992/regenie_v3.3.gz_x86_64_Centos7_mkl \
  --threads 15 \
  --step 2 \
  --bgen /well/ukbb-wtchg/v3/imputation/ukb_imp_chrX_v3.bgen \
  --ref-first \
  --sample /well/bag/wilson/ukb/ukb53100_imp_chrX_v3_s486645.sample \
  --phenoFile /well/bag/clme1992/saige_pipe_test/ukb41482.bd.gwas-pheno.05062023_hes_species_regenie.Streptococcus_pneumoniae.species.txt.gz \
  --covarFile /well/bag/clme1992/saige_pipe_test/ukb41482.bd.gwas-covar.05062023_hes_species_regenie.Streptococcus_pneumoniae.species.txt.gz \
  --bt \
  --firth \
  --approx \
  --pThresh 0.01 \
  --pred step1_05062023_hes_species_regenie_Streptococcus_pneumoniae.species_pred.list \
  --bsize 400 \
  --out /well/bag/clme1992/saige_pipe_test/step2.05062023_hes_species_regenie_Streptococcus_pneumoniae.species_chrX \
  --minMAC 100


#if [ "$CHR" = "1" ]; then
#  tail -n +1 step2.${STEM}_${STRAT}_chr${CHR}_pheno.regenie | gzip > step2.${STEM}_${STRAT}_chr${CHR}_test.gz
#else
#  tail -n +2 step2.${STEM}_${STRAT}_chr${CHR}_pheno.regenie | gzip > step2.${STEM}_${STRAT}_chr${CHR}_test.gz
#fi
#rm -f step2.${STEM}_${STRAT}_chr${CHR}_pheno.regenie
