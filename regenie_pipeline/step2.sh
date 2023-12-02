#!/bin/bash
#SBATCH -D ./
#SBATCH -A bag.prj
#SBATCH -p short
#SBATCH --cpus-per-task=15
#SBATCH --array=1-24

#TMPDIR=/well/bag/clme1992/tmp
#export TMPDIR

echo "job array ID: "${SLURM_ARRAY_TASK_ID}
CHR=${SLURM_ARRAY_TASK_ID}
STEM=$1
STRAT=$2
UKBDEDIR=$3
UKBDIR=$4
WKDIR=$5

SAMP=$UKBDEDIR/ukb53100_imp_any_chr_v3_s487296.sample

cd $WKDIR
if [ "$CHR" = "23" ]; then
	CHR="X"
	SAMP=$UKBDEDIR/ukb53100_imp_chrX_v3_s486645.sample
fi

if [ "$CHR" = "24" ]; then
	CHR="XY"
	SAMP=${UKBDEDIR}/ukb53100_imp_chrXY_v3_s486331.sample   
fi

/well/bag/clme1992/regenie_v3.2.8.gz_x86_64_Centos7_mkl \
  --threads 15 \
  --step 2 \
  --bgen ${UKBDIR}/v3/imputation/ukb_imp_chr${CHR}_v3.bgen \
  --ref-first \
  --sample ${SAMP} \
  --phenoFile ${WKDIR}/regenie_test.pheno \
  --covarFile ${WKDIR}/regenie_test.cov \
  --bt \
  --firth --approx --pThresh 0.01 \
  --pred step1_${STEM}_${STRAT}_pred.list \
  --bsize 400 \
  --out ${WKDIR}/step2.${STEM}_${STRAT}_chr${CHR} \
  --minMAC 30

if [ "$CHR" = "1" ]; then
  tail -n +2 step2.${STEM}_${STRAT}_chr${CHR}_pheno.regenie | gzip > step2.${STEM}_${STRAT}_chr${CHR}.gz
else
  tail -n +1 step2.${STEM}_${STRAT}_chr${CHR}_pheno.regenie | gzip > step2.${STEM}_${STRAT}_chr${CHR}.gz
fi
rm -f step2.${STEM}_${STRAT}_chr${CHR}_pheno.regenie
