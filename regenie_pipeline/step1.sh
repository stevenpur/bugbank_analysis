#!/bin/bash
#SBATCH -D ./
#SBATCH -A bag.prj
#SBATCH -p short

STEM=$1
STRAT=$2
NCORES=$3
WDIR=$4
UKBDIR=$5

cd $WKDIR
/well/bag/clme1992/regenie_v3.2.8.gz_x86_64_Centos7_mkl \
  --threads 13 \
  --step 1 \
  --bed ${UKBDIR}/ukb53100_s488264.imp_merged_QC \
  --phenoFile ${WDIR}/ukb41482.bd.gwas-pheno.${STEM}.${STRAT}.txt.gz \
  --covarFile ${WDIR}/ukb41482.bd.gwas-covar.${STEM}.${STRAT}.txt.gz \
  --bt \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix ${WDIR}/regenie_tmp_preds \
  --out ${WDIR}/step1_${STEM}_${STRAT}