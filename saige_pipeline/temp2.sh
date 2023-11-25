#!/bin/bash
#SBATCH -D ./
#SBATCH -A bag.prj
#SBATCH -p short
#SBATCH -J step2-saige
#SBATCH --cpus-per-task=4

#TMPDIR=/well/bag/clme1992/tmp
#export TMPDIR
echo "SGE job ID: "$SLURM_JOB_ID
echo "SGE task ID: 66"
echo "Running on host: "$SLURM_JOB_NODELIST
echo "Output file: "$SLURM_JOB_NAME".o"$SLURM_JOB_ID".66"
echo "Script argument: "$1 $2 $3 $5

# Limit multi-threading for the SGE environment
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Stem (read from argument)
STEM=$1
# Working directory
WDIR=$2/saige
# UKB derived data directory
UKBDEDIR=$3
# UKB data directory
UKBDIR=$4
# SAIGE container directory
SAIGEDIR=$5
# Path to the file that contains one column for IDs of samples in the dosage file
SAMP=${UKBDEDIR}/ukb41482.imp.bgen.eid.txt
# Model file from step 1
GMMAT=${WDIR}/step1.${STEM}.rda
# Variance ratio file from step 1
VARRAT=${WDIR}/step1.${STEM}.varianceRatio.txt
# Output files
OUT=${WDIR}/${STEM}.batch66.temp.txt
OUTGZ=${OUT}.gz
OUTSUM=${OUT}.sum.gz
TMP=/well/bag/clme1992/tmp

# Batch file specifies variants for analysis
BATCH=${UKBDEDIR}/batch/batch66.v3.txt
# Chromosome to analyse
SCHR=`cut -f1 $BATCH`
CHR=`expr ${SCHR} + 0`
XCHR="X"
PAR1="PAR1"
if [ "$SCHR" = "$XCHR" ]; then
	CHR="X"
	SAMP=${UKBDEDIR}/ukb41482.imp-X.bgen.eid.txt
fi
if [ "$SCHR" = "$PAR1" ]; then
	CHR="XY"
	SAMP=${UKBDEDIR}/ukb41482.imp-XY.bgen.eid.txt
fi
echo "Chromosome names, internal ${SCHR} and filename ${CHR}"
# Genotype data input file
BGEN=${UKBDIR}/v3/imputation/ukb_imp_chr${CHR}_v3.bgen

# Run SAIGE
# -B ${TMP} --env TMPDIR=${TMP} ${SAIGEDIR}
cd $WDIR
if [ "$SCHR" = "$PAR1" ]; then
	time singularity run -B ${UKBDEDIR} -B ${UKBDIR} -B ${WDIR} step2_SPAtests.R  --bgenFile=$BGEN  --bgenFileIndex=${BGEN}.bgi  --minMAF=0  --minMAC=0  --sampleFile=$SAMP  --GMMATmodelFile=$GMMAT  --varianceRatioFile=$VARRAT  --SAIGEOutputFile=$OUT  --numLinesOutput=100  --IsOutputAFinCaseCtrl=TRUE --chrom=${CHR} --LOCO=FALSE
elif [ "$SCHR" = "$XCHR" ]; then
	time singularity run -B ${UKBDEDIR} -B ${UKBDIR} -B ${WDIR}  ${SAIGEDIR} step2_SPAtests.R  --bgenFile=$BGEN  --bgenFileIndex=${BGEN}.bgi  --rangestoIncludeFile=$BATCH  --minMAF=0  --minMAC=0  --sampleFile=$SAMP  --GMMATmodelFile=$GMMAT  --varianceRatioFile=$VARRAT  --SAIGEOutputFile=$OUT  --numLinesOutput=100  --IsOutputAFinCaseCtrl=TRUE --chrom=${CHR} --LOCO=FALSE
else
	time singularity run -B ${UKBDEDIR} -B ${UKBDIR} -B ${WDIR}  ${SAIGEDIR} step2_SPAtests.R  --bgenFile=$BGEN  --bgenFileIndex=${BGEN}.bgi  --rangestoIncludeFile=$BATCH  --minMAF=0  --minMAC=0  --sampleFile=$SAMP  --GMMATmodelFile=$GMMAT  --varianceRatioFile=$VARRAT  --SAIGEOutputFile=$OUT  --numLinesOutput=100  --IsOutputAFinCaseCtrl=TRUE --chrom=${CHR}
fi

# Post-process
ONE="1"
if [ "66" = "$ONE" ]; then
	cut -f1,2,3,7,13 -d' ' $OUT | tail -n +1 | gzip > $OUTSUM
	cat $OUT | tail -n +1 | gzip > $OUTGZ
else
	cut -f1,2,3,7,13 -d' ' $OUT | tail -n +2 | gzip > $OUTSUM
	cat $OUT | tail -n +2 | gzip > $OUTGZ
fi
rm -f $OUT
