#!/bin/bash
#SBATCH -D ./
#SBATCH -A bag.prj
#SBATCH -p short
#SBATCH -J step1-saige
#SBATCH --cpus-per-task=1


echo "SGE job ID: "$SLURM_JOB_ID
echo "SGE task ID: "$SLURM_ARRAY_TASK_ID
echo "Running on host: "$SLURM_JOB_NODELIST
echo "Output file: "$SLURM_JOB_NAME".o"$SLURM_JOB_ID"."$SLURM_ARRAY_TASK_ID
echo "Script arguments: "
echo "stem:    "$1
echo "stratum: "$2
echo "working dir: "$3
echo "ncores:  "$4

#TMPDIR=/well/bag/clme1992/tmp
#export TMPDIR

# Arguments
STEM=$1
STRATUM=$2
WDIR=$3
NCORES=$4
UKBDERDIR=$5
UKBDIR=$6
SAIGEDIR=$7
OUTPREFIX=$8
TMP=/well/bag/clme1992/tmp

TAB=$( printf "\t" )
INFILE=${WDIR}/ukb41482.bd.gwas-bigsnpr-pca.${STEM}.${STRATUM}.txt
LOGFILE=${WDIR}/saige/step1.${STEM}.${STRATUM}.log
# Working directory
cd $WDIR

# Check that output files do not already exist
# If an error here, need to manually unlink the step1.files
if [ -f ${OUTPREFIX}.rda -o -f ${OUTPREFIX}.varianceRatio.txt -o -f ${OUTPREFIX}_30markers.SAIGE.results.txt ]; then
    echo "error: step1 files already exists!"
    exit 1
fi
#Check pca file exists
if [ ! -f ${INFILE}.gz ]; then 
    echo "error: infile does not exist!"
    exit 1
fi
# Unzip the input file
gunzip ${INFILE}.gz
# Check it worked
if [ ! -f ${INFILE} ]; then
    echo "error: gunzip failed!"
    exit 1
fi


COVAR=$(head -n 1 ${INFILE} | sed -E "s/eid|pheno//g"| sed "s/${TAB}/,/g"| sed -E "s/^,|,,|,$//g")
time singularity run -B ${UKBDERDIR} -B ${UKBDIR} -B ${WDIR} ${SAIGEDIR} step1_fitNULLGLMM.R --phenoFile=${INFILE} --outputPrefix=${OUTPREFIX} --covarColList=${COVAR} --traitType=binary --invNormalize=FALSE --nThreads=${NCORES} --sampleIDColinphenoFile=eid  --phenoCol=pheno --plinkFile=${UKBDERDIR}/ukb53100_s488264.imp_merged --isCovariateTransform=TRUE &> ${LOGFILE}
# Check it worked
if [ ! ${OUTPREFIX}.rda -o ! -f ${OUTPREFIX}.varianceRatio.txt -o ! -f ${OUTPREFIX}_30markers.SAIGE.results.txt ]; then
    echo "step1 result file not produced!"
    exit 1
fi
rm -f ${INFILE}

