#!/bin/bash
#SBATCH -D ./
#SBATCH -A bag.prj
#SBATCH -p short
#SBATCH -J step3-saige
#SBATCH --cpus-per-task=1

#TMPDIR=/well/bag/clme1992/tmp
#export TMPDIR

echo "SGE job ID: "$SLURM_JOB_ID
echo "SGE task ID: "$SLURM_ARRAY_TASK_ID
echo "Running on host: "$SLURM_JOB_NODELIST
echo "Output file: "$SLURM_JOB_NAME".o"$SLURM_JOB_ID"."$SLURM_ARRAY_TASK_ID

echo "Script argument: "$1 $2

# Stem (read from argument)
STEM=$1
# Working directory
WDIR=$2/saige
TRASH=$2/trash/${STEM}/
mkdir $TRASH

# First file
cd $WDIR
cat ${STEM}.batch1.txt.sum.gz > summary.${STEM}.txt.gz
# Loop over subsequent files
for i in {2..1325}
do
	cat ${STEM}.batch${i}.txt.sum.gz >> summary.${STEM}.txt.gz
	mv ${STEM}.batch${i}.txt.sum.gz ${TRASH}
done


