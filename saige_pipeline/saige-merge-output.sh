#!/bin/bash
#SBATCH -D ./
#SBATCH -A bag.prj
#SBATCH -p short
#SBATCH -J merge-saige
#SBATCH --cpus-per-task=1

#TMPDIR=/well/bag/clme1992/tmp
#export TMPDIR

echo "SGE job ID: "$SLURM_JOB_ID
echo "SGE task ID: "$SLURM_ARRAY_TASK_ID
echo "Running on host: "$SLURM_JOB_NODELIST
echo "Output file: "$SLURM_JOB_NAME".o"$SLURM_JOB_ID"."$SLURM_ARRAY_TASK_ID
echo "Script argument: "$1

# Stem (read from argument)
STEM=$1
# Output file
OUTFILE=merged.${STEM}.txt.gz
# Working directory
WDIR=$2/saige
# Trash dir
TRASH=$2/trash/${STEM}/

# Merge all zipped SAIGE imputation result files 1:1378
cd $WDIR
cat ${STEM}.batch1.txt.gz > $OUTFILE
mv ${STEM}.batch1.txt.gz $TRASH
# Loop over subsequent files
for i in {2..1325}
do
	cat ${STEM}.batch${i}.txt.gz >> $OUTFILE
	mv ${STEM}.batch${i}.txt.gz $TRASH
done
