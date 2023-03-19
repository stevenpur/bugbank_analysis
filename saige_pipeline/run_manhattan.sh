#!/bin/bash
#SBATCH -D ./
#SBATCH -A bag.prj
#SBATCH -p short
#SBATCH -J hgi-manhattan
#SBATCH --cpus-per-task=7

TMPDIR=/well/bag/clme1992/tmp
export TMPDIR

echo "SGE job ID: "$SLURM_JOB_ID
echo "SGE task ID: "$SLURM_ARRAY_TASK_ID
echo "Running on host: "$SLURM_JOB_NODELIST
echo "Output file: "$SLURM_JOB_NAME".o"$SLURM_JOB_ID"."$SLURM_ARRAY_TASK_ID
echo "Script arguments: "
echo "stem:    "$1
echo "stratum: "$2

# Limit multi-threading for the SGE environment
export OMP_NUM_THREADS=$NSLOTS

# Arguments
STEM=$1
STRATUM=$2
WDIR=$3
SRCDIR=$4

echo $1
echo $2
echo $3
echo $4
#echo $SRCDIR
# Working directory
#cd $WDIR

# Run command
module purge
module add R/3.6.2-foss-2019b
Rscript ${SRCDIR}/saige_pipeline/plot_manhattan.R $1 $2 #highlight
