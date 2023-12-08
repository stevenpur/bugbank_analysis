#!/bin/bash
#SBATCH -D ./
#SBATCH -A bag.prj
#SBATCH -p short
#SBATCH -J hgi-regenie-gwas
#SBATCH --cpus-per-task=1

echo "SGE job ID: "$SLURM_JOB_ID
echo "SGE task ID: "$SLURM_ARRAY_TASK_ID
echo "Running on host: "$SLURM_JOB_NODELIST
echo "Output file: "$SLURM_JOB_NAME".o"$SLURM_JOB_ID"."$SLURM_ARRAY_TASK_ID
echo "Script arguments: "
echo "stem:    "$1
echo "stratum: "$2
echo "ncores:  "$3
echo "npcs: "$4
echo "working dir: "$5

# Limit multi-threading for the SGE environment
export OMP_NUM_THREADS=$NSLOTS

# Arguments
STEM=$1
STRATUM=$2
NCORES=$3
NPCS=$4
WDIR=$5
SRCDIR=$6

# Working directory
cd $WDIR

# Run command
module purge
module add R/3.6.2-foss-2019b
Rscript ${SRCDIR}/run_pipeline_sub.R $1 $2 $3 $4

