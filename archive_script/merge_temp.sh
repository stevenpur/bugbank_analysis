#!/bin/bash
#SBATCH -D ./
#SBATCH -A bag.prj
#SBATCH -p short
#SBATCH --cpus-per-task=1

STEM=$1
STRAT=$2
WKDIR=$3

echo "stem:$STEM"
echo "strat:$STRAT"
echo "wkdir:$WKDIR"

cd $WKDIR
for i in {1..24}
do
    CHR=${i}
    echo "chromosome"$CHR
    if [ "$CHR" = "23" ]; then
        CHR="X"
    fi

    if [ "$CHR" = "24" ]; then
        CHR="XY"
    fi
        
    zcat step2.${STEM}_${STRAT}_chr${CHR}.gz | tail -n +2  >> summary.$STEM.$STRAT.txt
done