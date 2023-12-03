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
CHR=1
cat step2.${STEM}_${STRAT}_chr${CHR}.gz > summary.$STEM.$STRAT.txt.gz
for i in {2..24}
do
    echo "chromosome"$CHR
    CHR=${i}
    if [ "$CHR" = "23" ]; then
        CHR="X"
    fi

    if [ "$CHR" = "24" ]; then
        CHR="XY"
    fi
        
    cat step2.${STEM}_${STRAT}_chr${CHR}.gz >> summary.$STEM.$STRAT.txt.gz
done