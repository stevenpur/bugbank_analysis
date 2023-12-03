#!/bin/bash
#SBATCH -D ./
#SBATCH -A bag.prj
#SBATCH -p short
#SBATCH -J remove_files
#SBATCH --cpus-per-task=1

# clean the files from saige
stem=$1
strat=$2
result_dir=/well/bag/clme1992/saige_pipe_test/
# check if the Manhattan files that have been produced
stratum=${stem}.${strat}
manhattan_f=Manhattan.${stem}.${strat}_regenie.png
if [ ! -f ${result_dir}/$manhattan_f]; then
    echo "manhattan file does not exist for $stem.$stratum"
    exit 1
fi
# look to see if summary file exist
if [ ! -f ${result_dir}/summary.${stratum}.txt.gz ]; then
    echo "merged results does not exist for: ${stem}.${stratum}"
    exit 1
fi

f_size=$(ls -l ${result_dir}/summary.${stratum}.txt.gz | awk '{print $5}')
echo $stratum
if [ $f_size -lt 4000000000 ] && [ $f_size -gt 1000000000 ]; then
    echo "removing bigsnpr files..."
    rm -f /well/bag/clme1992/saige_pipe_test/ukb41482.bd.hgi-bigsnpr-pca.${stratum}.txt
    echo "removing step1 files..."
    rm -f /well/bag/clme1992/saige_pipe_test/step1_${stem}_${strat}*
    echo "removing unmerged files"
    rm -f /well/bag/clme1992/saige_pipe_test/step2.${stem}_${strat}*
else
    echo "suspicious merge files size for ${stratum}"
    exit 1 
fi