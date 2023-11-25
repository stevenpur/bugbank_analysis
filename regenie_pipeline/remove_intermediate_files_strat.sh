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
manhattan_f=Manhattan.${stem}.${strat}.png
if [ ! -f ${result_dir}/$manhattan_f]; then
    echo "manhattan file does not exist for $stem.$stratum"
    exit 1
fi
# look to see if the size of the merge file is reasonable
if [ ! -f ${result_dir}/merged.${stratum}.txt.gz ]; then
    echo "merged results does not exist for: ${stem}.${stratum}"
    exit 1
fi

f_size=$(ls -l ${result_dir}/merged.${stratum}.txt.gz | awk '{print $5}')
echo $stratum
if [ $f_size -lt 15000000000 ] && [ $f_size -gt 8000000000 ]; then
    # seems reasonable start deleting
    echo "removing trash directory..."
    #rm -rf /well/bag/clme1992/saige_pipe_test/trash/$stratum
    #echo "removing step2 batch files..."
    #find /well/bag/clme1992/saige_pipe_test/saige/ -name "${stratum}.batch*" -exec rm -f {} \;
    echo "removing bigsnpr files..."
    rm -f /well/bag/clme1992/saige_pipe_test/ukb41482.bd.hgi-bigsnpr-pca.${stratum}.txt
    echo "removing step1 files..."
    rm -f /well/bag/clme1992/saige_pipe_test/step1_${stem}_${strat}*
    echo "removing unmerged files"
    rm -f /well/bag/clme1992/saige_pipe_test/step2.${stem}_${strat}*
    #echo "removing removing step2 cout files..."
    #rm -f /well/bag/clme1992/saige_pipe_test/cout/hgi-saige-gwas/hgi-step2_${stratum}.*
    #echo "removing summary files..."
    #rm -f /well/bag/clme1992/saige_pipe_test/saige/summary.${stratum}.txt.gz
    #echo "moving merge files..."
    #data_folder=$(echo $stem | sed 's/bb/sgss/g' | sed 's/11122021_//g')
    #mv /well/bag/clme1992/saige_pipe_test/saige/merged.${stratum}.txt.gz /data/dingo/wilson/lin/${data_folder}_11122021/
else
    echo "suspicious merge files size for ${stratum}"
    exit 1 
fi