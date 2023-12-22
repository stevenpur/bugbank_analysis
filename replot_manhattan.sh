# This code is run when the regenie pipeline is finished, 
# but the user made adjustments to the plotting details of the manhattan plots afterwards 
# and wish to change the plots, but not the association analysis itself

# read the parameters
STEM=$1
srcdir=$HOME/bugbank/github2/regenie_pipeline

# go the the result diretory from home
cd ${HOME}/saige_pipe_test/
# get the list of summary files that need to be replotted
SUMMARY_FILES=$(ls summary.${STEM}.*.txt.gz)
# loop through the summary files and plot the result
for SUMMARY_FILE in ${SUMMARY_FILES}
do
    # get the strat name from the filename
    # replace head summary\..*\. to nothing
    # replace tail .txt.gz to nothing
    STRAT=$(echo ${SUMMARY_FILE} | sed "s/summary\.${STEM}\.//g" | sed 's/\.txt\.gz//g')
    # plot the manhattan plot
    echo "sbatch $HOME/bugbank/github2/regenie_pipeline/plot_manhattan.R $STEM $STRAT $srcdir"
done
