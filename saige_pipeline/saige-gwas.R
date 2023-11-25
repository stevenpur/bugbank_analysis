safari ##################
# Run SAIGE GWAS #
##################
# 22 August 2020
library(batchtools)

################
# Data logging #
################
# Every object in lg will be saved
lg <- list()

######################
# Configuration file #
######################
source("~/.saige_pipe.config")

##########################
# Command-line arguments #
##########################
help <- paste(
    "Usage: Rscript hgi-saige-gwas.R stem stratum ncores npcs",
    sep = "\n"
)
# Argument
lg$args <- commandArgs(trailingOnly = TRUE)
print("Arguments:")
print(lg$args)
if (length(lg$args) != 4) {
    cat(help, sep = "\n")
    stop("\nIncorrect usage\n")
}
# outfile prefix
lg$stem <- as.character(lg$args[1])
if (is.na(lg$stem) || lg$stem == "") stop("Nonsensical stem")
# stratification group
lg$stratum <- as.character(lg$args[2])
if (is.na(lg$stratum) || lg$stratum == "") stop("Nonsensical stratum")
# number of cores to use
lg$ncores <- as.numeric(lg$args[3])
if (is.na(lg$ncores) || lg$ncores < 1) stop("Nonsensical ncores")
if (lg$ncores > 40) stop("Not suitable for large ncores")
# number of PCs to use
lg$npcs <- as.numeric(lg$args[4])
if (is.na(lg$npcs) || lg$npcs < 1) stop("Nonsensical npcs")
if (lg$npcs > 50) stop("Not suitable for large numbers of PCs")

################
# Code version #
################
# Record GitLab version
print("User:")
print((lg$username <- Sys.info()["user"]))
print("Source directory:")
print((lg$srcdir <- config$srcdir))
print("Git repository version")
system(paste0("(cd ", lg$srcdir, " && git show --oneline -s)"))

##########################
# Input and output files #
##########################
setwd(config$wrkdir)
# Input file created by hgi-bigsnpr-pca.R
lg$infilename.gz <- paste0("ukb41482.bd.gwas-bigsnpr-pca.", lg$stem, ".", lg$stratum, ".txt.gz")
lg$infilename <- paste0("ukb41482.bd.gwas-bigsnpr-pca.", lg$stem, ".", lg$stratum, ".txt")
# Step 1 output files
lg$stem.out <- paste0(lg$stem, ".", lg$stratum)
lg$outputPrefix <- paste0(config$wrkdir, "/saige/step1.", lg$stem.out)
lg$log.rds.outfilename <- paste0(
    config$wrkdir,
    "/saige/log.ukb41482.gwas-saige-gwas.", lg$stem.out, ".rds"
)

# Enclose what follows in a tryCatch statement so the log is output
# even if R quits with an error
tryCatch(
    {
        ##################
        # PCA by bigsnpr #
        ##################
        jobname <- paste0("bigsnpr_", lg$stem, "_", lg$stratum)
        lg$bigsnpr.cmd <- paste(
            "sbatch",
            "-o", paste0(config$bigsnpr.stddir, "/", jobname, ".out"),
            "-J", jobname,
            paste0("--cpus-per-task=", 3 * (lg$ncores + 1)),
            paste0(config$srcdir, "/gwas_run_pca.sh"),
            lg$stem, lg$stratum, lg$npcs, lg$ncores, config$wrkdir, config$srcdir
        )

        # Submit and record jobid
        lg$bigsnpr.jobid <- readLines(pipe(lg$bigsnpr.cmd))
        print("bigsnpr jobid:")
        print(lg$bigsnpr.jobid)
        lg$bigsnpr.jobid <- unlist(
            strsplit(unlist(strsplit(lg$bigsnpr.jobid, " "))[4], ".", fixed = TRUE)
        )[1]
        print(lg$bigsnpr.jobid)

        ####################
        # Run SAIGE step 1 #
        ####################
        jobname <- paste0("saige-step1_", lg$stem, ".", lg$stratum)
        lg$step1.cmd <- paste(
            "sbatch",
            paste0("--dependency=afterok:", lg$bigsnpr.jobid),
            "-o", paste0(config$saige.stddir, "/", jobname),
            "-J", jobname,
            paste0("--cpus-per-task=", (lg$ncores * 3)),
            paste0(config$srcdir, "/saige-step1-runner.sh"),
            # shell script arguments to pass
            lg$stem, lg$stratum, config$wrkdir,
            lg$ncores * 3, config$ukb.derived.dir, config$ukbdir,
            config$saige.container.dir, lg$outputPrefix
        )

        # Submit and record jobid
        lg$step1.jobid <- readLines(pipe(lg$step1.cmd))
        print("Step 1 jobid:")
        print(lg$step1.jobid)
        lg$step1.jobid <- unlist(
            strsplit(unlist(strsplit(lg$step1.jobid, " "))[4], ".", fixed = TRUE)
        )[1]

        ####################
        # Run SAIGE Step 2 #
        ####################
        jobname <- paste0("saige-step2_", lg$stem, ".", lg$stratum)
        lg$step2.cmd <- paste(
            "sbatch",
            paste0("--dependency=afterok:", lg$step1.jobid),
            "-o", paste0(config$saige.stddir, "/", jobname, "_%a"),
            "-J", jobname,
            paste0(config$srcdir, "/saige-step2-runner.sh"),
            # shell script arguments to pass
            lg$stem.out, config$wrkdir, config$ukb.derived.dir,
            config$ukbdir, config$saige.container.dir
        )

        # Submit and record jobid
        lg$step2.jobid <- readLines(pipe(lg$step2.cmd))
        print("Step 2 jobid:")
        print(lg$step2.jobid)
        lg$step2.jobid <- unlist(
            strsplit(unlist(strsplit(lg$step2.jobid, " "))[4], ".", fixed = TRUE)
        )[1]

        #################################
        # Merge the SAIGE summary files #
        #################################
        jobname <- paste0("regenie-merge_", lg$stem, ".", lg$stratum)
        lg$step3.cmd <- paste(
            "sbatch",
            paste0("--dependency=afterok:", lg$step2.jobid),
            "-o", paste0(config$saige.stddir, "/", jobname),
            "-J", jobname,
            paste0(config$srcdir, "/merge_files.sh"),
            # shell script arguments to pass
            lg$stem, lg$stratum, config$wrkdir
        )

        # Submit and record jobid
        lg$step3.jobid <- readLines(pipe(lg$step3.cmd))
        print("Step 3 jobid:")
        print(lg$step3.jobid)
        lg$step3.jobid <- unlist(strsplit(lg$step3.jobid, " "))[4]

        ##############################
        # Create the Manhattan plots #
        ##############################
        jobname <- paste0("regenie-manhattan_", lg$stem, "_", lg$stratum)
        lg$step5.cmd <- paste(
            "sbatch",
            paste0("--dependency=afterok:", lg$step4.jobid),
            "-J", jobname,
            "-o", paste0(config$manhattan.stddir, "/", jobname),
            paste0(config$srcdir, "/run_manhattan.sh"),
            # shell script arguments to pass
            lg$stem, lg$stratum, config$wrkdir, config$srcdir
        )

        # Submit and record jobid
        lg$step5.jobid <- readLines(pipe(lg$step5.cmd))
        print("Step 5 jobid:")
        print(lg$step5.jobid)
        lg$step5.jobid <- unlist(strsplit(lg$step5.jobid, " "))[4]

        ###########################
        # Clean intermediate data #
        ###########################
        jobname <- paste0("cleaning_", lg$stem, "_", lg$stratum)
        lg$step6.cmd <- paste(
            "sbatch",
            paste0("--dependency=afterok:", lg$step5.jobid),
            "-J", jobname,
            "-o", paste0(config$manhattan.stddir, "/", jobname),
            paste0(config$srcdir, "/remove_intermediate_files_strat.sh"),
            lg$stem, lg$stratum
        )
        lg$step6.jobid <- readLines(pipe(lg$step6.cmd))
        print("Step 6 jobid:")
        print(lg$step6.jobid)
        lg$step6.jobid <- unlist(strsplit(lg$step6.jobid, " "))[4]
    },
    finally = {
        # On error or clean exit
        saveRDS(lg, file = lg$log.rds.outfilename)
    }
)
