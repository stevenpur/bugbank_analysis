# screen -x 154198
# R --args temp 10
############
# Run GWAS #
############
# 22 August 2020

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
    "Usage: Rscript hgi-run-gwas.R stem ncores npcs",
    sep = "\n"
)
# Argument
lg$args <- commandArgs(trailingOnly = TRUE)
print("Arguments:")
print(lg$args)
if (length(lg$args) != 3) {
    cat(help, sep = "\n")
    stop("\nIncorrect usage\n")
}
# outfile prefix
lg$stem <- as.character(lg$args[1])
if (is.na(lg$stem) | lg$stem == "") stop("Nonsensical stem")
# number of cores to use
lg$ncores <- as.numeric(lg$args[2])
if (is.na(lg$ncores) | lg$ncores < 1) stop("Nonsensical ncores")
if (lg$ncores > 40) stop("Not suitable for large ncores")
# number of PCs to use
lg$npcs <- as.numeric(lg$args[3])
if (is.na(lg$npcs) | lg$npcs < 1) stop("Nonsensical npcs")
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

# Read hgi-stratify file
print(config$wrkdir)
strat_logfile <- paste0("./log.ukb41482.hgi-stratify.", lg$stem, ".rds")
strat <- readRDS(strat_logfile)
lg$strat_names <- names(strat$n.cases[strat$n.cases >= 100])
print("Not running the following strata due to small case counts (<100):")
print(strat$n.cases[strat$n.cases < 100])

# Log file
lg$log.outfile <- paste0("log.ukb41482.bd.hgi-run-gwas.", lg$stem, ".rds")

################
# Run the GWAS #
################
lg$cmds <- paste(
    "env -i SGE_ROOT=/mgmt/uge/8.6.8 /mgmt/uge/8.6.8/bin/lx-amd64/qsub",
    "-o", config$saige.stddir,
    config$srcdir, "/pheno/hgi-saige-gwas.sh",
    lg$stem,
    lg$strat_names,
    lg$ncores,
    config$wrkdir,
    lg$npcs,
    config$srcdir
)
names(lg$cmds) <- lg$strat_names

for (i in 1:length(lg$cmds)) {
    jobid <- readLines(pipe(lg$cmds[i]))
    lg$jobid[[lg$strat_names[i]]] <- unlist(
        strsplit(unlist(strsplit(jobid, " "))[3],
            ".",
            fixed = TRUE
        )
    )[1]
}

saveRDS(lg, file = lg$log.outfile)
