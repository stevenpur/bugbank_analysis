
# For now, the analysis will consider only the phenotype of EUR ancestry

################
# Data logging #
################
# Every object in lg will be saved
lg <- list()

######################
# Congiruation file  #
######################
source("~/.saige_pipe.config")

##########################
# Command-line arguments #
##########################
help <- paste(
    "Usage: Rscript hgi-stratify.R stem",
    sep = "\n"
)
# Argument
lg$args <- commandArgs(trailingOnly = TRUE)
print("Arguments:")
print(lg$args)
if (length(lg$args) != 1) {
    cat(help, sep = "\n")
    stop("\nIncorrect usage\n")
}
# outfile prefix
lg$stem <- as.character(lg$args[1])
if (is.na(lg$stem) | lg$stem == "") stop("Nonsensical stem")

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
# Joint phenotype file
lg$joint_phenotypes_file <- paste0(config$wrkdir, "/ukb41482.bd.gwasdata.", lg$stem, ".txt")

# Output files
lg$log_rds_outfile <- paste0(config$wrkdir, "/log.ukb41482.gwas-stratify.", lg$stem, ".rds")
# phenotype: UK Biobank data field reference
lg$data_outfileprefix <- paste0(config$wrkdir, "/ukb41482.bd.gwas-stratify.", lg$stem)
lg$data_outfiles <- list()
tryCatch(
    {
        #############################
        # Load the joint phenotypes #
        #############################
        data <- read.csv(lg$joint_phenotypes_file, sep = "\t")

        pheno_cols <- setdiff(colnames(data), c("eid", "sex", "age", "ancestry"))
        print("Found phenotypes called:")
        print((lg$pheno_colnames <- pheno_cols))

        print("Post all filters pre stratification: phenotype counts")
        print((lg$table_filtered_pheno <- apply(as.matrix(data[, pheno_cols]), 2, table, useNA = "a")))

        print("NA counts per column in data:")
        print((lg$na_counts_columns_data <- colSums(is.na(data))))


        #########################################
        # Write quorate sub-populations to file #
        #########################################
        lg$n.cases <- list()
        for (pheno in pheno_cols) {
            stratum <- pheno
            lg$data_outfiles[[pheno]] <- paste0(lg$data_outfileprefix, ".", stratum, ".txt.gz")
            out <- data.frame("eid" = data$eid, "pheno" = data[, pheno], "age" = data[, "age"], "sex" = data[, "sex"])
            # WARNING: hacky approach here, need to update for long term usage
            # out$pheno[which(is.na(data$ancestry))] <- NA
            # out$pheno[which(data$ancestry != "eur")] <- NA
            lg$n.cases[[stratum]] <- sum(out$pheno == 1, na.rm = T)
            print(paste0("number of cases for ", stratum, ":", lg$n.cases[[stratum]]))
            write.table(out, file = gzfile(lg$data_outfiles[[pheno]]), row = FALSE, col = TRUE, quote = FALSE, sep = "\t")
        }
    },
    finally = {
        # On error or clean exit
        saveRDS(lg, file = lg$log_rds_outfile)
    }
)
