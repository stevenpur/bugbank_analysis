# screen -x 142818
# export OMP_NUM_THREADS=10

#######################
# Sample-specific PCA #
#######################
# 19 August 2020
library(bigsnpr)

# Every object in lg will be saved
lg <- list()

# Configuration file #
source("~/.saige_pipe.config")
setwd(config$wrkdir)

# Command-line arguments
help <- paste(
    "Usage: Rscript gwas-bigsnpr-pca.R stem stratum npcs ncores",
    sep = "\n"
)

# Argument
lg$args <- commandArgs(trailingOnly = TRUE)
print("Arguments:")
print(lg$args)
if ((length(lg$args) != 4)) {
    cat(help, sep = "\n")
    stop("\nIncorrect usage\n")
}

# outfile prefix
lg$stem <- as.character(lg$args[1])
if (is.na(lg$stem) || lg$stem == "") stop("Nonsensical stem")

# stratification group
lg$stratum <- as.character(lg$args[2])
if (is.na(lg$stratum) || lg$stratum == "") stop("Nonsensical stratum")

# number of PCs to calculate
lg$npcs <- as.numeric(lg$args[3])
if (is.na(lg$npcs) || lg$npcs < 1) stop("Nonsensical npcs")
if (lg$npcs > 50) stop("Not suitable for large numbers of PCs")

# number of cores to use
lg$ncores <- as.numeric(lg$args[4])
if (is.na(lg$ncores) || lg$ncores < 1) stop("Nonsensical ncores")
if (lg$ncores > 40) stop("Not suitable for large ncores")

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
# Created by hgi-bigsnpr-pca-prepare.R
lg$backingfile <- paste0(
    config$ukb.derived.dir,
    "/ukb53100_s488264.imp_merged_QC.bigsnpr.backingfile.rds"
)
# Created by hgi-stratify.R
lg$data_infile <- paste0("ukb41482.bd.gwas-stratify.", lg$stem, ".", lg$stratum, ".txt.gz")
# Output file
lg$data_covar_outfile <- paste0("ukb41482.bd.gwas-covar.", lg$stem, ".", lg$stratum, ".txt.gz")
lg$data_pheno_outfile <- paste0("ukb41482.bd.gwas-pheno.", lg$stem, ".", lg$stratum, ".txt.gz")

print(lg$data_infile)
print(lg$data_outfile)

#################
# Load the data #
#################
# Attach the plink file
plink <- snp_attach(lg$backingfile)

# Load the phenotype data
data <- read.delim(lg$data_infile, as.is = TRUE, stringsAsFactors = FALSE)

###################
# Compute the PCs #
###################
# Check no problems with overlap of eids for non-NA phenotypes
data2plink <- match(data$eid, plink$fam$sample.ID)
sum(!is.na(data2plink)) # 487295
sum(is.na(data2plink)) #  15210
# Issue a warning if any non-NA phenotypes have NA data2plink
stopifnot(!any(is.na(data2plink[!is.na(data$pheno)])))

# Define the eids for analysis
eid_include <- data$eid[apply(!is.na(data), 1, all)]
length(eid_include)
bed_include <- which(!is.na(match(plink$fam$sample.ID, eid_include)))
length(bed_include)

# Compute sample frequencies
# ct is a 4xL matrix of genotype counts in the order 0, 1, 2, NA
# where L is the number of SNPs
system.time((ct <- big_counts(plink$genotypes, ind.row = bed_include)))
#   user  system elapsed
# 29.007  34.244  63.385

# Sanity check that there are no missing genotype calls
any_miss <- (ct[4, ] > 0)
stopifnot(sum(any_miss) == 0)

# Identify any variants fixed in the sample
any_fixed <- (ct[1, ] == 0 & ct[2, ] == 0) |
    (ct[1, ] == 0 & ct[3, ] == 0) |
    (ct[2, ] == 0 & ct[3, ] == 0) |
    (ct[1, ] == 0 & ct[2, ] == 0 & ct[3, ] == 0)
print("Number of variants fixed in the sample:")
print(sum(any_fixed))

# Identify variants for analysis (all not missing and not fixed)
ind_col <- setdiff(1:ncol(plink$genotypes), which(any_miss | any_fixed)) # nolint: seq_linter.
print("Number of variants for analysis:")
print(length(ind_col))

# For replicability, set random seed to zero
set.seed(0)

# Assert no parallel BLAS is running (did not run in testing)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

# Run the 'fast' SVD
print("running bigsnpr ncores:")
print(lg$ncores)

system.time((
    svd <- big_randomSVD(
        plink$genotypes,
        snp_scaleBinom(),
        ind.row = bed_include,
        ind.col = ind_col,
        verbose = FALSE,
        ncores = lg$ncores,
        k = lg$npcs
    )
))
# For 176k non-NA phenotypes on 10 cores (in testing):
#  user  system elapsed
# 0.232  16.841 683.671

# Project individuals on to the PCs
system.time((pcs <- predict(svd)))
#  user  system elapsed
# 0.026   0.009   0.035

# extract the phenotype files
pheno <- data[, c("eid", "eid", "pheno")]
colnames(pheno) <- c("FID", "IID", "pheno")

# Incorporate the PCs into the covariate files
data_pcs <- as.data.frame(
    matrix(NA, nrow(data),
        lg$npcs,
        dimnames = list(NULL, paste0("pc", 1:lg$npcs))
    )
)
data_pcs[match(plink$fam$sample.ID[bed_include], data$eid), ] <- pcs
data <- cbind(data, data_pcs)
data <- data[, c("eid", "eid", "age", "sex", paste0("pc", 1:lg$npcs))]
colnames(data) <- c("FID", "IID", "age", "sex", paste0("pc", 1:lg$npcs))
data$sex <- map_int(data$sex, function(x) {
    if (x == "Male") {
        return(1)
    }
    if (x == "Female") {
        return(0)
    }
    return(NA)
})

# Write the files
write.table(
    pheno,
    file = gzfile(lg$data_pheno_outfile),
    row = FALSE,
    col = TRUE,
    quote = FALSE,
    sep = "\t"
)
write.table(
    data,
    file = gzfile(lg$data_covar_outfile),
    row = FALSE,
    col = TRUE,
    quote = FALSE,
    sep = "\t"
)
