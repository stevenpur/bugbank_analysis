# this code is to decide which pathogens should be feed into the SAIGE pipeline

library(data.tree)
library(tidyverse)
setwd("~/bugbank_data")

# input files
hes_diag_f <- "/well/bag/wilson/ukb/hes/hesin_diag.latest.txt.gz"
hes_tax_tree_f <- "./hes_pathogen_tax_tree_refined_13032023.rds"
ukb_f <- "/well/bag/wilson/ukb/ukb41482.ukb41376.fields.RData"

# load files
hes_diag <- read.csv(hes_diag_f, sep = "\t")
hes_tax_tree <- readRDS(hes_tax_tree_f)
load(ukb_f)
ukb_eids <- bd[, "f.eid"]

# for computational efficiency, I pre-compute the eids for each ICD10 code
all_icd10 <- unique(unlist(hes_tax_tree$Get(function(node) node$cum_icd10)))
hes_icd10_eid_lgl <- map(all_icd10, function(icd10) {
    cur_eids <- unique(hes_diag$eid[which(hes_diag$diag_icd10 == icd10)])
    # for computational efficiency down stream, store eid info into a boolin vector
    # the vector has the same order as ukb_eid
    # Note: this means that eid in hes_diag but not in ukb_f will be omitted
    return(!is.na(cur_eids[match(ukb_eids, cur_eids)]))
})
names(hes_icd10_eid_lgl) <- all_icd10

# assign eids boolin vector to each node
hes_tax_tree$Do(function(node) {
    cur_icd10s <- node$cum_icd10
    eid_lgl <- rep(FALSE, length(ukb_eids))
    for (icd10 in cur_icd10s) {
        eid_lgl <- eid_lgl | hes_icd10_eid_lgl[[icd10]]
    }
    node$cum_eid_lgl <- eid_lgl
})

# traverse the tree, if a node has
# 1. more than 110% case eids than any of its children nodes
# 2. more than 100 cases eids
# then we can include the node into analysis
in_analysis <- hes_tax_tree$Get(function(node) {
    if (sum(node$cum_eid_lgl) < 100) {
        return(FALSE)
    }
    # check if the eids of any child is more than 90% of the current node
    # if yes, then the current node should not be included in the analysis
    children_eids <- vector()
    for (child in node$children) {
        if (sum(node$cum_eid_lgl) * 0.95 < sum(child$cum_eid_lgl)) {
            return(FALSE)
        }
    }
    return(TRUE)
})

in_analysis <- names(in_analysis)[unlist(in_analysis)]

# assign eids to each node
# hes_tax_tree$Do(function(node) {
#    cur_icd10 <- node$cum_icd10
#    sub_hes_diag <- hes_diag %>% filter(diag_icd10 %in% cur_icd10)
#    node$cum_eids <- unique(sub_hes_diag$eid)
# })

# traverse the tree, if a node has
# 1. more than 110% case eids than any of its children nodes
# 2. more than 100 cases eids
# then we can include the node into analysis
in_analysis <- hes_tax_tree$Get(function(node) {
    if (length(node$cum_eids) < 100) {
        return(FALSE)
    }
    # for each child node, check if the eids of the child is more than 90% of the current node
    # if yes, then the current node should not be included in the analysis
    children_eids <- vector()
    for (child in node$children) {
        if (length(node$cum_eids) * 0.95 < length(child$cum_eids)) {
            return(FALSE)
        }
    }
    return(TRUE)
})

in_analysis <- names(in_analysis)[unlist(in_analysis)]

# write result to files for each taxa level
tax_class <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus")
for (tax_lev in tax_class) {
    message(tax_lev)
    cur_analysis <- in_analysis[grep(tax_lev, in_analysis)]
    cur_analysis <- strsplit(cur_analysis, ":")
    result <- do.call(rbind, cur_analysis)
    colnames(result) <- c("tax_lev", "name")
    wrt_file <- paste0("13032023_hes_", tax_lev, ".input.tsv")
    write.table(result, sep = "\t", col.names = T, row.names = F, quote = F)
}
