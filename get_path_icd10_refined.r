# This is the code to use the refined taxonomy tree to create a table mapping pathogen to its icd10
# unlike the raw file icd10_pathogen_description_13032023.tsv, this new map should not have the problem of a single ICD10 mapping to multiple organism

# Note: the result of this new mapping would result in some origin_name disappearing, the reason is that some of the origin_names are more specific than the species level
# in such instances, the origin_names would be collapse to the species levels

library(data.tree)
setwd("~/bugbank_data")
# input file
hes_tax_tree_f <- "./hes_pathogen_tax_tree_refined_13032023.rds"
# load file
hes_tax_tree <- readRDS(hes_tax_tree_f)

path_to_uni_icd10 <- list()

# travel the tree nodes, assign each uni_icd10 to the name of the node
path_to_uni_icd10 <- hes_tax_tree$Get(function(node) {
    paste0(node$uni_icd10, collapse = ",")
})
path_to_uni_icd10 <- path_to_uni_icd10[-which(path_to_uni_icd10 == "")]

path_name <- map_chr(strsplit(names(path_to_uni_icd10), ":"), 2)
path_tax_lev <- map_chr(strsplit(names(path_to_uni_icd10), ":"), 1)
result <- data.frame(
    tax_lev = path_tax_lev,
    org_name = path_name,
    icd10 = unlist(path_to_uni_icd10)
)

write.table(result, "./pathogen_to_unique_icd10.tsv", sep = "\t", col.names = T, row.names = F, quote = F)