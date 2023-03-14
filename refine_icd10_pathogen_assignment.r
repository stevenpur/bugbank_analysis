# In the raw pathogen assignment file, there are instances where a single ICD10 code is assinged to multiple pathogens
# This is mainly due to the same ICD10 code being assign to a species and then also to its genus
# This should be avoided so this code is used to adjust for this problem

# input file
path_desc_raw_f <- "./icd10_pathogen_description_13032023.tsv"
hes_taxtree_raw_f <- "./hes_pathogen_tax_tree_fromRaw_13032023.rds"

# load file
path_desc_raw <- read.csv(path_desc_raw_f, sep = "\t")
hes_taxtree_raw <- readRDS(hes_taxtree_raw_f)

# build a dictionary to map pathgoen to icd10
path_to_icd10 <- list()
for (i in 1:nrow(path_desc_raw)) {
    pathogen <- path_desc_raw$org_name[i]
    icd10s <- strsplit(path_desc_raw$UKB_code[i], ",")[[1]]

    # for path2icd10
    path_to_icd10[[pathogen]] <- unique(c(path_to_icd10[[pathogen]], icd10s))
}

# for each node of the species tree, remove icd10 that is presented in its children node
# use callback function to achieve this
remove_redundant_icd10 <- function(node) {
    # create a stash to store all of the ICD10s in the children nodes
    children_icd10_stash <- vector()

    if (node$isLeaf) {
        # get the icd10 code
        if (! grepl("origin_name:", node$name)) {
            stop(paste0("leaf node is not origin_name for ", node$name))
        }
        pathogen_name <- gsub("origin_name:", "", node$name)
        return(path_to_icd10[[pathogen_name]])
    }

    if (!node$isLeaf) {
        # get icd10 from all its children
        for (child in node$children) {
            children_icd10_stash <- c(
                children_icd10_stash,
                remove_redundant_icd10(child)
            )
        }

    }
}





# using the dictionary, add a column in the taxonomy table to include ICD10 info
hes_tax_raw$icd10 <- map(hes_tax_raw$origin_name, function(org_name) {
    path_to_icd10[[org_name]]
})

# starting from the species level, purge all duplicated icd10 codes from all parent taxonomy
for (tax_col in 8:1) {
    for ()
}