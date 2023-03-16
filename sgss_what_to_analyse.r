# this code is to decide which pathogens should be feed into the SAIGE pipeline

library(data.tree)
library(tidyverse)
setwd("~/bugbank_data")

# input files
bb_f <- "./ukb_sgss_extract_refined.csv"
sgss_tax_tree_f <- "./sgss_pathogen_tax_tree_13032023.rds"
ukb_f <- "/well/bag/wilson/ukb/ukb41482.ukb41376.fields.RData"

# load file
bb <- read.csv(bb_f, sep = "\t")
sgss_tax_tree <- readRDS(sgss_tax_tree_f)
load(ukb_f)
ukb_eids <- bd[, "f.eid"]


# travel the tree and get eids for each node (for each speciemn type)
sgss_tax_tree$Do(function(node) {
    # filter the bb file by infection relevant to the node
    cur_bb <- bb %>% filter(ORGANISM_SPECIES_NAME %in% node$cum_origin_name)
    specimens <- unique(cur_bb$SPECIMEN_TYPE_DESC)
    # collect eids for each specimen
    specimen_eid <- map(specimens, function(specimen) {
        # filter the bb further by specimen type
        cur_bb_sub <- cur_bb %>% filter(SPECIMEN_TYPE_DESC == specimen)
        return(unique(cur_bb_sub$UKB_EID))
    })
    names(specimen_eid) <- specimens
    # add an additional specimen "all" that pools all specimen together
    specimen_eid[[length(specimen_eid) + 1]] <- unique(unlist(specimen_eid))
    names(specimen_eid)[length(specimen_eid)] <- "all"

    node$specimen_eid <- specimen_eid
})

# traverse the tree, if a node-specimen combination has
# 1. more than 110% case eids than any of its children node-specimen combination
# 2. more than 100 cases eids
# then we can include the node-specimen combination into the analysis
in_analysis <- sgss_tax_tree$Get(function(node) {
    specimens <- names(node$specimen_eid)
    specimen_in_analysis <- map_lgl(specimens, function(specimen) {
        cur_eids <- node$specimen_eid[[specimen]]
        if (length(cur_eids) < 100) {
            return(FALSE)
        }
        # for each child, check the eids of the same node-specimen combination
        # if the eids of any child is more than 95% of the current node, the current node-specimen should not be included in the analysis
        for (child in node$children) {
            # get the eid of node-specimen combination for the child
            child_eids <- child$specimen_eid[[specimen]]
            if (length(cur_eids) * 0.95 < length(child_eids)) {
                return(FALSE)
            }
        }
        return(TRUE)
    })

    # name the result in the format of texa_level:texa:specimen
    names(specimen_in_analysis) <- paste0(node$name, ":", specimens)
    return(specimen_in_analysis)
})
in_analysis <- unlist(in_analysis)
in_analysis <- names(in_analysis)[in_analysis]
