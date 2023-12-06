# This R code is used to assign the taxonomy of the pathogen identified in SGSS (bugbank) and in HES.

library(tidyverse)
library(taxizedb)
library(data.tree)
setwd("~/bugbank_data/")
tax_class <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")

# define function
tax2tree <- function(tax_tb) {
    tax_path <- apply(tax_tb, 1, function(x) {
        row <- na.exclude(unlist(x))
        paste0("life/", paste(c(paste0(names(row), ":", row)), collapse = "/"))
    })
    tree_data <- tax_tb
    tree_data$pathString <- tax_path
    tax_tree <- as.Node(tree_data)
    return(tax_tree)
}


#############################
# Bugbank --Assign taxonomy #
#############################

# input file
species_tb_file <- "result/bugbank_species_cnt_16112021.tsv"
species_altname_file <- "./bb_pathogen_renamed.tsv"

# load input file
species_tb <- read.table(species_tb_file, sep = "\t", header = T, stringsAsFactors = F)
species_altname <- read.table(species_altname_file, sep = "\t", stringsAsFactors = F, header = F)
colnames(species_altname) <- c("origin_name", "alt_name")

# clean species name column
origin_pathogens <- species_tb$species_name
pathogens <- gsub("GROUP.*", "", origin_pathogens)
pathogens <- gsub(" \\(.*\\)", "", pathogens)
pathogens <- gsub(" OTHER.*", "", pathogens)
pathogens <- gsub("^ | $", "", pathogens)
pathogens <- gsub(" SP$| SUBSP.*", "", pathogens)
pathogens <- gsub(" COMPLEX$", "", pathogens)
pathogens <- map_chr(pathogens, function(pathogen) {
    alt_row <- which(species_altname$origin_name == pathogen)
    if (length(alt_row) != 0) {
        return(species_altname$alt_name[alt_row])
    }
    return(pathogen)
})

# search for taxonomy #
search_results <- list()
for (i in seq_len(length(pathogens))) {
    origin_pathogen <- origin_pathogens[i]
    pathogen <- pathogens[i]
    message(paste0(i, "/", length(pathogens)))
    if (is.na(pathogen)) {
        # in this instance the pathogen cannot be identified
        # the record should be "unknown", which still counts as an infection but with no organim assign
        search_result <- list(data.frame(name = NA, rank = tax_class))
        rownames(search_result[[1]]) <- tax_class
    } else {
        # search NCBI database for the taxonomy with some exception handling
        taxid <- name2taxid(pathogen, out_type = "summary")$id[1]
        search_result <- classification(taxid)
        # Note: I have checked instances where there are more than one id:
        # for all of these cases, the first id is the correct taxon
    }

    # error handling
    if (is.na(search_result)) {
        stop(paste0("can't find taxonomy for pathogen ", pathogen))
    } else if (length(search_result) != 1) {
        stop(paste0("there is a problem finding taxonomy for pathogen ", origin_pathogen))
    }

    # compiling results
    search_result[[1]] <- search_result[[1]][which(search_result[[1]]$rank %in% tax_class), ]
    rownames(search_result[[1]]) <- search_result[[1]]$rank
    search_results[[origin_pathogen]] <- c(search_result[[1]][tax_class, 1], origin_pathogen)
}

# compile all taxonomy results into a table
pathogen_taxonomy <- data.frame(do.call(rbind, search_results), stringsAsFactors = F)
colnames(pathogen_taxonomy) <- c(tax_class, "origin_name")

write.table(pathogen_taxonomy, "./bb_pathogen_taxonomy_13032023.tsv", sep = "\t", quote = F, col.names = T, row.names = F)


###################################
# Bugbank -- build a species tree #
###################################

# turn the taxonomy table into a tree (origin_name should not be included as a node in the species tree)
tree_data <- pathogen_taxonomy
tree_data$pathString <- map_chr(1:nrow(tree_data), function(i) {
    row <- tree_data[i, tax_class]
    row <- na.exclude(unlist(row))
    if (length(row) == 0) {
        # the organism is unknown (all NA in the taxonomy)
        return("life:all_org")
    }
    paste0("life:all_org/", paste(c(paste0(names(row), ":", row)), collapse = "/"))
})
tax_tree <- as.Node(tree_data)

# 1. add origin_names information to corresponding nodes
# 2. add cumulative origin_names (including origin_names from all descendants)
add_origin_name_info <- function(node) {
    cum_origin_name <- NULL
    if (node$pathString %in% tree_data$pathString) {
        # this node has at least one origin_name, add this info to the node
        origin_names <- tree_data$origin_name[which(tree_data$pathString == node$pathString)]
        node$origin_name <- origin_names
        node$cum_origin_name <- origin_names
    }
    if (!node$isLeaf) {
        # if internal node, than cum_origin_names is origin_names + children_origin_names
        for (child in node$children) {
            add_origin_name_info(child)
            node$cum_origin_name <- c(node$cum_origin_name, child$cum_origin_name)
        }
    }
}

add_origin_name_info(tax_tree)

# write result
saveRDS(tax_tree, "./sgss_pathogen_tax_tree_13032023.rds")





##########################
# HES -- assign taxonomy #
##########################

# input file
pathogen_icd10_desc_file <- "./icd10_pathogen_description_13032023.tsv"
# load input file
pathogen_icd10_desc <- read.csv(pathogen_icd10_desc_file, sep = "\t", stringsAsFactors = F)
# clean pathogens names
pathogens <- unique(pathogen_icd10_desc$org_name)

# start taxonomy search
search_results <- list()
for (i in 1:length(pathogens)) {
    pathogen <- pathogens[i]
    message(paste0(i, "/", length(pathogens)))
    if (is.na(pathogen)) {
        next
    }

    # search NCBI database for the taxonomy with some exceptions handling
    if (pathogen == "Spirillum minus") {
        # Spirillum minus (nor its synonyms) can't be found in NCBI
        # But it is appear repeatedly in literature, so I create the species
        taxid <- name2taxid("Spirillum")
        search_result <- classification(taxid)
        search_result[[1]] <- rbind(search_result[[1]], c("Spirillum minus", "species", "NA"))
    } else if (pathogen == "Treponema carateum") {
        # Treponema carateum (nor its synonyms) can't be found in NCBI
        # But it is appear repeatedly in literature, so I create the species
        taxid <- name2taxid("Spirillum")
        taxid <- name2taxid("Treponema")
        search_result <- classification(taxid)
        search_result[[1]] <- rbind(search_result[[1]], c("Treponema carateum", "species", "NA"))
    } else if (pathogen == "Prions") {
        # This is not a organism, we would not include this in our analysis
        next
    } else {
        taxid <- name2taxid(pathogen, out_type = "summary")$id[1]
        search_result <- classification(taxid)
        # Note: I have checked instances where there are more than one id:
        # for all of these cases, the first id is the correct taxon
    }

    # error handling
    if (is.na(search_result)) {
        stop(paste0("can't find taxonomy for pathogen ", pathogen))
    } else if (length(search_result) != 1) {
        stop(paste0("there is a problem finding taxonomy for pathogen ", pathogen))
    }

    # compiling results
    search_result[[1]] <- search_result[[1]][which(search_result[[1]]$rank %in% tax_class), ]
    rownames(search_result[[1]]) <- search_result[[1]]$rank
    search_results[[pathogen]] <- c(search_result[[1]][tax_class, 1], pathogen)
}

# compile all pathogen taxonomy into a table
pathogen_taxonomy <- data.frame(do.call(rbind, search_results), stringsAsFactors = F)
colnames(pathogen_taxonomy) <- c(tax_class, "origin_name")

# more exception handling
pathogen_taxonomy$species[which(pathogen_taxonomy$species == "unidentified")] <- NA

# write result to file
write.table(pathogen_taxonomy, "./hes_pathogen_taxonomy_fromRaw_13032023.tsv", sep = "\t", quote = F, col.names = T, row.names = F)

################################
# HES -- build a taxonomy tree #
################################

# turn the taxonomy table into a tree (origin_name should not be included as a node in the species tree)
tree_data <- pathogen_taxonomy
tree_data$pathString <- map_chr(1:nrow(tree_data), function(i) {
    row <- tree_data[i, tax_class]
    row <- na.exclude(unlist(row))
    if (length(row) == 0) {
        # the organism is unknown (taxonomy assignment is all NA an all levels)
        return("life:all_org")
    }
    paste0("life:all_org/", paste(c(paste0(names(row), ":", row)), collapse = "/"))
})
tax_tree <- as.Node(tree_data)

# write result
saveRDS(tax_tree, "./hes_pathogen_tax_tree_fromRaw_13032023.rds")


#################################
# HES -- species tree refinment #
#################################

# build a dictionary to map pathgoen to icd10
path_to_icd10 <- list()
for (i in 1:nrow(pathogen_icd10_desc)) {
    pathogen <- pathogen_icd10_desc$org_name[i]
    icd10s <- strsplit(pathogen_icd10_desc$UKB_code[i], ",")[[1]]
    path_to_icd10[[pathogen]] <- unique(c(path_to_icd10[[pathogen]], icd10s))
}

# we can travel and modify data.tree as in pointers in C
modify_tree <- function(node) {
    # add the following info to the tax tree:
    # 1. origin_name to the corresponding nodes
    # 2. cum_icd10: ICD10 codes that is the accumulation of the current node and all its descedents
    # 3. unique_icd10: ICD10 codes that is unique to the node (excluding ICD10 already present in its descendents)

    # add origin name
    if (node$pathString %in% tree_data$pathString) {
        origin_names <- tree_data$origin_name[which(tree_data$pathString == node$pathString)]
        node$origin_name <- origin_names
    }

    # add cumulative and unique ICD10 info
    node$cum_icd10 <- vector()
    node$uni_icd10 <- vector()
    if ("origin_name" %in% names(node)) {
        # if this node is mapped to any origin_name:
        # (a) add the icd10 of those origin_names into cum_icd10
        # (b) uni_icd10 should be the icd10 of those origin_names minus all children icd10
        for (origin_name in node$origin_name) {
            node$cum_icd10 <- c(node$cum_icd10, path_to_icd10[[origin_name]])
            node$uni_icd10 <- c(node$uni_icd10, path_to_icd10[[origin_name]])
        }
    }

    if (node$isLeaf) {
        # if it is a leaf node, add icd10 codes to both cumulative and unique icd10
        # exception handling
        if (!"origin_name" %in% names(node)) {
            stop(paste0("there is a leaf node that does not have at least one corresponding origin_name: ", node$pathString))
        }
    } else {
        # not a leaf node, add all icd10 from children to cum_icd10
        children_icd10 <- vector()
        for (child in node$children) {
            children_icd10 <- c(children_icd10, modify_tree(child))
        }
        node$cum_icd10 <- unique(c(node$cum_icd10, children_icd10))

        # set unique icd10
        node$uni_icd10 <- setdiff(node$uni_icd10, children_icd10)
    }
    return(node$cum_icd10)
}

modify_tree(tax_tree)

# write result
saveRDS(tax_tree, "./hes_pathogen_tax_tree_refined_13032023.rds")
