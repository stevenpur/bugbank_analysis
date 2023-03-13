# This R code is used to assign the taxonomy of the pathogen identified in bugbank and in HES.

library(tidyverse)
# library(taxize)
library(taxizedb)
setwd("~/bugbank_data/")
tax_class <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")


###############
# FOR BUGBANK #
###############

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
        next
    }
    # search NCBI database for the taxonomy with some exception handling
    taxid <- name2taxid(pathogen, out_type = "summary")$id[1]
    search_result <- classification(taxid)
    # Note: I have checked instances where there are more than one id:
    # for all of these cases, the first id is the correct taxon

    # error handling
    if (is.na(search_result)) {
        stop(paste0("can't find taxonomy for pathogen ", pathogen))
    } else if (length(search_result) != 1) {
        stop(paste0("there is a problem finding taxonomy for pathogen ", pathogen))
    }

    # compiling results
    search_result[[1]] <- search_result[[1]][which(search_result[[1]]$rank %in% tax_class), ]
    rownames(search_result[[1]]) <- search_result[[1]]$rank
    search_results[[pathogen]] <- c(search_result[[1]][tax_class, 1], origin_pathogen)
}
pathogen_taxonomy <- data.frame(do.call(rbind, search_results), stringsAsFactors = F)
colnames(pathogen_taxonomy) <- c(tax_class, "origin_name")
# write result
write.table(pathogen_taxonomy, "./bb_pathogen_taxonomy_13032023.tsv", sep = "\t", quote = F, col.names = T, row.names = F)

###########
# FOR HES #
###########

# input file
pathogen_icd10_desc_file <- "./icd10_pathogen_description_13032023.tsv"
# load input file
pathogen_icd10_desc <- read.csv(pathogen_icd10_desc_file, sep = "\t")
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
pathogen_taxonomy <- data.frame(do.call(rbind, search_results), stringsAsFactors = F)
colnames(pathogen_taxonomy) <- c(tax_class, "origin_name")
# write result
write.table(pathogen_taxonomy, "./hes_pathogen_taxonomy_13032023.tsv", sep = "\t", quote = F, col.names = T, row.names = F)