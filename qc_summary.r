# save working image
# save.image(file = "summary_table_image.RData")
# load working image
 load("summary_table_image.RData")

library(tidyverse)
bugbank_data_dir = "/well/bag/clme1992/bugbank_data/"
ukb_data_dir = "/well/bag/wilson/ukb/"
gwas_wrkdir = "/well/bag/clme1992/saige_pipe_test"
outdir = "/well/bag/bugbank_data/"


# At each step of QC, track the change in the infection cases

resolve.symlink <- function(x) {
    y <- Sys.readlink(x)
    if (any(is.na(y))) {
        stop("Could not resolve symlink ", x)
    }
    ifelse(y == "", x, file.path(dirname(x), y))
}

# configuration file
source("~/.saige_pipe.config")

lg <- list()
lg$bd_RDdata_file <- paste0(config$ukb.derived.dir, "/ukb41482.ukb41376.fields.RData")
lg$bd_not_lost2followup_file <- paste0(config$ukb.derived.dir, "/ukb41482.English-not-lost-to-followup-8-April-2020.txt")
lg$bed_sample_qc_file <- paste0(config$ukbdir, "/v2/qc/ukb_sqc_v2.txt")
lg$withdrawn_eid_file <- paste0(config$bbdatadir, "/w53100_2023-04-25.csv")
# Pre-computed eids for the bed-format genotypes
lg$bed_eid_file <- paste0(config$ukb.derived.dir, "/analysis.bed.eids.txt")
# Individuals with first degree relatives
lg$remrels_file <- paste0(config$ukb.derived.dir, "/ukb41482.English-remove-first-degree-relatives.eids.txt")
# panUKB ancestral files
lg$pan_ukb_file <- paste0(config$panukb.dir, "/Files for retman/all_pops_non_eur_pruned_within_pop_pc_covs.tsv")
lg$pan_ukb_bridge_file <- paste0(config$panukb.dir, "/ukb53100bridge31063.txt")

### load input files ###
system.time(load(lg$bd_RDdata_file))
all_eids <- bd[, "f.eid"]
bd_not_lost2followup <- scan(lg$bd_not_lost2followup_file, what = "logical") == "TRUE"
withdrawn_eid <- scan(lg$withdrawn_eid_file)
# Sample QC
bed_sample_qc <- read.csv(lg$bed_sample_qc_file, sep = " ")
# The corresponding eids
bed_eid <- scan(lg$bed_eid_file)
# Convert to bd_eid order
sample_qc <- bed_sample_qc[match(all_eids, bed_eid), ]
# Close (first degree) relatives
remrels <- scan(lg$remrels_file)
# load and match the panukb data
panukb <- read.csv(lg$pan_ukb_file, sep = "\t")[, c("s", "pop")]
panukb_bridge <- read.csv(lg$pan_ukb_bridge_file, sep = " ", header = F)
bridge_matched <- panukb_bridge[match(panukb$s, panukb_bridge[, 2]), ]
panukb$eid <- bridge_matched[, 1]
panukb_matched <- panukb[match(all_eids, panukb$eid), ]

### assign assessment centre data ###
f_assesscentre <- "f.54.0.0"
assess_centre_England <- c(
    11012, 11021, 11011, 11008, 11024, 11020, 11018, 11010, 11016, 
    11001, 11017, 11009, 11013, 11002, 11007, 11014, 10003, 11006, 
    11025, 11026, 11027, 11028 
)

# load SGSS
sgss_file = paste0(bugbank_data_dir, "ukb_sgss_extract_refined.csv")
sgss = read.csv(sgss_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# load HES
start_time = Sys.time()
hes_file = paste0(ukb_data_dir, "hes/hesin_diag.latest.txt.gz")
hes = read.table(hes_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
end_time = Sys.time()
print(paste0("Time taken to load hes: ", end_time - start_time))

# get hes infect related records
# load icd10 data 
icd10_desc_file <- paste0(bugbank_data_dir, "pathogen_to_unique_icd10.tsv")
icd10_desc <- read.table(icd10_desc_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# subset hes to only infection related
infect_icd10_codes <- unique(unlist(strsplit(icd10_desc$icd10, split = ",")))
hes_infect <- hes[hes$diag_icd10 %in% infect_icd10_codes, ]


filter_lst <- list()
filter <- rep(T, nrow(bd))
filter_lst[[1]] <- filter
names(filter_lst)[1] <- "all"

filter <- filter & (bd[, f_assesscentre] %in% assess_centre_England)
filter_lst[[2]] <- filter
names(filter_lst)[2] <- "in assessment centre England"

filter <- filter & bd_not_lost2followup
filter <- filter & !(all_eids %in% withdrawn_eid)
filter_lst[[3]] <- filter
names(filter_lst)[3] <- "not lost to followup"

filter <- filter & sample_qc$putative.sex.chromosome.aneuploidy == 0
filter_lst[[4]] <- filter
names(filter_lst)[4] <- "no aneuploidy in sex chromosome"

filter <- filter & sample_qc$Submitted.Gender == sample_qc$Inferred.Gender
filter_lst[[5]] <- filter
names(filter_lst)[5] <- "reported sex matches genetic sex"

filter <- filter & sample_qc$het.missing.outliers == 0
filter_lst[[6]] <- filter
names(filter_lst)[6] <- "no het missing outliers"

filter <- filter & sample_qc$excluded.from.kinship.inference == 0
filter_lst[[7]] <- filter
names(filter_lst)[7] <- "not excluded from kinship inference"

filter <- filter & sample_qc$excess.relatives == 0 & is.na(match(all_eids, remrels))
filter_lst[[8]] <- filter
names(filter_lst)[8] <- "no excess nor 1st degree relatives"

filter <- filter & sample_qc$in.Phasing.Input.chr1_22 == 1 & 
    sample_qc$in.Phasing.Input.chrX == 1 &
    sample_qc$in.Phasing.Input.chrXY == 1
filter_lst[[9]] <- filter
names(filter_lst)[9] <- "in phasing input"

#filter <- filter & all_eids %in% panukb_matched$eid[which(panukb_matched$pop == "EUR")]
#filter_lst[[9]] <- filter
#names(filter_lst)[9] <- "in panukb EUR"



sgss_filter <- sgss
hes_filter <- hes_infect
ukb_filter <- list()
for (i in 1:length(filter_lst)) {
    filter <- filter_lst[[i]]
    eids_filtered <- unique(all_eids[filter])
    eids_filtered <- eids_filtered[!is.na(eids_filtered)]
    sgss_filter <- sgss_filter[sgss_filter$UKB_EID %in% eids_filtered, ]
    hes_filter <- hes_filter[hes_filter$eid %in% eids_filtered, ]
    sgss_cnt <- length(unique(sgss_filter$UKB_EID))
    hes_cnt <- length(unique(hes_filter$eid))
    ukb_filter[[i]] <- c(length(eids_filtered), sgss_cnt, hes_cnt)
}
ukb_filter <- do.call(rbind, ukb_filter)
rownames(ukb_filter) <- names(filter_lst)
colnames(ukb_filter) <- c("record_cnt", "sgss_cnt", "hes_cnt")
ukb_filter

# save the table
ukb_outfile <- paste0(outdir, "ukb_filter_summary.tsv")
write.table(ukb_filter, file = "ukb_filter_summary.tsv", sep = "\t", quote = FALSE)

# for SGSS
sgss_tax_file <- paste0(bugbank_data_dir, "bb_pathogen_taxonomy_13032023.tsv")
sgss_tax <- read.table(sgss_tax_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# ------------------------------------------------------------
# create a dictionary from origin_name to most specific taxa
# ------------------------------------------------------------

# sanity check: no NA nor duplicated origin_name
sum(is.na(sgss_tax$origin_name))
sum(duplicated(sgss_tax$origin_name))
origin_name_to_tax <- list()

# for each origin_name, find the most specific tax

# get the tax levels
tax_levs <- colnames(sgss_tax)[-ncol(sgss_tax)]
# remove the last column (which is the origin_name) for later use
sgss_tax_only <- sgss_tax[, -ncol(sgss_tax)]

for (i in 1:nrow(sgss_tax)) {
    origin_name <- sgss_tax$origin_name[i]
    # remove any NA and the last column (which is the origin_name)
    tax <- sgss_tax_only[i, ]
    tax_lev <- tax_levs[!is.na(tax)]
    tax <- tax[!is.na(tax)]
    specific_tax <- tail(tax, 1)
    specific_tax_lev <- tail(tax_lev, 1)
    origin_name_to_tax[[origin_name]] <- c(specific_tax, specific_tax_lev)
}



# filter1: starting from all
filter_lst <- list()
filter <- rep(T, nrow(sgss))

filter_lst[[1]] <- filter
names(filter_lst)[1] <- "all"

# filter2: only include individuals that passed the QC
filter_lst[[2]] <- filter & sgss$UKB_EID %in% eids_filtered
names(filter_lst)[2] <- "individual that passed QC"

# filter3: only include infection that have a pathogen label in species level
sgss$tax_lev <- map_chr(sgss$ORGANISM_SPECIES_NAME, function(origin_name) {
    tax_lev <- origin_name_to_tax[[origin_name]][2]
})
filter_lst[[3]] <- filter_lst[[2]] & sgss$tax_lev == "species"
names(filter_lst)[3] <- "species level pathogen"

# filter4: only include pathogen that have infection cases > 100
sgss$tax <- map_chr(sgss$ORGANISM_SPECIES_NAME, function(origin_name) {
    species_name <- origin_name_to_tax[[origin_name]][1]
})
species_name_qc <- unique(sgss$tax[filter_lst[[3]]])
# get the frequency of each species by the number of unique individuals
species_name_freq <- map_int(species_name_qc, function(species_name) {
    cur_filter <- (sgss$tax == species_name) & filter_lst[[3]]
    cur_eids <- sgss$UKB_EID[cur_filter]
    return(length(unique(cur_eids)))
})
names(species_name_freq) <- species_name_qc
species_name_gwas <- names(species_name_freq)[which(species_name_freq > 100)]
filter_lst[[4]] <- filter_lst[[3]] & (sgss$tax %in% species_name_gwas)
names(filter_lst)[4] <- "pathogen with > 100 cases"

sgss_filter_table <- list()
for(i in 1:length(filter_lst)) {
    filter <- filter_lst[[i]]
    record_cnt <- sum(filter, na.rm = TRUE)
    ind_cnt <- length(unique(sgss$UKB_EID[filter]))
    tax_cnt <- length(unique(sgss$tax[filter]))
    sgss_filter_table[[i]] <- c(record_cnt, ind_cnt, tax_cnt)
}
sgss_filter_table <- do.call(rbind, sgss_filter_table)
rownames(sgss_filter_table) <- names(filter_lst)
colnames(sgss_filter_table) <- c("record_cnt", "ind_cnt", "tax_cnt")
sgss_filter_table

# save the table
sgss_outfile <- paste0(outdir, "sgss_filter_summary.tsv")
write.table(sgss_filter_table, file = "sgss_filter_summary.tsv", sep = "\t", quote = FALSE)

# sanity check 1
# see if the resulting number of tax is the same as the number of GWAS
gwas_wrkdir <- "/well/bag/clme1992/saige_pipe_test"
gwas_files <- list.files(gwas_wrkdir, pattern = "summary.05062023_sgss_species.*all.txt.gz", full.names = TRUE)
gwas_files_species <- gsub(".*summary.05062023_sgss_species.", "", gwas_files)
gwas_files_species <- gsub(".species.all.txt.gz", "", gwas_files_species)
gwas_files_species <- gsub("_", " ", gwas_files_species)
print(paste0("number of pathogen in the actual GWAS: ", length(gwas_files_species)))
print(paste0("number of pathogen that are qc here: ", length(species_name_gwas)))
print(setdiff(species_name_gwas, gwas_files_species))

# sanity check 2
# see if number of infection is the same for E.coli from the phenotype file
test_species <- names(species_name_freq)[which(species_name_freq > 100)]
for(cur_species in test_species){
    # get the eids from gwas phenotype file
    cur_species_fname <- gsub(" ", "_", cur_species)
    phe_file <- paste0(gwas_wrkdir, "/ukb41482.bd.gwas-pheno.05062023_sgss_species.", cur_species_fname, ".species.all.txt.gz")
    phe <- read.table(phe_file, sep = "\t", header = T, stringsAsFactors = F)
    phe_eids <- phe$IID[which(phe$pheno == 1)]

    # get the eids from sgss filtering
    cur_filter <- (sgss$tax == cur_species) & filter_lst[[3]]
    sgss_eids <- sgss$UKB_EID[cur_filter]

    # compare the two
    print(cur_species)
    print(length(setdiff(phe_eids, sgss_eids)))
}

#####################
# the main problem is that I forgot to filter by n>100 when running the GWAS, so the QC results here would be a subset of the actual GWAS. I should used the subset GWAS in the final analysis
#####################

# For HES

# load hes_diag
start_time = Sys.time()
hes_file = paste0(ukb_data_dir, "hes/hesin_diag.latest.txt.gz")
hes = read.table(hes_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
end_time = Sys.time()
print(paste0("Time taken to load hes: ", end_time - start_time))

# load icd10 data 
icd10_desc_file <- paste0(bugbank_data_dir, "pathogen_to_unique_icd10.tsv")
icd10_desc <- read.table(icd10_desc_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# subset hes to only infection related
infect_icd10_codes <- unique(unlist(strsplit(icd10_desc$icd10, split = ",")))
hes_infect <- hes[hes$diag_icd10 %in% infect_icd10_codes, ]

# create a dictionary of ICD-10 to pathogen mapping
icd_to_pathogen = list()
for (i in 1:nrow(icd10_desc)) {
  cur_icd10s = unlist(strsplit(icd10_desc$icd10[i], ","))
  for (icd10 in cur_icd10s) {
    # check if the icd10 is already in the dictionary
    # if so, report an error
    if (icd10 %in% names(icd_to_pathogen)) {
      print(paste0("Error: ", icd10, " already in the dictionary"))
    }
    icd_to_pathogen[[icd10]] = c(icd10_desc$org_name[i], icd10_desc$tax_lev[i])
  }
}

# filter1: starting from all
filter_lst <- list()
filter <- rep(T, nrow(hes_infect))
filter_lst[[1]] <- filter
names(filter_lst)[1] <- "all"

# filter2: only include individuals that passed the QC
filter_lst[[2]] <- filter & hes_infect$eid %in% eids_filtered
names(filter_lst)[2] <- "individual that passed QC"

#filter3: only include infection that have a pathogen label in species level
hes_infect$tax_lev <- map_lgl(hes_infect$diag_icd10, function(icd10) {
  tax_lev <- icd_to_pathogen[[icd10]][2]
  return(tax_lev == "species")
})
filter_lst[[3]] <- filter_lst[[2]] & hes_infect$tax_lev
names(filter_lst)[3] <- "species level pathogen"

# filter4: only include pathogen that have infection cases > 100
hes_infect$tax <- map_chr(hes_infect$diag_icd10, function(icd10) {
    species_name <- icd_to_pathogen[[icd10]][1]
    return(species_name)  
})
species_name_qc <- unique(hes_infect$tax[filter_lst[[3]]])
# get the frequency of each species by the number of unique individuals
species_name_freq <- map_int(species_name_qc, function(species_name) {
    cur_filter <- (hes_infect$tax == species_name) & filter_lst[[3]]
    cur_eids <- hes_infect$eid[cur_filter]
    return(length(unique(cur_eids)))
})
names(species_name_freq) <- species_name_qc
species_name_gwas <- names(species_name_freq)[which(species_name_freq > 100)]
filter_lst[[4]] <- filter_lst[[3]] & (hes_infect$tax %in% species_name_gwas)
names(filter_lst)[4] <- "pathogen with > 100 cases"

# gather the filter into a table
hes_filter_table <- list()
for(i in 1:length(filter_lst)) {
    filter <- filter_lst[[i]]
    record_cnt <- sum(filter, na.rm = TRUE)
    ind_cnt <- length(unique(hes_infect$eid[filter]))
    tax_cnt <- length(unique(hes_infect$tax[filter]))
    hes_filter_table[[i]] <- c(record_cnt, ind_cnt, tax_cnt)
}
hes_filter_table <- do.call(rbind, hes_filter_table)
rownames(hes_filter_table) <- names(filter_lst)
colnames(hes_filter_table) <- c("record_cnt", "ind_cnt", "tax_cnt")
hes_filter_table

# sanity check, see if the resulting number of tax is the same as the number of GWAS
gwas_wrkdir <- "/well/bag/clme1992/saige_pipe_test"
gwas_files <- list.files(gwas_wrkdir, pattern = "summary.05062023_hes_species.*.txt.gz", full.names = TRUE)
gwas_files_species <- gsub(".*summary.05062023_hes_species.regenie.", "", gwas_files)
gwas_files_species <- gsub(".species.txt.gz", "", gwas_files_species)
gwas_files_species <- gsub("_", " ", gwas_files_species)
print(length(gwas_files_species))
print(setdiff(gwas_files_species, species_name_gwas))
print(setdiff(species_name_gwas, gwas_files_species))



