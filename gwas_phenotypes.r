library(tidyverse)

# configuration file
source("~/.saige_pipe.config")

# command-line arguments
lg <- list()
help <- paste(
    "Usage: Rscript hgi-phenotypes-bugbank.R stem pathogen_file source",
    sep = "\n"
)

# Argument
lg$args <- commandArgs(trailingOnly = TRUE)
print("Arguments:")
print(lg$args)
if (!(length(lg$args) %in% c(2, 3))) {
    cat(help, sep = "\n")
    stop("\nIncorrect usage\n")
}

lg$stem <- as.character(lg$args[1])
if (is.na(lg$stem) | lg$stem == "") stop("Nonsensical stem")
lg$pathogen_file <- as.character(lg$args[2])
if (is.na(lg$pathogen_file) | lg$pathogen_file == "") stop("Nonsensical pathogen file")
lg$source <- as.character(lg$args[3])
if (is.na(lg$source) | lg$source == "") stop("Please specify data source (hes/sgss)")


### Input and output files ###

# input files

lg$bugbank_file <- paste0(config$bbdatadir, "/ukb_sgss_extract_refined.csv")
if (lg$source == "sgss") {
    lg$pathogen_taxonomy_file <- paste0(config$bbdatadir, "/bb_pathogen_taxonomy_13032023.tsv")
} else if (lg$source == "hes") {
    lg$pathogen_taxonomy_file <- paste0(config$bbdatadir, "/hes_pathogen_taxonomy_fromRaw_13032023.tsv")
}


lg$hesin_diag_file <- paste0(config$ukb.derived.dir, "/hes/hesin_diag.latest.txt.gz")
lg$pathogen_icd10_file <- paste0(config$bbdatadir, "/pathogen_to_unique_icd10.tsv")
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

# output files
lg$log_rds_outfile <- paste0(config$wrkdir, "/log.ukb41482.bd.gwasdata.", lg$stem, ".rds")
lg$data_outfile <- paste0(config$wrkdir, "/ukb41482.bd.gwasdata.", lg$stem, ".txt")

tryCatch(
    {
        ### load input files ###
        pathogen_tb <- read.csv(lg$pathogen_file, sep = "\t")
        pathogen_taxonomy <- read.csv(lg$pathogen_taxonomy_file, sep = "\t")
        pathogen_icd10 <- read.csv(lg$pathogen_icd10_file, sep = "\t")
        bugbank_data <- read.csv(lg$bugbank_file, sep = "\t")
        hes_diag <- read.csv(lg$hesin_diag_file, sep = "\t")
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


        ### get phenotype ###

        # find cases eid for each pathogen-specimen combination
        if (lg$source == "sgss") {
            cases_eids <- map(1:nrow(pathogen_tb), function(i) {
                tax_lev <- pathogen_tb$tax_lev[i]
                name <- pathogen_tb$name[i]
                specimen <- pathogen_tb$specimen[i]
                pathogens <- pathogen_taxonomy$origin_name[which(pathogen_taxonomy[[tax_lev]] == name)]
                if (tolower(specimen) == "all") {
                    sub_bb_data <- bugbank_data %>% filter(ORGANISM_SPECIES_NAME %in% pathogens)
                    return(unique(sub_bb_data$UKB_EID))
                } else {
                    sub_bb_data <- bugbank_data %>% filter(
                        ORGANISM_SPECIES_NAME %in% pathogens &
                            SPECIMEN_GROUP_DESC == specimen
                    )
                    return(unique(sub_bb_data$UKB_EID))
                }
            })
            names(cases_eids) <- gsub(" ", "_", paste0(pathogen_tb$name, ".", pathogen_tb$tax_lev, ".", pathogen_tb$specimen))
        } else if (lg$source == "hes") {
            cases_eids <- map(1:nrow(pathogen_tb), function(i) {
                tax_lev <- pathogen_tb$tax_lev[i]
                name <- pathogen_tb$name[i]
                # to get the icd10s, we need the taxonomy level name instead of origin name
                icd10s <- pathogen_icd10$icd10[which(pathogen_icd10$org_name %in% name)]
                icd10s <- unique(unlist(strsplit(icd10s, ",")))
                result_eids <- unique(hes_diag$eid[which(hes_diag$diag_icd10 %in% icd10s)])
            })
            names(cases_eids) <- gsub(" ", "_", paste0(pathogen_tb$name, ".", pathogen_tb$tax_lev))
        }

        ### find the confounding individuals (anyone who's ever been infected) ###

        # for hes
        infect_icd10 <- unique(unlist(strsplit(pathogen_icd10$icd10, ",")))
        hes_diag_infect <- hes_diag %>% filter(diag_icd10 %in% infect_icd10)
        hes_infect_eid <- unique(hes_diag_infect$eid)
        # for sgss
        sgss_infect_eid <- unique(bugbank_data$UKB_EID)
        # all together
        infect_eid <- union(hes_infect_eid, sgss_infect_eid)



        ### assign assessment centre data ###
        f.assesscentre <- "f.54.0.0"
        assess_centre_England <- c(
            11012, # 	Barts
            11021, # 	Birmingham
            11011, # 	Bristol
            11008, # 	Bury
            # 11003	Cardiff
            11024, # 	Cheadle (revisit)
            11020, # 	Croydon
            # 11005	Edinburgh
            # 11004	Glasgow
            11018, # 	Hounslow
            11010, # 	Leeds
            11016, # 	Liverpool
            11001, # 	Manchester
            11017, # 	Middlesborough
            11009, # 	Newcastle
            11013, # 	Nottingham
            11002, # 	Oxford
            11007, # 	Reading
            11014, # 	Sheffield
            10003, # 	Stockport (pilot)
            11006, # 	Stoke
            # 11022	Swansea
            # 11023	Wrexham
            11025, # 	Cheadle (imaging)
            11026, # 	Reading (imaging)
            11027, # 	Newcastle (imaging)
            11028 # 	Bristol (imaging)
        )


        ### individual filtering ###
        # filter out samples not suitable to include in analysis
        f_assesscentre <- "f.54.0.0"
        filter <- bd[, f_assesscentre] %in% assess_centre_England &
            !(all_eids %in% withdrawn_eid) &
            bd_not_lost2followup &
            sample_qc$het.missing.outliers == 0 &
            sample_qc$putative.sex.chromosome.aneuploidy == 0 &
            sample_qc$Submitted.Gender == sample_qc$Inferred.Gender &
            sample_qc$excluded.from.kinship.inference == 0 &
            sample_qc$excess.relatives == 0 &
            sample_qc$in.Phasing.Input.chr1_22 == 1 &
            sample_qc$in.Phasing.Input.chrX == 1 &
            sample_qc$in.Phasing.Input.chrXY == 1 &
            is.na(match(all_eids, remrels))
        filter[is.na(filter)] <- FALSE


        ### assign phenotype ###
        pheno <- matrix(0, nrow = length(all_eids), ncol = nrow(pathogen_tb))
        colnames(pheno) <- names(cases_eids)
        for (pheno_name in colnames(pheno)) {
            # remove individuals with infection history from controls
            pheno[which(all_eids %in% infect_eid), pheno_name] <- NA
            # assign cases
            pheno[which(all_eids %in% cases_eids[[pheno_name]]), pheno_name] <- 1
        }
        pheno[which(!filter), ] <- NA

        # handle the weird character that appear in specimen name
        colnames(pheno) <- gsub("\\(|\\)", "", gsub("/|-", "_", colnames(pheno)))

        # construct pheno file
        data <- data.frame(eid = all_eids, pheno)
        data$sex <- as.character(bd[, f.sex])
        data$age <- (2021 - bd[, f.birthyear])
        data$ancestry <- c("oth", "eur")[1 + sample_qc$in.white.British.ancestry]
        data$ancestry_panukb <- panukb_matched$pop

        # add other info to phenotype table
        print("NA counts per column in data:")
        print((lg$na.counts.columns.data <- colSums(is.na(data))))
        # Count number of cases
        print("number of cases:")
        print((lg$cases.counts <- colSums(data[, colnames(pheno)], na.rm = T)))
        # Count number of cases
        print("number of cases:")
        print((lg$cases.counts <- colSums(data[, colnames(pheno)], na.rm = T)))

        write.table(data, file = lg$data_outfile, row.names = F, quote = F, sep = "\t")
    },
    finally = {
        # On error or clean exit
        saveRDS(lg, file = lg$log_rds_outfile)
    }
)
