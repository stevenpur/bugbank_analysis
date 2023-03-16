# this code is to add a column to the bugbank raw file so that the EID can be mapped to the EID of UKB_53100

bb_file <- "./ukb_sgss_extract_20211115.csv"
bb_bridge_file <- "./sgss_bridge.csv"

bb <- read.csv(bb_file, stringsAsFactors = F, header = F)
bb_bridge <- read.csv(bb_bridge_file, header = T)
colnames(bb) <- c(
    "SPECIMEN_NUMBER",
    "UKB_EID",
    "SPECIMEN_DATE",
    "LAB_REPORT_DATE",
    "REPORTING_LAB_NAME",
    "LAB_GEOG_NAME_CURRENT",
    "LOCAL_AUTHORITY_NAME",
    "SPECIMEN_GROUP_DESC",
    "SPECIMEN_TYPE_DESC",
    "ORGANISM_CATEGORY_DESC",
    "ORGANISM_GENUS_NAME",
    "ORGANISM_SPECIES_NAME",
    "ORGANISM_SUBSPECIES_NAME"
)

map <- match(bb$UKB_EID, bb_bridge$id_sgss)
bb$UKB_EID <- bb_bridge$eid_53100[map]

write.table(bb, "ukb_sgss_extract_refined.csv", sep = "\t", col.names = T, row.names = F, quote = F)
