# perform LD clumping to determine individual GWAS hits
import os
import sys
import pandas as pd
import logging

# check input arguments
if len(sys.argv) != 2:
    logging.error("Usage: python get_independent_hit.py stem")
    sys.exit(0)
stem = sys.argv[1]

# set parameters
sig_snp_dir = os.path.expanduser("~/bugbank_data/postgwas_regenie/sig_snp/")
ld_dir = os.path.expanduser("~/bugbank_data/postgwas_regenie/sig_snp_ld/")
ld_threshold = 0.2
# make a dictionary and read LD file for each chromosome
# the LD file was pre-calculated using plink --r2 with the filtered SNPs
ld_dict = {}
for chrom in range(1, 23):
    if chrom == 23:
        chrom = "X"
    ld_f = os.path.join(ld_dir, "ld." + stem + ".filtered_snps.chr" + str(chrom) + ".vcor")
    if os.path.exists(ld_f):
        ld_dict[str(chrom)] = pd.read_csv(ld_f, sep="\t")
        # print status
        print("Read LD file for chromosome " + str(chrom))

# iterate through filtered SNP files, which filtered GWAS results by p-value, MAF, and INFO score
# for each SNP file, perform LD clumping using the pre-calculated LD between filtered SNPs

# list filtered snp files
filter_snp_f_pattern = "filtered." + stem
filter_snp_fs = [os.path.join(sig_snp_dir, f) for f in os.listdir(sig_snp_dir) if f.startswith(filter_snp_f_pattern)]
print("number of files to perform clumping: " + str(len(filter_snp_fs)) + "\n")

# iterate through filtered snp files
for filter_snp_f in filter_snp_fs:
    print("start clumping for file: " + filter_snp_f)
    # read the file as df, and extract the SNP IDs
    filter_snp_df = pd.read_csv(filter_snp_f, sep=" ", comment="#")
    print("number of snps in the file: " + str(filter_snp_df.shape[0]))
    # start LD clumping
    # record index snps in a new df
    index_snps_df = pd.DataFrame(columns=filter_snp_df.columns)
    while filter_snp_df.shape[0] > 0:
        # record the index snp
        min_p_row = filter_snp_df.loc[[filter_snp_df["LOG10P"].idxmax()]]
        index_snps_df = pd.concat([index_snps_df, min_p_row])

        # check if there is any LD to be considered
        if str(min_p_row["CHROM"].iloc[0]) not in ld_dict.keys():
            # LD file does not exists for the chromosome of the index snp
            filter_snp_df = filter_snp_df.loc[filter_snp_df["ID"] != min_p_row["ID"].iloc[0]]
            continue
        # get list of snps that fits the LD threshold
        ld_df = ld_dict[str(min_p_row["CHROM"].iloc[0])]
        # extract relevant rows from ID_A and ID_B
        ld_df = ld_df.loc[(ld_df["ID_A"] == min_p_row["ID"].iloc[0]) | (ld_df["ID_B"] == min_p_row["ID"].iloc[0])]
        # check if there is any SNPs in LD with the index SNP
        if ld_df.shape[0] == 0:
            # if not, remove the index snp
            filter_snp_df = filter_snp_df.loc[filter_snp_df["ID"] != min_p_row["ID"].iloc[0]]
            continue
        # after the check, get the snps that fits the LD threshold
        ld_df = ld_df.loc[ld_df["PHASED_R2"] >= ld_threshold]
        ld_snps = ld_df["ID_A"].tolist() + ld_df["ID_B"].tolist()
        ld_snps = list(set(ld_snps))
        ld_snps.append(min_p_row["ID"].iloc[0])
        # remove the index snp and the snps that fits the LD threshold
        filter_snp_df = filter_snp_df.loc[~filter_snp_df["ID"].isin(ld_snps)]
        # print number of rows left
        print("number of snps left: " + str(filter_snp_df.shape[0]))
    # save the index snps
    logging.info("number of index snps: " + str(index_snps_df.shape[0]))
    output_filename = filter_snp_f.split("/")[-1].replace("filtered", "clump")
    index_snps_df.to_csv(ld_dir + output_filename, sep="\t", index=False)
