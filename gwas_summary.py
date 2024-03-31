#!/usr/bin/env python
# coding: utf-8

import os
import glob
import numpy as np
import pandas as pd

gwas_dir = os.path.expanduser("~/saige_pipe_test/")
data_dir = os.path.expanduser("~/bugbank_data/postgwas_regenie/")

# list GWAS summary files
pattern = os.path.join(gwas_dir, "summary.05062023*txt.gz")
gwas_files = glob.glob(pattern)

# split the file name to get the metadata for each GWAS
def break_down_file_name(file_name):
    file_meta = file_name.split("/")[-1].split(".")[1:5]
    # replace "_" with " " in the pathogen name
    file_meta[1] = file_meta[1].replace("_", " ")
    if "hes" in file_meta[0]:
        file_meta.append("HES")
        file_meta[3] = None
    elif "sgss" in file_meta[0]:
        file_meta.append("SGSS")
    file_meta.append(file_name)
    return file_meta

def get_hit_file(stem, pathogen, tax_lev, specimen):
    hit_file_dir = data_dir + "sig_snp_ld/"
    if specimen != None:
        string_lst = ["clump", stem, pathogen.replace(" ", "_"), tax_lev, specimen, "snp"]
    else:
        string_lst = ["clump", stem, pathogen.replace(" ", "_"), tax_lev, "snp"]
    hit_file = ".".join(string_lst)
    return hit_file_dir + hit_file

def get_hit_tb(stem, pathogen, tax_lev, specimen):
    hit_file = get_hit_file(stem, pathogen, tax_lev, specimen)
    hit_tb = pd.read_csv(hit_file, sep="\t")
    return hit_tb

def get_pheno_file(stem, pathogen, tax_lev, specimen):
    pheno_file_dir = gwas_dir + "ukb41482.bd.gwas-pheno."
    if specimen != None:
        string_lst = [stem, pathogen.replace(" ", "_"), tax_lev, specimen, "txt.gz"]
    else:
        string_lst = [stem, pathogen.replace(" ", "_"), tax_lev, "txt.gz"]
    pheno_file = pheno_file_dir + ".".join(string_lst)
    return pheno_file

def get_case_n(stem, pathogen, tax_lev, specimen):
    pheno_file = get_pheno_file(stem, pathogen, tax_lev, specimen)
    pheno_tb = pd.read_csv(pheno_file, sep="\t")
    case_n = pheno_tb["pheno"].value_counts().get(1, 0)
    # get the 
    return case_n


files_meta = [break_down_file_name(file_name) for file_name in gwas_files]
# make the list into df
files_meta = pd.DataFrame(files_meta, columns=["stem", "pathogen", "tax_lev", "specimen", "source", "file_name"])
# get the hit number for each GWAS
files_meta["hit_num"] = files_meta.apply(lambda x: get_hit_tb(x.stem, x.pathogen, x.tax_lev, x.specimen).shape[0], axis=1)
files_meta["case_n"] = files_meta.apply(lambda x: get_case_n(x.stem, x.pathogen, x.tax_lev, x.specimen), axis=1)
files_meta.to_csv(data_dir + "gwas_summary_meta.csv", index=False)


# In[71]:


# find the GWAS hits for pathogen presented in both SGSS and HES
# get the pathogen list
pathogen_sgss = files_meta[files_meta.source == "SGSS"].pathogen.unique().tolist()
pathogen_hes = files_meta[files_meta.source == "HES"].pathogen.unique().tolist()
pathogen_both = list(set(pathogen_sgss).intersection(set(pathogen_hes)))

# get the GWAS summary files for pathogen presented in both SGSS and HES
files_meta_both_sgss = files_meta[files_meta.pathogen.isin(pathogen_both) & (files_meta.source == "SGSS") & (files_meta.specimen == "all")]
files_meta_both_hes = files_meta[files_meta.pathogen.isin(pathogen_both) & (files_meta.source == "HES")]
# join the two df
files_meta_both = pd.merge(files_meta_both_sgss, files_meta_both_hes, on=["pathogen", "tax_lev"], suffixes=("_sgss", "_hes"))
files_meta_both[["pathogen", "case_n_sgss", "case_n_hes", "hit_num_sgss", "hit_num_hes"]].sort_values(by = "case_n_sgss", ascending=False)

                                


# In[63]:


# find the GWAS hits for pathogen presented in SGSS only
pathogen_sgss = list(set(pathogen_sgss) - set(pathogen_both))
files_meta_sgss = files_meta[files_meta.pathogen.isin(pathogen_sgss) & (files_meta.source == "SGSS") & (files_meta.specimen == "all")]
files_meta_sgss[["pathogen", "hit_num"]].sort_values(by="hit_num", ascending=False)


# In[64]:


# find the GWAS hits for pathogen presented in HES only
pathogen_hes = list(set(pathogen_hes) - set(pathogen_both))
files_meta_hes = files_meta[files_meta.pathogen.isin(pathogen_hes) & (files_meta.source == "HES")]
files_meta_hes[["pathogen", "hit_num"]].sort_values(by="hit_num", ascending=False)


# In[ ]:




