import numpy as np
import pandas as pd
import scipy.stats

import sys
import time
import tqdm
import json
import os
import time

# Import python package
import core.extern
import core.utils

# Arg parser
import argparse


# Argument parser
parser = argparse.ArgumentParser()
parser.add_argument("config_path")
args = parser.parse_args()

# Load config file
CONFIG_PATH = args.config_path
config = json.load(open(CONFIG_PATH, "r"))

DATA_PATH = config["data_path"]
DESCRIPTION_PATH = config["description_path"]
# INTERACTION_PATH = config["interaction_path"]
OUTPUT_DIR_PATH = config["output_dir_path"]

REFERENCE_GROUP = config["reference_group"]
EXPERIMENTAL_GROUP = config["experimental_group"]

CORRELATION = config["correlation"]
ALTERNATIVE = config["alternative"]
REPEATS_NUMBER = config["repeats_number"]
PROCESS_NUMBER = config["process_number"]

FDR_THRESHOLD = config["fdr_treshold"]

core.utils.check_directory_existence(OUTPUT_DIR_PATH)

# Main part
data_df = pd.read_csv(DATA_PATH, sep=",", index_col=0)
description_df = pd.read_csv(DESCRIPTION_PATH, sep=",")

if (CORRELATION == "spearman"):
    correlation = core.extern.spearmanr
elif (CORRELATION == "pearson"):
    correlation = core.extern.pearsonr
elif (CORRELATION == "spearman_test"):
    CORRELATION = "spearman"
    correlation = core.extern.spearmanr_test

interaction_df=None
source_indexes=None
target_indexes=None

reference_indexes = description_df.loc[
    description_df["Group"] == REFERENCE_GROUP,
    "Sample"
].to_list()
experimental_indexes = description_df.loc[
    description_df["Group"] == EXPERIMENTAL_GROUP,
    "Sample"
].to_list()

# Test mode
# data_df = data_df.iloc[:1000]

print("Pipeline")
start = time.time()

ref_corrs, ref_pvalues, \
exp_corrs, exp_pvalues, \
stat, pvalue, boot_pvalue = \
core.extern.ztest_pipeline(
    data_df,
    reference_indexes,
    experimental_indexes,
    source_indexes,
    target_indexes,
    correlation=CORRELATION,
    alternative=ALTERNATIVE,
    repeats_num=REPEATS_NUMBER,
    process_num=PROCESS_NUMBER,
    correlation_alternative="two-sided"
)
print("Computational time: {:.3f}".format(time.time() - start))

adjusted_pvalue = pvalue * len(pvalue) / \
    scipy.stats.rankdata(pvalue)
adjusted_pvalue[adjusted_pvalue > 1] = 1
adjusted_pvalue = adjusted_pvalue.flatten()

indexes = np.arange(len(adjusted_pvalue))

ref_corrs = ref_corrs[indexes]
ref_pvalues = ref_pvalues[indexes]

exp_corrs = exp_corrs[indexes]
exp_pvalues = exp_pvalues[indexes]

stat = stat[indexes]
pvalue = pvalue[indexes]
boot_pvalue = boot_pvalue[indexes]

adjusted_pvalue = adjusted_pvalue[indexes]
df_indexes = data_df.index.to_numpy()

np.save(OUTPUT_DIR_PATH.rstrip("/") + \
    "/{}_ztest_pipeline_stat.npy".format(CORRELATION),
    stat
)
np.save(OUTPUT_DIR_PATH.rstrip("/") + \
    "/{}_ztest_pipeline_pvalue.npy".format(CORRELATION),
    pvalue
)
np.save(OUTPUT_DIR_PATH.rstrip("/") + \
    "/{}_ztest_pipeline_bootpv.npy".format(CORRELATION),
    boot_pvalue
)
np.save(OUTPUT_DIR_PATH.rstrip("/") + \
    "/{}_ztest_pipeline_fdr.npy".format(CORRELATION),
    adjusted_pvalue
)

# Test mode
# indexes = np.where(adjusted_pvalue < FDR_THRESHOLD)[0]

ref_corrs = ref_corrs[indexes]
ref_pvalues = ref_pvalues[indexes]

exp_corrs = exp_corrs[indexes]
exp_pvalues = exp_pvalues[indexes]

stat = stat[indexes]
pvalue = pvalue[indexes]
boot_pvalue = boot_pvalue[indexes]
adjusted_pvalue = adjusted_pvalue[indexes]
df_indexes = data_df.index.to_numpy()


# Save mode
print("Creating FDR/pvalue array")
FDR_pvalue = np.core.records.fromarrays(
    [adjusted_pvalue, pvalue],
    names='FDR, pvalue'
)

print("Sorting FDR/pvalue array")
sorted_indexes = np.argsort(FDR_pvalue, order=('FDR','pvalue'))
del FDR_pvalue
print("FDR/pvalue array is sorted")

df_template = pd.DataFrame(columns=[
    "Source", "Target", "RefCorr", "RefPvalue", 
    "ExpCorr", "ExpPvalue", "Statistic",
    "Pvalue", "Bootpv", "FDR"
])
df_columns = [
    ref_corrs, ref_pvalues,
    exp_corrs, exp_pvalues, stat,
    pvalue, boot_pvalue, adjusted_pvalue
]

path_to_file = OUTPUT_DIR_PATH.rstrip("/") + "/{}_ztest_pipeline.csv".format(CORRELATION)
core.utils.save_by_chunks(
    sorted_indexes, 
    df_indexes, df_template, df_columns,
    path_to_file
)
