import numpy as np
import pandas as pd
import scipy.stats

import sys
import time
import tqdm
import json

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

DATA_PATH, DESCRIPTION_PATH, OUTPUT_DIR_PATH, INTERACTION_PATH, \
REFERENCE_GROUP, EXPERIMENTAL_GROUP, CORRELATION, \
ALTERNATIVE, REPEATS_NUMBER, PROCESS_NUMBER, FDR_THRESHOLD = \
core.utils.read_json(CONFIG_PATH)

core.utils.check_directory_existence(OUTPUT_DIR_PATH)

### ------------------------------------ ###

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

if INTERACTION_PATH:
    interaction_df = pd.read_csv(INTERACTION_PATH, sep=",")

    data_molecules = set(data_df.index.to_list())
    interaction_df = interaction_df[interaction_df["Source"].isin(data_molecules)]
    interaction_df = interaction_df[interaction_df["Target"].isin(data_molecules)]

    source_indexes = interaction_df["Source"]
    target_indexes = interaction_df["Target"]

    data_molecules = data_molecules.intersection(
        set(source_indexes) | set(target_indexes)
    )
    data_df = data_df.loc[data_df.index.isin(data_molecules)]
else:
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

print("Bootstrap phase")
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

# Remove low fdr interactions 
indexes = np.where(adjusted_pvalue < FDR_THRESHOLD)[0]

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

source_indexes = []
target_indexes = []

for ind in sorted_indexes:
    s, t = core.extern.paired_index(ind, len(df_indexes))
    source_indexes.append(df_indexes[s])
    target_indexes.append(df_indexes[t])


output_df = pd.DataFrame(data = {
                                "Source": source_indexes,
                                "Target": target_indexes,
                                "RefCorr": ref_corrs[sorted_indexes], 
                                "RefPvalue": ref_pvalues[sorted_indexes], 
                                "ExpCorr": exp_corrs[sorted_indexes], 
                                "ExpPvalue": exp_pvalues[sorted_indexes], 
                                "Statistic": stat[sorted_indexes],
                                "Pvalue": pvalue[sorted_indexes], 
                                "Bootpv": boot_pvalue[sorted_indexes], 
                                "FDR": adjusted_pvalue[sorted_indexes]})



'''
# Generate report
print("Report phase")
output_df = pd.DataFrame(interaction_df)

output_df["RefCorr"] = ref_corrs 
output_df["RefPvalue"] = ref_pvalues 

output_df["ExpCorr"] = exp_corrs 
output_df["ExpPvalues"] = exp_pvalues

output_df["Statistic"] = stat
output_df["Pvalue"] = pvalue
output_df["Bootpv"] = boot_pvalue
output_df["FDR"] = adjusted_pvalue
output_df = output_df.sort_values(["FDR", "Pvalue"])
'''
output_df.to_csv(
    OUTPUT_DIR_PATH.rstrip("/") + \
    "/{}_ztest.csv".format(CORRELATION),
    sep=",",
    index=None
)
