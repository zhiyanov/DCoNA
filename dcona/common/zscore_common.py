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

sources, scores, pvalues = \
core.extern.score_pipeline(
    data_df,
    reference_indexes,
    experimental_indexes,
    source_indexes,
    target_indexes,
    correlation=CORRELATION,
    score=SCORE,
    alternative=ALTERNATIVE,
    repeats_num=REPEATS_NUMBER,
    process_num=PROCESS_NUMBER
)
print("Computational time: {:.3f}".format(time.time() - start))

adjusted_pvalue = pvalues * len(pvalues) / \
    scipy.stats.rankdata(pvalues)
adjusted_pvalue[adjusted_pvalue > 1] = 1
adjusted_pvalue = adjusted_pvalue.flatten()

# Generate report
print("Report phase")
output_df = pd.DataFrame()
output_df["Source"] = data_df.index.to_numpy()[sources]
output_df["Score"] = scores
output_df["Pvalue"] = pvalues
output_df["FDR"] = adjusted_pvalue
output_df = output_df.sort_values(["FDR", "Pvalue"])

output_df.to_csv(
    OUTPUT_DIR_PATH.rstrip("/") + \
    "/{}_{}_{}_zscore.csv".format(CORRELATION, SCORE, ALTERNATIVE),
    sep=",",
    index=None
)
