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
config = json.load(open(CONFIG_PATH, "r"))

DATA_PATH = config["data_path"]
DESCRIPTION_PATH = config["description_path"]
# INTERACTION_PATH = config["interaction_path"]
# OUTPUT_DIR_PATH = config["output_dir_path"]

REFERENCE_GROUP = config["reference_group"]
EXPERIMENTAL_GROUP = config["experimental_group"]

CORRELATION = config["correlation"]
ALTERNATIVE = config["alternative"]
SCORE = config["score"]

# REPEATS_NUMBER = config["repeats_number"]
REPEATS_NUMBER = 1
PROCESS_NUMBER = config["process_number"]

FDR_THRESHOLD = config["fdr_treshold"]

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
# data_df = data_df.iloc[:100]

print("Pipeline")
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
end = time.time()
print("Computational time: {:.3f}".format(end - start))

# Expected computation time
def score_comlexity(size):
    if SCORE == "mean":
        return size
    else:
        return size * np.log(size)

def correlation_complexity(sample_size):
    if CORRELATION == "spearman":
        return sample_size
    else:
        return sample_size * np.log(sample_size)

warmup_computation_complexity = \
    correlation_complexity(len(reference_indexes)) * len(data_df)**2 // 2 + \
    correlation_complexity(len(experimental_indexes)) * len(data_df)**2 // 2 + \
    len(data_df) * score_comlexity(len(data_df))

computation_time = []
for i in tqdm.tqdm(range(100, len(data_df), 100)):
    pairs_number = i**2 // 2
    computation_complexity = \
        correlation_complexity(len(reference_indexes)) * pairs_number + \
        correlation_complexity(len(experimental_indexes)) * pairs_number + \
        i * score_comlexity(i)
    
    recommended_repeats_number = i / FDR_THRESHOLD
    
    time_complexity = (end - start) * \
        computation_complexity / warmup_computation_complexity * \
        recommended_repeats_number / REPEATS_NUMBER / 360

    computation_time.append(time_complexity)

print(recommended_repeats_number)

import matplotlib.pyplot as plt

computation_time = np.array(computation_time)
computation_time = computation_time[computation_time < 24]
plt.plot(
    list(range(100, len(data_df), 100))[:len(computation_time)],
    computation_time
)
# plt.yscale("log")
# plt.xscale("log")

plt.ylabel("Time (hours)")
plt.xlabel("Number of genes")
plt.savefig("warmup.png")



