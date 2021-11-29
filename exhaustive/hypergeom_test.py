import numpy as np
import pandas as pd
import scipy.stats

import sys
import time
import tqdm
import json

# Import python package
import core.extern

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
INTERACTION_PATH = config["interaction_path"]
OUTPUT_DIR_PATH = config["output_dir_path"]

REFERENCE_GROUP = config["reference_group"]
EXPERIMENTAL_GROUP = config["experimnetal_group"]

CORRELATION = config["correlation"]
ALTERNATIVE = config["alternative"]
PROCESS_NUMBER = config["process_number"]

FDR_THRESHOLD = config["fdr_treshold"]

def get_occurrence(array):
    occurrence = {}
    for elem in array:
        if elem in occurrence:
            occurrence[elem] += 1
        else:
            occurrence[elem] = 1
    
    return occurrence


# Main part
report_stat = np.load(OUTPUT_DIR_PATH.rstrip("/") + \
    "/{}_report_stat.npy".format(CORRELATION)
)
report_pvalue = np.load(OUTPUT_DIR_PATH.rstrip("/") + \
    "/{}_report_pvalue.npy".format(CORRELATION)
)
report_fdr = np.load(OUTPUT_DIR_PATH.rstrip("/") + \
    "/{}_report_fdr.npy".format(CORRELATION)
)
molecules = pd.read_csv(DATA_PATH, sep=",", index_col=0).index.to_numpy()

indexes = np.where(report_fdr < FDR_THRESHOLD)[0]
report_stat = report_stat[indexes]
report_pvalue = report_pvalue[indexes]
report_fdr = report_fdr[indexes]
source_indexes = []
target_indexes = []
for ind in tqdm.tqdm(indexes):
    s, t = core.extern.paired_index(ind, len(molecules))
    source_indexes.append(molecules[s])
    target_indexes.append(molecules[t])
source_indexes = np.array(source_indexes, dtype=np.str)
target_indexes = np.array(source_indexes, dtype=np.str)

# Report indexing
ll = list(source_indexes)
ll.extend(list(target_indexes))
report_occurrence = get_occurrence(ll)

report_interaction_number = 0
for index in report_occurrence:
    report_interaction_number += report_occurrence[index]

initial_interaction_number = len(molecules) * (len(molecules) - 1) // 2

# Analisis
output_df = pd.DataFrame()
output_df["Molecule"] = [index for index in report_occurrence]
output_df["Diff"] = [report_occurrence[molecule] for molecule in output_df["Molecule"]]
output_df["Total"] = [len(molecules) - 1 for molecule in output_df["Molecule"]]
output_df["Proportion"] = output_df["Diff"] / output_df["Total"]

output_df["Pvalue"] = 1 - scipy.stats.hypergeom.cdf(
    output_df["Diff"] - 1,
    initial_interaction_number, # M
    report_interaction_number, # n
    output_df["Total"] # N
)

adjusted_pvalue = np.array(output_df["Pvalue"] * len(output_df["Pvalue"]) / \
    scipy.stats.rankdata(output_df["Pvalue"]))
adjusted_pvalue[adjusted_pvalue > 1] = 1
adjusted_pvalue = adjusted_pvalue.flatten()
output_df["FDR"] = adjusted_pvalue
output_df = output_df[output_df["FDR"] < FDR_THRESHOLD]

output_df["NegProportion"] = 1 - output_df["Proportion"]
output_df = output_df.sort_values(["FDR", "Pvalue", "NegProportion"])
output_df = output_df.drop(columns=["NegProportion"])

output_df.to_csv(
    OUTPUT_DIR_PATH.rstrip("/") +
    "/{}_hyper.csv".format(CORRELATION),
    sep=",",
    index=None
)
