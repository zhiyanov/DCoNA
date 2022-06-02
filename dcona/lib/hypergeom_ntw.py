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
parser.add_argument("-o", "--oriented", action="store_true")
args = parser.parse_args()

# Load config file
CONFIG_PATH = args.config_path
config = json.load(open(CONFIG_PATH, "r"))

DATA_PATH = config["data_path"]
DESCRIPTION_PATH = config["description_path"]
INTERACTION_PATH = config["interaction_path"]
OUTPUT_DIR_PATH = config["output_dir_path"]

REFERENCE_GROUP = config["reference_group"]
EXPERIMENTAL_GROUP = config["experimental_group"]

CORRELATION = config["correlation"]
ALTERNATIVE = config["alternative"]
PROCESS_NUMBER = config["process_number"]

FDR_THRESHOLD = config["fdr_treshold"]

ORIENTED = args.oriented

def get_occurrence(array):
    occurrence = {}
    for elem in array:
        if elem in occurrence:
            occurrence[elem] += 1
        else:
            occurrence[elem] = 1
    
    return occurrence


# Main part
data_df = pd.read_csv(DATA_PATH, sep=",", index_col=0)

report_df = pd.read_csv(
    OUTPUT_DIR_PATH.rstrip("/") + 
    # "/{}_report.csv".format(CORRELATION),
    "/{}_ztest.csv".format(CORRELATION),
    sep=","
)
if ALTERNATIVE == "less":
    report_df = report_df[
        (report_df["FDR"] < FDR_THRESHOLD) &
        (report_df["Statistic"] < 0)
    ]
elif ALTERNATIVE == "greater":
    report_df = report_df[
        (report_df["FDR"] < FDR_THRESHOLD) &
        (report_df["Statistic"] > 0)
    ]
else:
    # ALTERNATIVE == "two-sided"
    report_df = report_df[report_df["FDR"] < FDR_THRESHOLD]

# Report indexing
if ORIENTED:
    report_occurrence = get_occurrence(list(report_df["Source"]))
else:
    ll = list(report_df["Source"])
    ll.extend(list(report_df["Target"]))
    report_occurrence = get_occurrence(ll)

report_interaction_number = 0
for index in report_occurrence:
    report_interaction_number += report_occurrence[index]

# Description indexing
interaction_df = pd.read_csv(INTERACTION_PATH, sep=",")

data_molecules = set(data_df.index.to_list())
interaction_df = interaction_df[interaction_df["Source"].isin(data_molecules)]
interaction_df = interaction_df[interaction_df["Target"].isin(data_molecules)]

if ORIENTED:
    initial_occurrence = get_occurrence(list(interaction_df["Source"]))
else:
    ll = list(interaction_df["Source"])
    ll.extend(list(interaction_df["Target"]))
    initial_occurrence = get_occurrence(ll)

initial_interaction_number = 0
for index in initial_occurrence:
    initial_interaction_number += initial_occurrence[index]

# Analisis
output_df = pd.DataFrame()
output_df["Molecule"] = [index for index in report_occurrence]
output_df["Diff"] = [report_occurrence[molecule] for molecule in output_df["Molecule"]]
output_df["Total"] = [initial_occurrence[molecule] for molecule in output_df["Molecule"]]
output_df["Proportion"] = output_df["Diff"] / output_df["Total"]
output_df["Pvalue"] = 1 - scipy.stats.hypergeom.cdf(
    output_df["Diff"] - 1,
    initial_interaction_number,
    report_interaction_number,
    output_df["Total"]
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
    "/{}_{}_hypergeom.csv".format(CORRELATION, ALTERNATIVE),
    sep=",",
    index=None
)
