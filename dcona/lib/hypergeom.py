import numpy as np
import pandas as pd
import scipy.stats

from . import dump

FDR_THRESHOLD = 0.05


def hypergeom(
    ztest_df,
    alternative="two-sided",
    oriented=True,
    output_dir=None
):
    if alternative == "less":
        report_df = ztest_df[
            (ztest_df["AdjPvalue"] < FDR_THRESHOLD) &
            (ztest_df["Statistic"] < 0)
        ]
    elif alternative == "greater":
        report_df = ztest_df[
            (ztest_df["AdjPvalue"] < FDR_THRESHOLD) &
            (ztest_df["Statistic"] > 0)
        ]
    else:
        report_df = ztest_df[ztest_df["AdjPvalue"] < FDR_THRESHOLD]
    
    if oriented:
        report_occurrence = get_occurrence(list(report_df["Source"]))
        initial_occurrence = get_occurrence(list(ztest_df["Source"]))
    else:
        ll = list(report_df["Source"])
        ll.extend(list(report_df["Target"]))
        report_occurrence = get_occurrence(ll)

        ll = list(data_df["Source"])
        ll.extend(list(data_df["Target"]))
        initial_occurrence = get_occurrence(ll)

    report_interaction_number = 0
    for index in report_occurrence:
        report_interaction_number += report_occurrence[index]

    initial_interaction_number = 0
    for index in initial_occurrence:
        initial_interaction_number += initial_occurrence[index]

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
    output_df["AdjPvalue"] = adjusted_pvalue
    
    if output_dir:
        dump.check_directory_existence(output_dir)
        path_to_file = output_dir.rstrip("/") + \
            f"/{alternative}_hypergeom.csv"
        output_df.to_csv(
            path_to_file,
            sep=",",
            index=None
        )

        print(f"File saved at: {path_to_file}")
        return None
    else:
        return output_df

    return output_df

def get_occurrence(array):
    occurrence = {}
    for elem in array:
        if elem in occurrence:
            occurrence[elem] += 1
        else:
            occurrence[elem] = 1
    
    return occurrence
