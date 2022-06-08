import numpy as np
import pandas as pd
import scipy.stats

from ..core import extern
from . import utils
from . import dump


def zscore(
    data_df, description_df,
    reference_group, experimental_group,
    correlation="spearman", score="mean", alternative="two-sided",
    interaction=None,
    repeats_number=None,
    output_dir=None,
    process_number=None 
):
    if process_number is None:
        process_number = 1

    # if gene names are in dataframe column, relocate them in df.index
    # if not pd.api.types.is_number(data_df.iloc[0, 0]):
    #     data_df = data_df.copy()
    #     data_df.set_index(data_df.columns[0], inplace=True)
        
    if interaction is not None:
        if isinstance(interaction, pd.core.frame.DataFrame):
            interaction_df = interaction
        else:
            interaction_df = utils.generate_pairs(interaction)
        
        if repeats_number is None:
            repeats_number = int(len(interaction_df) / 0.05)
    else:
        interaction_df = None
        
        if repeats_number is None:
            repeats_number = int(len(data_df) / 0.05)

    data_df, sources, scores, \
    pvalues, adjusted_pvalue = \
    _zscore(
        data_df, description_df, interaction_df, \
        reference_group, experimental_group, \
        correlation, score, alternative, \
        repeats_number, process_number
    )

    output_df = pd.DataFrame(data={
        "Source": data_df.index.to_numpy()[sources],
        "Score": scores,
        "Pvalue": pvalues, 
        "AdjPvalue": adjusted_pvalue, 
    })
    output_df = output_df.sort_values(["AdjPvalue", "Pvalue"])
                            
    if output_dir:
        dump.check_directory_existence(output_dir)
        path_to_file = output_dir.rstrip("/") + \
            f"/{correlation}_{score}_{alternative}_zscore.csv"
        output_df.to_csv(
            path_to_file,
            sep=",",
            index=None
        )

        print(f"File saved at: {path_to_file}")
        return None
    else:
        return output_df

def _zscore(
    data_df, description_df, interaction_df,
    reference_group, experimental_group,
    correlation, score, alternative,
    repeats_number, process_number
):
    if (correlation != "spearman"):
        correlation = "pearson"

    if interaction_df is not None:
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
        interaction_df = None
        source_indexes = None
        target_indexes = None 

    reference_indexes = description_df.loc[
        description_df["Group"] == reference_group,
        "Sample"
    ].to_list()
    experimental_indexes = description_df.loc[
        description_df["Group"] == experimental_group,
        "Sample"
    ].to_list()

    print("Z-score computation")
    sources, scores, pvalues = \
    extern.score_pipeline(
        data_df,
        reference_indexes,
        experimental_indexes,
        source_indexes,
        target_indexes,
        correlation=correlation,
        score=score,
        alternative=alternative,
        repeats_num=repeats_number,
        process_num=process_number
    )

    print("Adjusted p-value computation")
    adjusted_pvalue = pvalues * len(pvalues) / \
        scipy.stats.rankdata(pvalues)
    adjusted_pvalue[adjusted_pvalue > 1] = 1
    adjusted_pvalue = adjusted_pvalue.flatten()
    
    return data_df, sources, scores, \
        pvalues, adjusted_pvalue

