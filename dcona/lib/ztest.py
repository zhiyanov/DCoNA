import numpy as np
import pandas as pd
import scipy.stats
from multiprocessing import cpu_count

from ..core import extern
from . import utils
from . import dump


def ztest(
    data_df, description_df,
    reference_group, experimental_group,
    correlation="spearman", alternative="two-sided",
    interaction=None,
    repeats_number=None,
    output_dir=None,
    process_number=None
):
    if process_number is None:
        process_number = cpu_count()

    # If gene names are in dataframe column, relocate them to df.index
    if not pd.api.types.is_number(data_df.iloc[0, 0]):
        data_df = data_df.copy()
        data_df.set_index(data_df.columns[0], inplace=True)
        
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
            repeats_number = 0
    
    sorted_indexes, df_indexes, \
    ref_corrs, ref_pvalues, exp_corrs, exp_pvalues, \
    stat, pvalue, adjusted_pvalue, \
    boot_pvalue = _ztest(
        data_df, description_df, interaction_df, \
        reference_group, experimental_group, \
        correlation, alternative, \
        repeats_number, process_number
    )
    
    if output_dir:
        dump.check_directory_existence(output_dir)
        
        if repeats_number > 0: 
            df_template = pd.DataFrame(columns=[
                "Source", "Target", "RefCorr", "RefPvalue", 
                "ExpCorr", "ExpPvalue", "Statistic",
                "Pvalue", "AdjPvalue", "PermutePvalue"
            ])

            df_columns = [
                ref_corrs, ref_pvalues,
                exp_corrs, exp_pvalues, stat,
                pvalue, adjusted_pvalue, boot_pvalue
            ]

        else:
            df_template = pd.DataFrame(columns=[
                "Source", "Target", "RefCorr", "RefPvalue", 
                "ExpCorr", "ExpPvalue", "Statistic",
                "Pvalue", "AdjPvalue"
            ])

            df_columns = [
                ref_corrs, ref_pvalues,
                exp_corrs, exp_pvalues, stat,
                pvalue, adjusted_pvalue
            ]


        path_to_file = output_dir.rstrip("/") + f"/{correlation}_{alternative}_ztest.csv"
        dump.save_by_chunks(
            sorted_indexes,
            df_indexes, df_template, df_columns,
            path_to_file
        )
        
        print(f"File saved at: {path_to_file}")
        return None
    
    source_indexes = []
    target_indexes = []

    if (isinstance(df_indexes, tuple)) and (len(df_indexes) == 2):
        source_indexes, target_indexes = df_indexes
    else:
        for ind in sorted_indexes:
            s, t = extern.paired_index(ind, len(df_indexes))
            source_indexes.append(df_indexes[s])
            target_indexes.append(df_indexes[t])
    
    if repeats_number > 0:
        output_df = pd.DataFrame(data={
            "Source": source_indexes,
            "Target": target_indexes,
            "RefCorr": ref_corrs[sorted_indexes], 
            "RefPvalue": ref_pvalues[sorted_indexes], 
            "ExpCorr": exp_corrs[sorted_indexes], 
            "ExpPvalue": exp_pvalues[sorted_indexes], 
            "Statistic": stat[sorted_indexes],
            "Pvalue": pvalue[sorted_indexes], 
            "AdjPvalue": adjusted_pvalue[sorted_indexes],
            "PermutePvalue": boot_pvalue[sorted_indexes] 
        })
    else:
        output_df = pd.DataFrame(data={
            "Source": source_indexes,
            "Target": target_indexes,
            "RefCorr": ref_corrs[sorted_indexes], 
            "RefPvalue": ref_pvalues[sorted_indexes], 
            "ExpCorr": exp_corrs[sorted_indexes], 
            "ExpPvalue": exp_pvalues[sorted_indexes], 
            "Statistic": stat[sorted_indexes],
            "Pvalue": pvalue[sorted_indexes], 
            "AdjPvalue": adjusted_pvalue[sorted_indexes]
        })

    return output_df

def _ztest(
    data_df, description_df, interaction_df,
    reference_group, experimental_group,
    correlation, alternative,
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

    print("Z-test computation")
    ref_corrs, ref_pvalues, \
    exp_corrs, exp_pvalues, \
    stat, pvalue, boot_pvalue = \
    extern.ztest_pipeline(
        data_df,
        reference_indexes,
        experimental_indexes,
        source_indexes,
        target_indexes,
        correlation=correlation,
        alternative=alternative,
        repeats_num=repeats_number,
        process_num=process_number,
        correlation_alternative="two-sided"
    )

    print("Adjusted p-value computation")
    adjusted_pvalue = pvalue * len(pvalue) / \
        scipy.stats.rankdata(pvalue)
    adjusted_pvalue[adjusted_pvalue > 1] = 1
    adjusted_pvalue = adjusted_pvalue.flatten()
    
    if interaction_df is None:
        df_indexes = data_df.index.to_numpy()
    else:
        df_indexes = (source_indexes, target_indexes)

    fdr_pvalue = np.core.records.fromarrays(
        [adjusted_pvalue, pvalue],
        names='fdr, pvalue'
    )

    sorted_indexes = np.argsort(fdr_pvalue, order=('fdr', 'pvalue'))
    del fdr_pvalue

    return sorted_indexes, df_indexes, \
        ref_corrs, ref_pvalues, exp_corrs, exp_pvalues, \
        stat, pvalue, adjusted_pvalue, \
        boot_pvalue
