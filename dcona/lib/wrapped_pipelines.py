import numpy as np
import pandas as pd
import time
import scipy.stats
from multiprocessing import cpu_count

# Import python package
import dcona.core.extern
import dcona.core.utils


def wrapped_ztest(
                  data_df, description_df, interaction_df, \
                  REFERENCE_GROUP, EXPERIMENTAL_GROUP, CORRELATION, \
                  ALTERNATIVE, SCORE, REPEATS_NUMBER, PROCESS_NUMBER
):

    if (CORRELATION == "spearman"):
        correlation = dcona.core.extern.spearmanr
    elif (CORRELATION == "pearson"):
        correlation = dcona.core.extern.pearsonr
    elif (CORRELATION == "spearman_test"):
        CORRELATION = "spearman"
        correlation = dcona.core.extern.spearmanr_test

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
    dcona.core.extern.ztest_pipeline(
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

    # Filtering results
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
    print("FDR/pvalue array is sorted", flush=True)

    return sorted_indexes, indexes, df_indexes, ref_corrs, \
           ref_pvalues, exp_corrs, exp_pvalues, stat, pvalue, \
           boot_pvalue, adjusted_pvalue


def wrapped_zscore(
                  data_df, description_df, interaction_df, \
                  REFERENCE_GROUP, EXPERIMENTAL_GROUP, CORRELATION, \
                  ALTERNATIVE, SCORE, REPEATS_NUMBER, PROCESS_NUMBER
):

    if (CORRELATION == "spearman"):
        correlation = dcona.core.extern.spearmanr
    elif (CORRELATION == "pearson"):
        correlation = dcona.core.extern.pearsonr
    elif (CORRELATION == "spearman_test"):
        CORRELATION = "spearman"
        correlation = dcona.core.extern.spearmanr_test

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
    dcona.core.extern.score_pipeline(
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
    
    return data_df, sources, scores, pvalues, adjusted_pvalue


def zscore_to_df(
                 data_df, sources, scores, pvalues, adjusted_pvalue
):
    output_df = pd.DataFrame()
    output_df["Source"] = data_df.index.to_numpy()[sources]
    output_df["Score"] = scores
    output_df["Pvalue"] = pvalues
    output_df["FDR"] = adjusted_pvalue
    output_df = output_df.sort_values(["FDR", "Pvalue"])
    
    return output_df


def ztest_to_df(
                sorted_indexes, indexes, df_indexes,
                ref_corrs, ref_pvalues,
                exp_corrs, exp_pvalues, stat,
                pvalue, boot_pvalue, adjusted_pvalue        
):
    source_indexes = []
    target_indexes = []

    for ind in sorted_indexes:
        s, t = dcona.core.extern.paired_index(ind, len(df_indexes))
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
                                    "FDR": adjusted_pvalue[sorted_indexes]
                                    })
    return output_df


def ztest_to_file(
                  sorted_indexes, indexes, df_indexes,
                  ref_corrs, ref_pvalues,
                  exp_corrs, exp_pvalues, stat,
                  pvalue, boot_pvalue, adjusted_pvalue,
                  CORRELATION
):
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

    path_to_file = OUTPUT_DIR_PATH.rstrip("/") + "/{}_ztest.csv".format(CORRELATION)
    dcona.core.utils.save_by_chunks(
        sorted_indexes, 
        df_indexes, df_template, df_columns,
        path_to_file,
        # index_transform=None
        index_transform=indexes
    )


def zscore(data_df, description_df, \
           reference_group, experimental_group, \
           output_dir=None, interaction=None, \
           correlation="spearman", alternative="two-sided", score="mean", \
           repeats=1000, process_number=cpu_count()
):

    # If gene names are in dataframe column, relocate them in df.index
    if not pd.api.types.is_number(data_df.iloc[0, 0]):
        data_df = data_df.copy()
        data_df.set_index(data_df.columns[0], inplace = True)
        
    if interaction is not None:
        if isinstance(interaction, pd.core.frame.DataFrame):
            interaction_df = interaction
        else:
            interaction_df = dcona.core.utils.form_gene_pairs(interaction)
    else:
        interaction_df = None

    result = wrapped_zscore(data_df, description_df, interaction_df, \
                            reference_group, experimental_group, correlation, \
                            alternative, score, repeats, process_number)
                            
    output_df = zscore_to_df(*result)
    
    if output_dir:
        dcona.core.utils.check_directory_existence(output_dir)
        path_to_file = output_dir.rstrip("/") + \
                       "/{}_{}_{}_zscore.csv".format(correlation, score, alternative)
        output_df.to_csv(
            path_to_file,
            sep=",",
            index=None
        )
        print(f"File saved at: {path_to_file}")
        
        return
    else:
        return output_df


def ztest(data_df, description_df, \
          reference_group, experimental_group, \
          output_dir=None, interaction=None, \
          correlation="spearman", alternative="two-sided", score="mean", \
          repeats=1000, process_number=cpu_count()
):

    # If gene names are in dataframe column, relocate them in df.index
    if not pd.api.types.is_number(data_df.iloc[0, 0]):
        data_df = data_df.copy()
        data_df.set_index(data_df.columns[0], inplace = True)
        
        
    if interaction is not None:
        if isinstance(interaction, pd.core.frame.DataFrame):
            interaction_df = interaction
        else:
            interaction_df = dcona.core.utils.form_gene_pairs(interaction)
    else:
        interaction_df = None

    sorted_indexes, indexes, \
    df_indexes, ref_corrs, \
    ref_pvalues, exp_corrs, \
    exp_pvalues, stat, pvalue, \
    boot_pvalue, adjusted_pvalue = \
    wrapped_ztest(data_df, description_df, interaction_df, \
                  reference_group, experimental_group, correlation, \
                  alternative, score, repeats, process_number)
    

    if output_dir:
        dcona.core.utils.check_directory_existence(output_dir)
        
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

        path_to_file = output_dir.rstrip("/") + "/{}_ztest.csv".format(correlation)
        dcona.core.utils.save_by_chunks(
            sorted_indexes, 
            df_indexes, df_template, df_columns,
            path_to_file,
            # index_transform=None
            index_transform=indexes
        )
        print(f"File saved at: {path_to_file}")
        return
    else:
        source_indexes = []
        target_indexes = []

        for ind in sorted_indexes:
            s, t = dcona.core.extern.paired_index(ind, len(df_indexes))
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
                                        "FDR": adjusted_pvalue[sorted_indexes]
                                        })
        return output_df
