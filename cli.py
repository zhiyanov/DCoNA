import pandas as pd
import argparse

import dcona.core.extern
import dcona.core.utils
import dcona.common.wrapped_pipelines as pipelines


def ztest_cli(config_path):
    DATA_PATH, DESCRIPTION_PATH, OUTPUT_DIR_PATH, INTERACTION_PATH, \
    reference_group, experimental_group, correlation, \
    alternative, score, repeats, process_number = \
    dcona.core.utils.read_json(config_path)
    
    dcona.core.utils.check_directory_existence(OUTPUT_DIR_PATH)
    
    data_df = pd.read_csv(DATA_PATH, sep=",", index_col=0)
    description_df = pd.read_csv(DESCRIPTION_PATH, sep=",")
    if INTERACTION_PATH:
        interaction_df = pd.read_csv(INTERACTION_PATH, sep=",")
    else:
        interaction_df = None
    
    sorted_indexes, indexes, \
    df_indexes, ref_corrs, \
    ref_pvalues, exp_corrs, \
    exp_pvalues, stat, pvalue, \
    boot_pvalue, adjusted_pvalue = \
    pipelines.wrapped_ztest(data_df, description_df, interaction_df, \
                            reference_group, experimental_group, correlation, \
                            alternative, score, repeats, process_number)
    
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

    path_to_file = OUTPUT_DIR_PATH.rstrip("/") + "/{}_ztest.csv".format(correlation)
    dcona.core.utils.save_by_chunks(
            sorted_indexes, 
            df_indexes, df_template, df_columns,
            path_to_file,
            # index_transform=None
            index_transform=indexes
    )
    print(f"File saved at: {path_to_file}")
    
    
def zscore_cli(config_path):    
    DATA_PATH, DESCRIPTION_PATH, OUTPUT_DIR_PATH, INTERACTION_PATH, \
    reference_group, experimental_group, correlation, \
    alternative, score, repeats, process_number = \
    dcona.core.utils.read_json(config_path)
    
    dcona.core.utils.check_directory_existence(OUTPUT_DIR_PATH)
    
    data_df = pd.read_csv(DATA_PATH, sep=",", index_col=0)
    description_df = pd.read_csv(DESCRIPTION_PATH, sep=",")
    if INTERACTION_PATH:
        interaction_df = pd.read_csv(INTERACTION_PATH, sep=",")
    else:
        interaction_df = None
    
    result = pipelines.wrapped_zscore(data_df, description_df, interaction_df, \
                                      reference_group, experimental_group, correlation, \
                                      alternative, score, repeats, process_number)
                            
    output_df = pipelines.zscore_to_df(*result)
    
    path_to_file = OUTPUT_DIR_PATH.rstrip("/") + \
                   "/{}_{}_{}_zscore.csv".format(correlation, score, alternative)
    output_df.to_csv(
            path_to_file,
            sep=",",
            index=None
    )
    print(f"File saved at: {path_to_file}")


def main():
    parser = argparse.ArgumentParser(prog='dcona',
                                     description='DCoNA (Differential Correlation Network Analysis)',
                                     epilog='https://github.com/zhiyanov/DCoNA')
    parser.add_argument("tool", choices=['ztest', 'zscore', 'hypergeom'],
                        help="One of DCoNA tools")
    parser.add_argument("config_path",
                        help="Path to JSON config file")
    args = parser.parse_args()
    if args.tool=="ztest":
        ztest_cli(args.config_path)
    elif args.tool=="zscore":
        ztest_cli(args.config_path)

if __name__=="__main__":
    main()
