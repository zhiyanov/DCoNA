def ztest_cli(config_path):
    import pandas as pd
    
    from .core import utils as cutils
    import .lib

    data_path, description_path, output_dir_path, interaction_path, \
    reference_group, experimental_group, correlation, \
    alternative, score, repeats_number, process_number = \
    cutils.read_json(config_path)
    
    cutils.check_directory_existence(output_dir_path)
    
    data_df = pd.read_csv(data_path, sep=",", index_col=0)
    description_df = pd.read_csv(description_path, sep=",")
    if interaction_path:
        interaction_df = pd.read_csv(interaction_path, sep=",")
    else:
        interaction_df = None
    
    result = lib.ztest(
        data_df, description_df,
        reference_group, experimental_group,
        correlation, alternative,
        interaction=interaction_df,
        repeats_number=repeats_number,
        output_dir=output_dir_path,
        process_number=process_number
    )

    if not (result is None):        
        path_to_file = output_dir_path.rstrip("/") + "/{}_ztest.csv".format(correlation)
        result.to_csv(path_to_file, sep=",", index=None)
        print(f"File saved at: {path_to_file}") 
    
def zscore_cli(config_path):   
    import pandas as pd
    
    from .core import utils as cutils
    import .lib

    data_path, description_path, output_dir_path, interaction_path, \
    reference_group, experimental_group, correlation, \
    alternative, score, repeats, process_number = \
    dcona.core.utils.read_json(config_path)
    
    dcona.core.utils.check_directory_existence(output_dir_path)
    
    data_df = pd.read_csv(DATA_PATH, sep=",", index_col=0)
    description_df = pd.read_csv(description_path, sep=",")
    if interaction_path:
        interaction_df = pd.read_csv(interaction_path, sep=",")
    else:
        interaction_df = None
    
    result = lib.zscore(data_df, description_df, interaction_df, \
                                      reference_group, experimental_group, correlation, \
                                      alternative, score, repeats, process_number)
                            
    output_df = pipelines.zscore_to_df(*result)
    
    path_to_file = output_dir_path.rstrip("/") + \
                   "/{}_{}_{}_zscore.csv".format(correlation, score, alternative)
    output_df.to_csv(
            path_to_file,
            sep=",",
            index=None
    )
    print(f"File saved at: {path_to_file}")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        prog='dcona',
        description='DCoNA (Differential Correlation Network Analysis)',
        epilog='https://github.com/zhiyanov/DCoNA'
    )
    parser.add_argument(
        "tool", choices=['ztest', 'zscore', 'hypergeom'],
        help="One of DCoNA tools"
    )
    parser.add_argument(
        "config_path",
        help="Path to JSON config file"
    )

    args = parser.parse_args()
    if args.tool=="ztest":
        ztest_cli(args.config_path)
    elif args.tool=="zscore":
        zscore_cli(args.config_path)
    else:
        hypergeom_cli(args.config_path)

if __name__=="__main__":
    main()
