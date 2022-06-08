def ztest_cli(config_path):
    import pandas as pd 
    from . import lib

    data_path, description_path, \
    reference_group, experimental_group, \
    correlation, alternative, score, \
    interaction_path, \
    repeats_number, \
    output_dir_path, \
    process_number = \
    lib.utils.read_config(config_path)
    
    lib.dump.check_directory_existence(output_dir_path)
    
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
        path_to_file = output_dir_path.rstrip("/") + \
            f"/{correlation}_{alternative}_ztest.csv"
        result.to_csv(path_to_file, sep=",", index=None)
        print(f"File saved at: {path_to_file}") 
    
def zscore_cli(config_path):   
    import pandas as pd
    from . import lib

    data_path, description_path, \
    reference_group, experimental_group, \
    correlation, alternative, score, \
    interaction_path, \
    repeats_number, \
    output_dir_path, \
    process_number = \
    lib.utils.read_config(config_path)
    
    lib.dump.check_directory_existence(output_dir_path)
    
    data_df = pd.read_csv(data_path, sep=",", index_col=0)
    description_df = pd.read_csv(description_path, sep=",")
    if interaction_path:
        interaction_df = pd.read_csv(interaction_path, sep=",")
    else:
        interaction_df = None
    
    result = lib.zscore(
        data_df, description_df,
        reference_group, experimental_group,
        correlation=correlation,
        score=score,
        alternative=alternative,
        interaction=interaction_df,
        repeats_number=repeats_number,
        output_dir=output_dir_path,
        process_number=process_number
    )
                            
    if not (result is None): 
        path_to_file = output_dir_path.rstrip("/") + \
            f"/{correlation}_{score}_{alternative}_zscore.csv"
        result.to_csv(path_to_file, sep=",", index=None)
        print(f"File saved at: {path_to_file}") 

def hypergeom_cli(config_path):   
    import pandas as pd
    from . import lib

    data_path, description_path, \
    reference_group, experimental_group, \
    correlation, alternative, score, \
    interaction_path, \
    repeats_number, \
    output_dir_path, \
    process_number = \
    lib.utils.read_config(config_path)
    
    lib.dump.check_directory_existence(output_dir_path)

    path_to_file = output_dir_path.rstrip("/") + \
        f"/{correlation}_{alternative}_ztest.csv"
    ztest_df = pd.read_csv(path_to_file, sep=",")
    
    if interaction_path:
        oriented = True
    else:
        oriented = False

    result = lib.hypergeom(
        ztest_df,
        alternative=alternative,
        oriented=oriented,
        output_dir=output_dir_path
    )
                            
    if not (result is None): 
        path_to_file = output_dir.rstrip("/") + \
            f"/{alternative}_hypergeom.csv"
        result.to_csv(
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
