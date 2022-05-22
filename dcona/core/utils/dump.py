import tqdm
import sys
import os
import json

from .. import extern

CHUNK_LENGTH = 10**5


def check_directory_existence(
    directory_path,
    message_isdir="Output directory does not exist",
    message_access="Permission denied (output directory)"
):
    if os.path.isdir(directory_path) == False:
        print(message_isdir)
        sys.exit()
    if os.access(directory_path, os.W_OK | os.X_OK) == False:
        print(message_access)
        sys.exit()

def save_by_chunks(
    indexes,
    df_indexes, df_template, df_columns, 
    path_to_file,
    index_transform=None,
    chunk_length=CHUNK_LENGTH
):
    MODE, HEADER = 'w', True
    
    # Splitting rows into chunks
    chunk_number = (len(indexes) + chunk_length - 1) // chunk_length

    for i in tqdm.tqdm(range(chunk_number), desc="Save progress"):
        start = i * chunk_length
        end = (i + 1) * chunk_length
        if i >= chunk_number - 1:
            end = len(indexes)
        
        dump_indexes = indexes[start : end]
        
        source_indexes = []
        target_indexes = []
        for ind in dump_indexes:
            if not (index_transform is None):
                ind = index_transform[ind]

            s, t = extern.paired_index(ind, len(df_indexes))
            source_indexes.append(df_indexes[s])
            target_indexes.append(df_indexes[t])
        
        output_df = df_template.copy()
        output_df["Source"] = source_indexes
        output_df["Target"] = target_indexes
        
        for i in range(len(df_columns)):
            output_df.iloc[:, i + 2] = df_columns[i][dump_indexes]

        output_df.to_csv(
            path_to_file,
            sep=",",
            index=None,
            mode=MODE,
            header=HEADER
        )
        
        MODE, HEADER = 'a', False

def read_json(
    path_to_json
):
    config = json.load(open(path_to_json, "r"))
    DATA_PATH = config["data_path"]
    DESCRIPTION_PATH = config["description_path"]
    OUTPUT_DIR_PATH = config["output_dir_path"]
    REFERENCE_GROUP = config["reference_group"]
    EXPERIMENTAL_GROUP = config["experimental_group"]
    
    CORRELATION = config["correlation"]
    ALTERNATIVE = config["alternative"]
    SCORE = config["score"]
    
    REPEATS_NUMBER = config["repeats_number"]
    PROCESS_NUMBER = config["process_number"]
        
    if ("interaction_path" in config) and (config["interaction_path"] != ""):
        INTERACTION_PATH = config["interaction_path"]
    else:
        INTERACTION_PATH = None
        
    if ("score" in config) and (config["score"] != ""):
        SCORE = config["score"]
    else:
        SCORE = None
    
    return DATA_PATH, DESCRIPTION_PATH, OUTPUT_DIR_PATH, INTERACTION_PATH, \
           REFERENCE_GROUP, EXPERIMENTAL_GROUP, CORRELATION, \
           ALTERNATIVE, SCORE, REPEATS_NUMBER, PROCESS_NUMBER