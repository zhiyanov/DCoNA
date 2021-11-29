import tqdm
import sys
import os

import core.extern

CHUNK_LENGTH = 10**6


def check_directory_existence(
    directory_path,
    message="Output directory does not exist"
):
    if os.path.isdir(directory_path) == False:
        print(message)
        sys.exit()

def save_by_chunks(
    indexes,
    df_indexes, df_template, df_columns, 
    path_to_file,
    chunk_length=CHUNK_LENGTH,
    index_transform=None
):
    df_inds = []
    rows = len(indexes)
    MODE, HEADER = 'w', True
    
    # Splitting rows into chunks
    chunk_number = (len(indexes) + chunk_length - 1) // chunk_length
    for i in tqdm.tqdm(range(chunk_number), desc="Save progress"):
        start = i * chunk_length
        end = (i + 1) * chunk_length
        if i == chunk_number - 1:
            end = i * len(indexes)
        
        dump_indexes = indexes[start : end]
        
        source_indexes = []
        target_indexes = []
        for ind in dump_indexes:
            s, t = core.extern.paired_index(ind, len(df_indexes))
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
