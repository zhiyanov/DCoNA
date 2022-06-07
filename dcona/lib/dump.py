import tqdm
import sys
import os

from ..core import extern as cextern

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

    for i in tqdm.tqdm(range(chunk_number), desc="Save progress: ", ascii=True):
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

            s, t = cextern.paired_index(ind, len(df_indexes))
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

