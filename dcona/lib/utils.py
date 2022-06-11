import itertools
import pandas as pd
import numpy as np
import json


def generate_pairs(objects):    
    if isinstance(objects, (list, set, pd.core.series.Series, np.ndarray)):
        pairs = list(itertools.combinations(objects, 2))
    elif isinstance(objects, str):
        pairs = objects.replace(",", " ").split(" ")
        pairs = list(filter(lambda a: a != "", pairs))
        pairs = list(itertools.combinations(pairs, 2))
    else:
        try:
            pairs = list(itertools.combinations(objects, 2))
        except TypeError:
            print("Use appropriate data type: list, np.array, string, etc.")
            raise
        except:
            raise
    
    pairs_df = pd.DataFrame(pairs, columns=["Source", "Target"])
    return pairs_df

# Reads config from a json file
def read_config(
    path_to_config
):
    config = json.load(open(path_to_config, "r"))
    data_path = config["data_path"]
    description_path = config["description_path"]
    reference_group = config["reference_group"]
    experimental_group = config["experimental_group"]
    
    correlation = config["correlation"]
    alternative = config["alternative"]
    score = config["score"]
    
    output_dir_path = config["output_dir_path"]
    repeats_number = config["repeats_number"]
    process_number = config["process_number"]
        
    if ("interaction_path" in config) and (config["interaction_path"] != ""):
        interaction_path = config["interaction_path"]
    else:
        interaction_path = None
        
    if ("score" in config) and (config["score"] != ""):
        score = config["score"]
    else:
        score = None
   
    return data_path, description_path, \
        reference_group, experimental_group, \
        correlation, alternative, score, \
        interaction_path, \
        repeats_number, \
        output_dir_path, \
        process_number
