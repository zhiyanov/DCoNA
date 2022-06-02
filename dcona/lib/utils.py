import itertools
import pandas as pd
import numpy as np


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
