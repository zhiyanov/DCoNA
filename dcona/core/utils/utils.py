import pandas as pd
import numpy as np
import math
import tqdm
import sys
import os
import itertools

BOOTSTRAP_REPEATS = 10**3


def get_num_ind(indexes, *args):
    index_hash = {
        ind: num for num, ind in enumerate(indexes)
    }
    
    result = []
    for arg in args:
        result.append([
            index_hash[ind] for ind in arg
        ])
    
    if (len(result) == 1):
        return result[0]
    
    return result

def bound(array, left, right):
    array = np.array(array)
    array[array < left] = left
    array[array > right] = right
    return array

def bootstrap_sample(
    *args, statistic=None,
    bootstrap_repeats=BOOTSTRAP_REPEATS
):
    for i in range(bootstrap_repeats):
        indexes = np.random.choice(
            np.arange(len(args[0])),
            len(args[0]),
            replace=True
        )

        samples = []
        for arg in np.array(args):
            samples.append(arg[indexes])

        if (statistic != None):
            yield statistic(*samples)
        else:
            yield sample 
