import numpy as np
import pandas as pd
import scipy.stats

from ..core import extern
from ..core import utils


def hypergeom(
    one, two, three
):
    
    return None

def get_occurrence(array):
    occurrence = {}
    for elem in array:
        if elem in occurrence:
            occurrence[elem] += 1
        else:
            occurrence[elem] = 1
    
    return occurrence
