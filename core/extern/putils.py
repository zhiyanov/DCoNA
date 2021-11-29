import numpy as np

from .utils import UNDEFINED_INDEX
from .utils import \
    _unary_index, \
    _paired_index, \
    _unary_vector, \
    _unary_matrix, \
    _quadrate, \
    _reorder, \
    _reorder_data


def paired_index(index, base):
    return _paired_index(index, base)

def unary_index(first, second, base):
    result = _unary_index(first, second, base)
    
    if result == UNDEFINED_INDEX:
        return None
    
    return result

def unary_array(index, base):
    index = np.array(index, dtype="int32")
    return _unary_vector(index, base)

def unary_matrix(index, base):
    index = index.astyp("int32")
    return _unary_matrix(indx, base)

def quadrate(flatten_array, index, base):
    index = np.array(index, dtype="int32")
    flatten_array = np.array(array, dtype="float32") 
    return _quadrate(flatten_array, index, base)

def reorder(source_indexes, target_indexes, data=None):
    if data:
        _reorder_data(source_indexes, target_indexes, data)
    else:
        _reorder(source_indexes, target_indexes)
