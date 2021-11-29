import numpy as np

from .tests import \
    _ztest_sized, \
    _ztest_unsized


def ztest(
    first_rs, first_size,
    second_rs, second_size,
    correlation="spearman",
    alternative="two-sided",
    process_num=1 
):
    first_rs = first_rs.astype("float32")
    if hasattr(first_size, "__iter__"):
        first_size = np.array(first_size, dtype="int32") 
    
    second_rs = second_rs.astype("float32")
    if hasattr(second_size, "__iter__"):
        second_size = np.array(second_size, dtype="int32") 
    
    if hasattr(first_size, "__iter__") and \
            hasattr(second_size, "__iter__"):
        stat, pvalue = _ztest_sized(
            first_rs,
            first_size,
            second_rs,
            second_size,
            correlation,
            alternative,
            process_num
        )
    else:
        stat, pvalue = _ztest_unsized(
            first_rs,
            first_size,
            second_rs,
            second_size,
            correlation,
            alternative,
            process_num
        )
    
    return stat, pvalue
