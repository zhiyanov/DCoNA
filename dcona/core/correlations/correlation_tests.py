"""Tests of equal correlation hypothesis

The module contains methods testing hypothesis of
correlation equality
"""


import numpy as np
import scipy.stats
from .correlation_utils import \
    pearsonr_mean, pearsonr_std

def ztest(
    first_rs, first_size,
    second_rs, second_size,
    correlation="spearman",
    alternative="two-sided"
):
    """Check the hypothesis that pearson correlations
    are equal

    Parameters
    ----------
    first_rs: numerical value or list
        A sequence (potentially with one element)
        of correlations.
    first_size: numerical value or list
        Sizes of samples that were used to compute
        "first_rs" correlation(s).
    second_rs: numerical value or list
        A sequence (potentially with one element)
        of correlations.
    second_size: numerical value or list
        Sizes of samples that were used to compute
        "second_rs" correlation(s).
    alternative: "two-sided" (default), "less", "greater"
        Computes the probability of the following events:
        "two-sided" |arctanh(x1) - arctanh(x2)| >
        |arctanh(first_rs) - arctanh(second_rs)|,
        "less" arctanh(x1) - arctanh(x2) <=
        arctanh(first_rs) - arctanh(second_rs) ,
        "greater" arctanh(x1) - arctanh(x2) >
        arctanh(first_rs) - arctanh(second_rs).
    
    Returns
    -------
    pair of numerical values or numpy.arrays respectively to the input
        Contains statistic and pvalue.
    """

    if len(first_rs) != len(second_rs):
        return None
    result_len = len(first_rs)

    first_rs = np.array(first_rs)
    first_size = np.array(first_size)
    
    second_rs = np.array(second_rs)
    second_size = np.array(second_size)
    
    bound_indexes = (np.abs(first_rs + second_rs) == 2) | \
        (first_rs == None) | (second_rs == None)
    bound_indexes = ~bound_indexes

    first_rs = first_rs[bound_indexes]
    second_rs = second_rs[bound_indexes]

    first_ss = pearsonr_std(first_rs, first_size)
    second_ss = pearsonr_std(second_rs, second_size)
    
    if (correlation=="spearman"):
        first_ss *= np.sqrt(1.5)
        second_ss *= np.sqrt(1.5)

    stat = np.arctanh(first_rs) - np.arctanh(second_rs)
    std = np.sqrt(first_ss**2 + second_ss**2)
    
    pvalue = None
    
    if (alternative == "less"):
        pvalue = scipy.stats.norm.cdf(stat, scale=std)
    elif (alternative == "greater"):
        pvalue = 1 - scipy.stats.norm.cdf(stat, scale=std)
    elif (alternative == "two-sided"):
        pvalue = 2 * scipy.stats.norm.cdf(-np.abs(stat), scale=std)

    stat_result = np.zeros(result_len, dtype="float32")
    pvalue_result = np.zeros(result_len, dtype="float32")
    
    stat_result[~bound_indexes] = None
    pvalue_result[~bound_indexes] = None
    
    stat_result[bound_indexes] = stat
    pvalue_result[bound_indexes] = pvalue

    return stat_result, pvalue_result
