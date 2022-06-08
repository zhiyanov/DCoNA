import numpy as np
from scipy.stats import rankdata
from scipy.stats import t
from scipy.stats import norm

from .correlations import UNDEFINED_CORR_VALUE
from .correlations import \
    _correlation_indexed, \
    _correlation_exhaustive

TWO_SIDED = "two-sided"
LESS = "less"
GREATER = "greater"

SPEARMAN = "spearman"
PEARSON = "pearson"

LEFT_CORR_BOUND = -0.99
RIGHT_CORR_BOUND = 0.99

def correlation_test(
    corrs,
    size,
    correlation=SPEARMAN,
    alternative=TWO_SIDED
):
    corrs[corrs < LEFT_CORR_BOUND] = LEFT_CORR_BOUND
    corrs[corrs > RIGHT_CORR_BOUND] = RIGHT_CORR_BOUND

    tcorrs = np.arctanh(corrs)
    if correlation == SPEARMAN:
        tcorrs = tcorrs * np.sqrt((size - 3) / (1 + corrs**2 / 2))
    elif correlation == PEARSON:
        tcorrs = tcorrs * np.sqrt(size - 3)

    if alternative == TWO_SIDED:
        pvalues = 2 * norm.cdf(-np.abs(tcorrs))
    elif alternative == LESS: 
        pvalues = norm.cdf(tcorrs)
    elif alternative == GREATER:
        pvalues = 1 - norm.cdf(tcorrs)

    return pvalues

def get_num_ind(indexes, *args):
    index_hash = {
        ind: num for num, ind in enumerate(indexes)
    }
    
    result = []
    for arg in args:
        result.append([
            index_hash[ind] for ind in arg
        ])

    return result 

def spearmanr(
    df,
    source_indexes=None,
    target_indexes=None,
    alternative=None,
    process_num=1,
    numerical_index=False
): 
    data = df.to_numpy(copy=True).astype("float32")

    if np.all(source_indexes != None) and np.all(target_indexes != None):
        if not numerical_index:
            source_num_indexes, target_num_indexes = \
                get_num_ind(
                    df.index.to_list(),
                    source_indexes,
                    target_indexes
                )
        else:
            source_num_indexes = source_indexes
            target_num_indexes = target_indexes

        source_num_indexes = np.array(
            source_num_indexes
        ).astype("int32")
        target_num_indexes = np.array(
            target_num_indexes
        ).astype("int32")

        corrs = _correlation_indexed(
            data,
            source_num_indexes,
            target_num_indexes,
            "spearman",
            process_num
        ) 
    else:
        corrs = _correlation_exhaustive(
            data,
            "spearman",
            process_num
        )

    corrs[corrs == UNDEFINED_CORR_VALUE] = None
    
    if alternative:
        pvalues = correlation_test(
            corrs,
            data.shape[1],
            correlation=SPEARMAN,
            alternative=alternative
        )

        return corrs, pvalues
    
    return corrs

def pearsonr(
    df,
    source_indexes=None,
    target_indexes=None,
    alternative=None,
    process_num=1,
    numerical_index=False
):    
    data = df.to_numpy(copy=True).astype("float32")
    
    if np.all(source_indexes != None) and np.all(target_indexes != None):
        if not numerical_index:
            source_num_indexes, target_num_indexes = \
                get_num_ind(
                    df.index.to_list(),
                    source_indexes,
                    target_indexes
                )
        else:
            source_num_indexes = source_indexes
            target_num_indexes = target_indexes
     
        source_num_indexes = np.array(
            source_num_indexes
        ).astype("int32")
        target_num_indexes = np.array(
            target_num_indexes
        ).astype("int32")
    
        corrs = _correlation_indexed(
            data,
            source_num_indexes,
            target_num_indexes,
            "pearson",
            process_num
        )
    else:
        corrs = _correlation_exhaustive(
            data,
            "pearson",
            process_num
        )

    corrs[corrs == UNDEFINED_CORR_VALUE] = None

    if alternative:
        pvalues = correlation_test(
            corrs,
            data.shape[1],
            correlation=PEARSON,
            alternative=alternative
        )

        return corrs, pvalues
    
    return corrs

def spearmanr_test(
    df,
    source_indexes=None,
    target_indexes=None,
    alternative=None,
    process_num=1,
    numerical_index=False
): 
    data = df.to_numpy(copy=True).astype("float32")

    if np.all(source_indexes != None) and np.all(target_indexes != None):
        if not numerical_index:
            source_num_indexes, target_num_indexes = \
                get_num_ind(
                    df.index.to_list(),
                    source_indexes,
                    target_indexes
                ) 
        else:
            source_num_indexes = source_indexes
            target_num_indexes = target_indexes

        source_num_indexes = np.array(
            source_num_indexes
        ).astype("int32")
        target_num_indexes = np.array(
            target_num_indexes
        ).astype("int32")

        indexes = set(source_num_indexes)
        indexes.update(target_num_indexes)
        indexes = list(indexes)
        
        data[indexes] = rankdata(
            data[indexes],
            axis=1,
            method="ordinal"
        ).astype("float32")
        
        corrs = _correlation_indexed(
            data,
            source_num_indexes,
            target_num_indexes,
            "pearson",
            process_num
        ) 
    else:
        data = rankdata(
            data,
            axis=1,
            method="ordinal"
        ).astype("float32")

        corrs = _correlation_exhaustive(
            data,
            "pearson",
            process_num
        )

    corrs[corrs == UNDEFINED_CORR_VALUE] = None

    if alternative:
        pvalues = correlation_test(
            corrs,
            data.shape[1],
            correlation=SPEARMAN,
            alternative=alternative
        )

        return corrs, pvalues
    
    return corrs
