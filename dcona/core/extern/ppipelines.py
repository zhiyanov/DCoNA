import numpy as np

from .pipelines import \
    _ztest_pipeline_indexed, \
    _score_pipeline_indexed, \
    _ztest_pipeline_exhaustive, \
    _score_pipeline_exhaustive

from .pcorrelations import \
    correlation_test \

from .putils import \
    reorder


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

def ztest_pipeline( 
    df,
    reference_indexes,
    experimental_indexes,
    source_indexes=None,
    target_indexes=None,
    correlation="spearman",
    correlation_alternative=False,
    alternative="two-sided",
    repeats_num=1000,
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
            ref_num_indexes, exp_num_indexes = \
                get_num_ind(
                    df.columns.to_list(),
                    reference_indexes,
                    experimental_indexes
                )
        else:
            source_num_indexes = source_indexes
            target_num_indexes = target_indexes
            
            ref_num_indexes = reference_indexes
            exp_num_indexes = experimental_indexes

        source_num_indexes = np.array(
            source_num_indexes
        ).astype("int32")
        target_num_indexes = np.array(
            target_num_indexes
        ).astype("int32")

        ref_num_indexes = np.array(
            ref_num_indexes
        ).astype("int32")
        exp_num_indexes = np.array(
            exp_num_indexes
        ).astype("int32")

        ref_corrs, exp_corrs, \
        stat, pvalue, bootstrap_pvalue = \
            _ztest_pipeline_indexed(
                data,
                source_num_indexes,
                target_num_indexes,
                ref_num_indexes,
                exp_num_indexes,
                correlation,
                alternative,
                repeats_num,
                process_num
            ) 
    else:
        if not numerical_index:
            ref_num_indexes, exp_num_indexes = \
                get_num_ind(
                    df.columns.to_list(),
                    reference_indexes,
                    experimental_indexes
                )
        else:
            ref_num_indexes = reference_indexes
            exp_num_indexes = experimental_indexes
        
        ref_num_indexes = np.array(
            ref_num_indexes
        ).astype("int32")
        exp_num_indexes = np.array(
            exp_num_indexes
        ).astype("int32")

        ref_corrs, exp_corrs, \
        stat, pvalue, bootstrap_pvalue = \
            _ztest_pipeline_exhaustive(
                data,
                ref_num_indexes,
                exp_num_indexes,
                correlation,
                alternative,
                repeats_num,
                process_num
            ) 
    
    if correlation_alternative:
        ref_pvalues = correlation_test(
            ref_corrs,
            len(ref_num_indexes),
            correlation=correlation,
            alternative=correlation_alternative
        )

        exp_pvalues = correlation_test(
            exp_corrs,
            len(exp_num_indexes),
            correlation=correlation,
            alternative=correlation_alternative
        )

        return ref_corrs, ref_pvalues, exp_corrs, \
                exp_pvalues, stat, pvalue, bootstrap_pvalue

    return ref_corrs, exp_corrs, \
            stat, pvalue, bootstrap_pvalue

def score_pipeline( 
    df,
    reference_indexes,
    experimental_indexes,
    source_indexes=None,
    target_indexes=None,
    correlation="spearman",
    score="mean",
    alternative="two_sided",
    repeats_num=1000,
    process_num=1,
    numerical_index=False,
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
            ref_num_indexes, exp_num_indexes = \
                get_num_ind(
                    df.columns.to_list(),
                    reference_indexes,
                    experimental_indexes
                )
        else:
            source_num_indexes = source_indexes
            target_num_indexes = target_indexes
            
            ref_num_indexes = reference_indexes
            exp_num_indexes = experimental_indexes

        source_num_indexes = np.array(
            source_num_indexes
        ).astype("int32")
        target_num_indexes = np.array(
            target_num_indexes
        ).astype("int32")

        ref_num_indexes = np.array(
            ref_num_indexes
        ).astype("int32")
        exp_num_indexes = np.array(
            exp_num_indexes
        ).astype("int32")
        
        reorder(source_num_indexes, target_num_indexes)
        
        indexes, scores, pvalues = \
            _score_pipeline_indexed(
                data,
                source_num_indexes,
                target_num_indexes,
                ref_num_indexes,
                exp_num_indexes,
                correlation,
                score,
                alternative,
                repeats_num,
                process_num
            )
    else:
        if not numerical_index:
            ref_num_indexes, exp_num_indexes = \
                get_num_ind(
                    df.columns.to_list(),
                    reference_indexes,
                    experimental_indexes
                )
        else:
            ref_num_indexes = reference_indexes
            exp_num_indexes = experimental_indexes
        
        ref_num_indexes = np.array(
            ref_num_indexes
        ).astype("int32")
        exp_num_indexes = np.array(
            exp_num_indexes
        ).astype("int32")
        
        scores, pvalues = \
            _score_pipeline_exhaustive(
                data,
                ref_num_indexes,
                exp_num_indexes,
                correlation,
                score,
                alternative,
                repeats_num,
                process_num
            )

        indexes = np.arange(data.shape[0], dtype="int32")
    
    return indexes, scores, pvalues
