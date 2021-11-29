import numpy as np

from .scores import \
    _score_indexed, \
    _score_exhaustive

from .putils import \
    reorder 


def score_indexed(
    data,
    source_indexes,
    target_indexes,
    score="mean",
    process_num=1 
):
    data = data.astype("float32")
    source_indexes = source_indexes.astype("int32")
    target_indexes = target_indexes.astype("int32")
    
    reorder(source_indexes, target_indexes, data)
    
    sources, scores, pvalues = _score_indexed(
        data,
        source_indexes,
        score,
        process_num
    )

    return sources, scores, pvalues

def score_exhaustive(
    data,
    sources_size,
    score="mean",
    process_num=1 
):
    data = data.astype("float32")
    scores, pvalues = _score_exhaustive(
        data,
        sources_size,
        score,
        process_num
    )

    sources = np.arange(sources_size, dtype="int32")

    return sources, scores, pvalues
