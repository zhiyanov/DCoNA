from .putils import \
    unary_index, \
    paired_index, \
    unary_array, \
    unary_matrix, \
    quadrate, \
    reorder

from .pcorrelations import \
    spearmanr, \
    spearmanr_test, \
    pearsonr
 
from .ptests import \
    ztest

from .pscores import \
    score_indexed, \
    score_exhaustive 

from .ppipelines import \
    ztest_pipeline, \
    score_pipeline
