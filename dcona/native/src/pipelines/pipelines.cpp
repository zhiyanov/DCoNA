#include <string>
#include <utility>
#include <iostream>

#include "../correlations/correlations.h"
#include "../tests/tests.h"
#include "../scores/scores.h"
#include "pipelines.h"

int ztest_pipeline(
    float *data_ptr,
    int sample_size,
    int *source_ind_ptr,
    int *target_ind_ptr,
    int start_ind,
    int end_ind,
    int index_size,
    int *ref_ind_ptr,
    int *exp_ind_ptr,
    int ref_ind_size,
    int exp_ind_size,
    float *ref_corrs_ptr,
    float *exp_corrs_ptr,
    float *stat_ptr,
    float *pvalue_ptr,
    const std::string correlation,
    const std::string alternative
) {
    // Can be used in exhaustive and
    // interaction modes

    // std::cout << "Correlation computations\n";
    if (correlation == SPEARMAN) {
        pearsonr(
            data_ptr,
            sample_size,
            source_ind_ptr,
            target_ind_ptr,
            ref_corrs_ptr,
            start_ind,
            end_ind,
            index_size,
            ref_ind_ptr,
            ref_ind_size
        );
        
        pearsonr(
            data_ptr,
            sample_size,
            source_ind_ptr,
            target_ind_ptr,
            exp_corrs_ptr,
            start_ind,
            end_ind,
            index_size,
            exp_ind_ptr,
            exp_ind_size
        );
    } else {
        pearsonr(
            data_ptr,
            sample_size,
            source_ind_ptr,
            target_ind_ptr,
            ref_corrs_ptr,
            start_ind,
            end_ind,
            index_size,
            ref_ind_ptr,
            ref_ind_size
        );

        pearsonr(
            data_ptr,
            sample_size,
            source_ind_ptr,
            target_ind_ptr,
            exp_corrs_ptr,
            start_ind,
            end_ind,
            index_size,
            exp_ind_ptr,
            exp_ind_size
        );
    }

    // std::cout << "Z-test computations\n";
    ztest_unsized(
        ref_corrs_ptr, ref_ind_size,
        exp_corrs_ptr, exp_ind_size,
        stat_ptr, pvalue_ptr,
        start_ind, end_ind,
        correlation,
        alternative    
    );
    
    return 0;
}

int score_pipeline_indexed(
    float *data_ptr,
    int sample_size,
    int *source_ind_ptr,
    int *target_ind_ptr,
    int index_size,
    int *ref_ind_ptr,
    int *exp_ind_ptr,
    int ref_ind_size,
    int exp_ind_size,
    float *ref_corrs_ptr,
    float *exp_corrs_ptr,
    float *stat_ptr,
    float *pvalue_ptr,
    int *starts_ind_ptr,
    int *ends_ind_ptr,
    int start_ind,
    int end_ind,
    float *score_ptr,
    const std::string correlation,
    const std::string alternative,
    const std::string score
) { 
    int start = starts_ind_ptr[start_ind];
    int end = ends_ind_ptr[end_ind - 1];
    
    // std::cout << "Correlation computations\n";
    if (correlation == SPEARMAN) {
        pearsonr(
            data_ptr,
            sample_size,
            source_ind_ptr,
            target_ind_ptr,
            ref_corrs_ptr,
            start,
            end,
            index_size,
            ref_ind_ptr,
            ref_ind_size
        );
        
        pearsonr(
            data_ptr,
            sample_size,
            source_ind_ptr,
            target_ind_ptr,
            exp_corrs_ptr,
            start,
            end,
            index_size,
            exp_ind_ptr,
            exp_ind_size
        );
    } else {
        pearsonr(
            data_ptr,
            sample_size,
            source_ind_ptr,
            target_ind_ptr,
            ref_corrs_ptr,
            start,
            end,
            index_size,
            ref_ind_ptr,
            ref_ind_size
        );

        pearsonr(
            data_ptr,
            sample_size,
            source_ind_ptr,
            target_ind_ptr,
            exp_corrs_ptr,
            start,
            end,
            index_size,
            exp_ind_ptr,
            exp_ind_size
        );
    }

    // std::cout << "Z-test computations\n";
    ztest_unsized(
        ref_corrs_ptr,
        ref_ind_size,
        exp_corrs_ptr,
        exp_ind_size,
        stat_ptr,
        pvalue_ptr,
        start,
        end,
        correlation,
        TWO_SIDED    
    );

    bool absolute = false;
    if (alternative == TWO_SIDED) {
        absolute = true;
    }

    // std::cout << "Score computations\n";
    if (score == MEAN) {
        _mean(
            stat_ptr,
            starts_ind_ptr,
            ends_ind_ptr,
            start_ind,
            end_ind,
            score_ptr,
            absolute
        ); 
    } else if (score == MEDIAN) {
        _quantile( 
            stat_ptr,
            starts_ind_ptr,
            ends_ind_ptr,
            start_ind,
            end_ind,
            score_ptr,
            QMEDIAN,
            absolute
        );
    }
    
    return 0;
}

int score_pipeline_exhaustive(
    float *stat_ptr,
    int index_size,
    int start_ind,
    int end_ind,
    float *score_ptr,
    const std::string score,
    const std::string alternative
) {
    bool absolute = false;
    if (alternative == TWO_SIDED) {
        absolute = true;
    }

    if (score == MEAN) {
        __mean(
            stat_ptr,
            index_size,
            start_ind,
            end_ind,
            score_ptr,
            absolute
        ); 
    } else if (score == MEDIAN) {
        __quantile( 
            stat_ptr,
            index_size,
            start_ind,
            end_ind,
            score_ptr,
            QMEDIAN,
            absolute
        );
    }
    
    return 0;
}
