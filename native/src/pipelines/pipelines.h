#ifndef PIPELINES_H
#define PIPELINES_H

#include <string>

#include "../correlations/correlations.h"
#include "../tests/tests.h"
#include "../scores/scores.h"

const int REPEATS_NUMBER = 1000;


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
    const std::string correlation=SPEARMAN,
    const std::string alternative=TWO_SIDED
);

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
    const std::string correlation=SPEARMAN,
    const std::string alternative=TWO_SIDED,
    const std::string score=MEAN
);

int score_pipeline_exhaustive(
    float *stat_ptr,
    int index_size,
    int start_ind,
    int end_ind,
    float *score_ptr,
    const std::string score=MEAN,
    const std::string alternative=TWO_SIDED
);

#endif 
