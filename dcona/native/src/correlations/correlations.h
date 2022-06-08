#ifndef CORRS_H
#define CORRS_H

#include <utility>
#include <vector>
#include <string>

const float UNDEFINED_CORR_VALUE = -2;

const std::string SPEARMAN = "spearman";
const std::string PEARSON = "pearson";

int spearmanr(
    float *data_ptr,
    int sample_size,
    int *source_ind_ptr,
    int *target_ind_ptr,
    float *corrs_ptr,
    int start_ind,
    int end_ind,
    int index_size=-1,
    int *sample_ind_ptr=nullptr,
    int sample_ind_size=-1
);

int pearsonr(
    float *data_ptr,
    int sample_size,
    int *source_ind_ptr,
    int *target_ind_ptr,
    float *corrs_ptr,
    int start_ind,
    int end_ind,
    int index_size=-1,
    int *sample_ind_ptr=nullptr,
    int sample_ind_size=-1
);

#endif

