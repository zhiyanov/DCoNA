#ifndef SCORES_H
#define SCORES_H

const std::string MEAN = "mean";
const std::string MEDIAN = "median";

const float QMEDIAN = 0.5;


float mean(
    float *data_ptr,
    int *index_ptr,
    int start_ind,
    int end_ind,
    bool absolute=false
);

int _mean(
    float *data_ptr,
    int *starts_ind_ptr,
    int *ends_ind_ptr,
    int start_ind,
    int end_ind,
    float *scores_ptr,
    bool absolute=false
);

int __mean(
    float *data_ptr,
    int sources_size,
    int start_ind,
    int end_ind,
    float *scores_ptr,
    bool absolute=false
);

float quantile(
    float *data_ptr,
    int *index_ptr,
    int start_ind,
    int end_ind,
    float q=QMEDIAN,
    bool absolute=false
);

int _quantile(
    float *data_ptr,
    int *starts_ind_ptr,
    int *ends_ind_ptr,
    int start_ind,
    int end_ind,
    float *scores_ptr,
    float q=QMEDIAN,
    bool absolute=false
);

int __quantile(
    float *data_ptr,
    int sources_size,
    int start_ind,
    int end_ind,
    float *scores_ptr,
    float q=QMEDIAN, 
    bool absolute=false
);

#endif 
