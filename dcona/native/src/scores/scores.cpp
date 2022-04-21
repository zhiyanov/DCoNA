#include <cmath>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include <iostream>

#include "scores.h"
#include "../utils/utils.h"


float mean(
    float *data_ptr,
    int *index_ptr,
    int start_ind,
    int end_ind,
    bool absolute
) {
    float mean = 0;
    float iter_num = 0;

    int index = 0;
    for (int i = start_ind; i < end_ind; ++i) {
        if (!index_ptr) {
            index = i;
        } else {
            index = index_ptr[i];
        }

        if (index < 0) {
            continue;
        }

        mean += std::abs(data_ptr[index]);
        iter_num += 1;
    }

    if (iter_num > 0) {
        mean /= iter_num;
    }

    return mean;
}

int _mean(
    float *data_ptr,
    int *starts_ind_ptr,
    int *ends_ind_ptr,
    int start_ind,
    int end_ind,
    float *scores_ptr,
    bool absolute
) {
    for (int i = start_ind; i < end_ind; ++i) {
        scores_ptr[i] = mean(
            data_ptr,
            (int *) nullptr,
            starts_ind_ptr[i],
            ends_ind_ptr[i],
            absolute
        );  
    }

    return 0;
}

int __mean(
    float *data_ptr,
    int sources_size,
    int start_ind,
    int end_ind,
    float *scores_ptr,
    bool absolute
) {
    for (int i = start_ind; i < end_ind; ++i) {
        std::vector<int> targets = unary_vector(i, sources_size);
        scores_ptr[i] = mean(
            data_ptr,
            targets.data(),
            0,
            sources_size,
            absolute
        );  
    }

    return 0;
}

float quantile(
    float *data_ptr,
    int *index_ptr,
    int start_ind,
    int end_ind,
    float q,
    bool absolute
) {
    std::vector<float> values;
    int iter_num = 0;

    int index = 0;

    for (int i = start_ind; i < end_ind; ++i) {
        if (!index_ptr) {
            index = i;
        } else {
            index = index_ptr[i];
        }

        if (index < 0) {
            continue;
        }
        
        iter_num += 1;
        if (absolute) {
            values.push_back(std::abs(data_ptr[index]));
        } else {
            values.push_back(data_ptr[index]);
        }
    }
 
    if (iter_num > 0) {
        std::sort(values.begin(), values.end());
        return values[(int) iter_num * q];
    }

    return 0;
}

int _quantile(
    float *data_ptr,
    int *starts_ind_ptr,
    int *ends_ind_ptr,
    int start_ind,
    int end_ind,
    float *scores_ptr,
    float q,
    bool absolute
) {
    for (int i = start_ind; i < end_ind; ++i) {
        scores_ptr[i] = quantile(
            data_ptr,
            (int *) nullptr,
            starts_ind_ptr[i],
            ends_ind_ptr[i],
            q,
            absolute
        );  
    }

    return 0;
}

int __quantile(
    float *data_ptr,
    int sources_size,
    int start_ind,
    int end_ind,
    float *scores_ptr,
    float q,
    bool absolute
) {
    for (int i = start_ind; i < end_ind; ++i) {
        std::vector<int> targets = unary_vector(i, sources_size);
        scores_ptr[i] = quantile(
            data_ptr,
            targets.data(),
            0,
            sources_size,
            q,
            absolute
        );  
    }

    return 0;
}
