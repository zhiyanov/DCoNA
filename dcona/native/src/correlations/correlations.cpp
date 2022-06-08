#include <cmath>
#include <string>
#include <utility>
#include <queue>
#include <algorithm>

#include "correlations.h"
#include "../utils/utils.h"


int inverse(int *arr, int *reverse, int size) {
    for (int i = 0; i < size; ++i) {
        reverse[arr[i]] = i;
    }
    return 0;
}    

int spearmanr(
    float *data_ptr,
    int sample_size,
    int *source_ind_ptr,
    int *target_ind_ptr,
    float *corrs_ptr,
    int start_ind,
    int end_ind,
    int index_size,
    int *sample_ind_ptr,
    int sample_ind_size
) {
    if (!sample_ind_ptr) {
        sample_ind_size = sample_size;
    }

    int *source_ranks = new int[sample_ind_size];
    int *target_ranks = new int[sample_ind_size];    
    int *reverse = new int[sample_ind_size];
    
    for (int i = start_ind; i < end_ind; ++i) {
        int source_index, target_index;
        if (!source_ind_ptr || !target_ind_ptr) {
            std::pair<int, int> paired_ind =
                paired_index(i, index_size);
            source_index = paired_ind.first; 
            target_index = paired_ind.second;    
        } else {
            source_index = source_ind_ptr[i]; 
            target_index = target_ind_ptr[i];
        }

        float correlation = 0;
        float source_mean = 0, target_mean = 0;
        float source_var  = 0, target_var  = 0;    
        
        range(source_ranks, sample_ind_size);
        range(target_ranks, sample_ind_size);
        
        std::sort(
            source_ranks,
            source_ranks + sample_ind_size,
            [
                data_ptr,
                source_index,
                sample_size,
                sample_ind_ptr
            ](int i1, int i2) {
                if (!sample_ind_ptr) {
                    i1 = sample_size * source_index + i1;
                    i2 = sample_size * source_index + i2;
                } else {
                    i1 = sample_size * source_index + sample_ind_ptr[i1];
                    i2 = sample_size * source_index + sample_ind_ptr[i2];
                }
                return data_ptr[i1] < data_ptr[i2];
            }
        );
        inverse(source_ranks, reverse, sample_ind_size);
        std::swap(source_ranks, reverse);

        std::sort(
            target_ranks,
            target_ranks + sample_ind_size,
            [
                data_ptr,
                target_index,
                sample_size,
                sample_ind_ptr
            ](int i1, int i2) {
                if (!sample_ind_ptr) {
                    i1 = sample_size * target_index + i1;
                    i2 = sample_size * target_index + i2;
                } else {
                    i1 = sample_size * target_index + sample_ind_ptr[i1];
                    i2 = sample_size * target_index + sample_ind_ptr[i2];
                }
                return data_ptr[i1] < data_ptr[i2];
            }
        );
        inverse(target_ranks, reverse, sample_ind_size);
        std::swap(target_ranks, reverse);

        for (int j = 0; j < sample_ind_size; ++j) {
            correlation += (float) (source_ranks[j] * target_ranks[j]);
            
            source_mean += (float) source_ranks[j];
            source_var  += (float) (source_ranks[j] * source_ranks[j]);
            
            target_mean += (float) target_ranks[j];
            target_var  += (float) (target_ranks[j] * target_ranks[j]); 
        }

        source_mean /= sample_ind_size;
        target_mean /= sample_ind_size;

        correlation = correlation / sample_ind_size - 
            source_mean * target_mean;

        source_var = source_var / sample_ind_size -
            source_mean * source_mean; 
        target_var = target_var / sample_ind_size -
            target_mean * target_mean; 
        
        if (source_var == 0 || target_var == 0) { 
            correlation = UNDEFINED_CORR_VALUE;
        } else {
            correlation /= std::sqrt(source_var * target_var);
        }

        corrs_ptr[i] = correlation;
    }
    
    delete[] source_ranks;
    delete[] target_ranks;
    delete[] reverse;

    return 0;
}

int pearsonr(
    float *data_ptr,
    int sample_size,
    int *source_ind_ptr,
    int *target_ind_ptr,
    float *corrs_ptr,
    int start_ind,
    int end_ind,
    int index_size,
    int *sample_ind_ptr,
    int sample_ind_size
) {
    if (!sample_ind_ptr) {
        sample_ind_size = sample_size;
    }

    for (int i = start_ind; i < end_ind; ++i) {
        int source_index, target_index;
        if (!source_ind_ptr || !target_ind_ptr) {
            std::pair<int, int> paired_ind =
                paired_index(i, index_size);
            source_index = paired_ind.first; 
            target_index = paired_ind.second;    
        } else {
            source_index = source_ind_ptr[i]; 
            target_index = target_ind_ptr[i];
        }

        float correlation = 0;
        float source_mean = 0, target_mean = 0;
        float source_var  = 0, target_var  = 0;

        for (int j = 0; j < sample_ind_size; ++j) {
            int jj = (sample_ind_ptr == nullptr) ? j : sample_ind_ptr[j];
            correlation += data_ptr[source_index * sample_size + jj] *
                data_ptr[target_index * sample_size + jj];
            
            source_mean += data_ptr[source_index * sample_size + jj];
            source_var  += data_ptr[source_index * sample_size + jj] *
                           data_ptr[source_index * sample_size + jj];
        
            target_mean += data_ptr[target_index * sample_size + jj];
            target_var  += data_ptr[target_index * sample_size + jj] *
                           data_ptr[target_index * sample_size + jj]; 
        }

        source_mean /= sample_ind_size;
        target_mean /= sample_ind_size;

        correlation = correlation / sample_ind_size - 
            source_mean * target_mean;

        source_var = source_var / sample_ind_size -
            source_mean * source_mean; 
        target_var = target_var / sample_ind_size -
            target_mean * target_mean; 
        
        if (source_var == 0 || target_var == 0) { 
            correlation = UNDEFINED_CORR_VALUE;
        } else {
            correlation /= std::sqrt(source_var * target_var);
        }

        corrs_ptr[i] = correlation;
    }

    return 0;
}
