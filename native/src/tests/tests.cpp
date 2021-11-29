#include <cmath>
#include <string>
#include <utility>
#include <queue>
#include <tuple>
#include <algorithm>
#include <random>
#include <vector>

#include "tests.h"


float norm_cdf(float x) {
    return std::erfc(-x / std::sqrt(2)) / 2;
}

int ztest_unsized(
    float *first_rs_ptr, int first_size,
    float *second_rs_ptr, int second_size,
    float *stat_ptr, float *pvalue_ptr,
    int start_ind, int end_ind,
    const std::string correlation,
    const std::string alternative
) {
    for (int ind = start_ind; ind < end_ind; ++ind) {
        float first_rs = first_rs_ptr[ind];
        float second_rs = second_rs_ptr[ind];
        
        // Bound corrs
        if (first_rs < LEFT_CORR_BOUND) {
            first_rs = LEFT_CORR_BOUND;
        }

        if (second_rs < LEFT_CORR_BOUND) {
            second_rs = LEFT_CORR_BOUND;
        }
        
        if (first_rs > RIGHT_CORR_BOUND) {
            first_rs = RIGHT_CORR_BOUND;
        }

        if (second_rs > RIGHT_CORR_BOUND) {
            second_rs = RIGHT_CORR_BOUND;
        }

        float stat = (std::atanh(first_rs) -
                std::atanh(second_rs));    
        
        // Note: we can use student distribution
        // in the case of spearman correlation
        
        // This block is a copy of 
        // core.correlation_utils.pearson_std
        float first_ss = 1 / std::sqrt(first_size - 3);
        float second_ss = 1 / std::sqrt(second_size - 3);
        
        // It is a maximal bound of the statistic
        // variation: sqrt(1 + r^2 / 2)
        if (correlation == SPEARMAN) {
            first_ss *= std::sqrt(1.5);
            second_ss *= std::sqrt(1.5);
        }

        float std = std::sqrt(first_ss * first_ss +
                second_ss * second_ss);    
        
        if (pvalue_ptr != nullptr) {
            float pvalue = UNDEFINED_CORR_DIFF_TEST_VALUE;
            if (alternative == LESS) {
                pvalue = norm_cdf(stat / std);    
            } else if (alternative == GREATER) {
                pvalue = 1 - norm_cdf(stat / std);
            } else if (alternative == TWO_SIDED) {
                pvalue = 2 * norm_cdf(-std::abs(stat) / std);            
            }
            
            pvalue_ptr[ind] = pvalue;
        }
   
        stat_ptr[ind] = stat;
    }

    return 0;
}

int ztest_sized(
    float *first_rs_ptr, int *first_size_ptr,
    float *second_rs_ptr, int *second_size_ptr,
    float *stat_ptr, float *pvalue_ptr,
    int start_ind, int end_ind,
    const std::string correlation,
    const std::string alternative
) {
    for (int ind = start_ind; ind < end_ind; ++ind) {
          float first_rs = first_rs_ptr[ind]; 
          float second_rs = second_rs_ptr[ind]; 
        
        // Bound corrs
        if (first_rs < LEFT_CORR_BOUND) {
            first_rs = LEFT_CORR_BOUND;
        }

        if (second_rs < LEFT_CORR_BOUND) {
            second_rs = LEFT_CORR_BOUND;
        }
        
        if (first_rs > RIGHT_CORR_BOUND) {
            first_rs = RIGHT_CORR_BOUND;
        }

        if (second_rs > RIGHT_CORR_BOUND) {
            second_rs = RIGHT_CORR_BOUND;
        }

        float stat = (std::atanh(first_rs) -
                std::atanh(second_rs));    

        // This block is a copy of 
        // core.correlation_utils.pearson_std
        float first_ss = 1 / std::sqrt(first_size_ptr[ind] - 3);
        float second_ss = 1 / std::sqrt(second_size_ptr[ind] - 3);
        if (correlation == SPEARMAN) {
            first_ss *= std::sqrt(1.5);
            second_ss *= std::sqrt(1.5);
        }

        float std = std::sqrt(first_ss * first_ss +
                second_ss * second_ss);    
           
        if (pvalue_ptr != nullptr) {
            float pvalue = UNDEFINED_CORR_DIFF_TEST_VALUE;
            if (alternative == LESS) {
                pvalue = norm_cdf(stat / std);    
            } else if (alternative == GREATER) {
                pvalue = 1 - norm_cdf(stat / std);
            } else if (alternative == TWO_SIDED) {
                pvalue = 2 * norm_cdf(-std::abs(stat) / std);            
            }
            
            pvalue_ptr[ind] = pvalue;
        }

        stat_ptr[ind] = stat;
    }

    return 0;
}
