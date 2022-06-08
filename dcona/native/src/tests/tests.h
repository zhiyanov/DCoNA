#ifndef TESTS_H
#define TESTS_H

#include <string>

#include "../correlations/correlations.h"

const float UNDEFINED_CORR_DIFF_TEST_VALUE = -2;

const float LEFT_CORR_BOUND = -0.99;
const float RIGHT_CORR_BOUND = 0.99;

const std::string TWO_SIDED = "two-sided";
const std::string LESS = "less";
const std::string GREATER = "greater";


int ztest_unsized(
    float *first_rs_ptr, int first_size,
    float *second_rs_ptr, int second_size,
    float *stat_ptr, float *pvalue_ptr,
    int start_ind, int end_ind,
    const std::string correlation=SPEARMAN,
    const std::string alternative=TWO_SIDED
);

int ztest_sized(
    float *first_rs_ptr, int *first_size_ptr,
    float *second_rs_ptr, int *second_size_ptr,
    float *stat_ptr, float *pvalue_ptr,
    int start_ind, int end_ind,
    const std::string correlation=SPEARMAN,
    const std::string alternative=TWO_SIDED
);

#endif
