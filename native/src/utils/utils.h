#ifndef UTILS_H
#define UTILS_H

#include <utility>
#include <vector>

const int UNDEFINED_INDEX = -1;


int unary_index(int first, int second, int base);

std::pair<int, int> paired_index(int index, int base);

std::vector<int> unary_vector(int index, int base);

int range(int *arr, int size);

int rank_data(
    float *data_ptr,
    float *rank_ptr,
    int sample_size,
    int index_size,
    int *sample_ind_ptr=nullptr,
    int sample_ind_size=-1,
    int process_num=1
);

int reorder(
    int *source_ind_ptr,
    int *target_ind_ptr,
    int size,
    float *data_ptr=nullptr
);

#endif
