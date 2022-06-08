#include <cmath>
#include <string>
#include <utility>
#include <thread>
#include <algorithm>
#include <queue>

#include "utils.h"

int unary_index(int first, int second, int base) {
    /* This funcito computes an indexed of
     * a pair ("first", "second") in the alphabetically
     * ordered raw that was made by all unique pairs
     * from "0" to "base" */

    if (first == second) {
        return UNDEFINED_INDEX;
    }
    
    if (first > second) {
        int tmp = second;
        second = first;
        first = tmp;
    }

    int unary_index = (2 * base - first - 1) * (first) / 2;
    unary_index += second - first - 1;

    return unary_index;
}

std::pair<int, int> paired_index(int index, int base) {
    /* This function is inverse to "unary_index" */

    int i = std::floor(((2 * base - 1) - std::sqrt((2 * base - 1) * 
            (2 * base - 1) - 8 * index)) / 2);
    int j = (index % base + ((i + 2) * (i + 1) / 2) % base) % base; 
    
    return std::pair<int, int>(i, j);
}

std::vector<int> unary_vector(int index, int base) {
    /* Computes all indexes in the alphabetically
     * ordered raw with "index" in the first position
     * of the paired index */

    std::vector<int> paired_array(base);
    for (int j = 0; j < base; ++j) {
        paired_array[j] = unary_index(index, j, base);    
    }
    
    return paired_array;
}

int inverse(
    int *arr,
    float *reverse,
    int *index_ptr,
    int size
) {
    int index;
    for (int i = 0; i < size; ++i) {
        if (!index_ptr) {
            index = arr[i];
        } else {
            index = index_ptr[arr[i]];
        }

        reverse[index] = (float) i;
    }

    return 0;
}    

int range(int *arr, int size) {
    for (int i = 0; i < size; ++i) {
        arr[i] = i;
    }

    return 0;
}

int _rank_data(
    float *data_ptr,
    float *rank_ptr,
    int sample_size,
    int start_ind, 
    int end_ind,
    int *sample_ind_ptr,
    int sample_ind_size
) {
    if (!sample_ind_ptr) {
        sample_ind_size = sample_size;
    }
    
    std::vector<int> indexes(sample_ind_size);

    for (int i = start_ind; i < end_ind; ++i) {
        range(indexes.data(), sample_ind_size);
        std::sort(
            indexes.begin(),
            indexes.end(),
            [
                data_ptr,
                i,
                sample_size,
                sample_ind_ptr
            ](int i1, int i2) {
                if (!sample_ind_ptr) {
                    i1 = sample_size * i + i1;
                    i2 = sample_size * i + i2;
                } else {
                    i1 = sample_size * i + sample_ind_ptr[i1];
                    i2 = sample_size * i + sample_ind_ptr[i2];
                }

                return data_ptr[i1] < data_ptr[i2];
            }
        );

        inverse(
            indexes.data(),
            rank_ptr + sample_size * i,
            sample_ind_ptr,
            sample_ind_size
        );
    }

    return 0;
}

int rank_data(
    float *data_ptr,
    float *rank_ptr,
    int sample_size,
    int index_size,
    int *sample_ind_ptr,
    int sample_ind_size,
    int process_num
) {
    if (!sample_ind_ptr) {
        sample_ind_size = sample_size;
    }

    std::queue<std::thread> threads;
    int batch_size = index_size / process_num;
    for (int i = 0; i < process_num; ++i) {
        int left_border = i * batch_size;
        int right_border = (i + 1) * batch_size;
        if (i == process_num - 1) {
            right_border = index_size;
        }

        std::thread thr(_rank_data,
            data_ptr,
            rank_ptr,
            sample_size,
            left_border,
            right_border,
            sample_ind_ptr,
            sample_ind_size
        );

        threads.push(move(thr));
    }

    while (!threads.empty()) {
        threads.front().join();
        threads.pop();
    }

    return 0;
}

int swap(
    int *data,
    int *indexes,
    int size
) {
    std::vector<int> swaper(size);
    for (int i = 0; i < size; ++i) {
        swaper[i] = data[indexes[i]];
    }

    for (int i = 0; i < size; ++i) {
        data[i] = swaper[i];
    }

    return 0;
}

int swap(
    float *data,
    int *indexes,
    int size
) {
    std::vector<float> swaper(size);
    for (int i = 0; i < size; ++i) {
        swaper[i] = data[indexes[i]];
    }

    for (int i = 0; i < size; ++i) {
        data[i] = swaper[i];
    }

    return 0;
}

int reorder(
    int *source_ind_ptr,
    int *target_ind_ptr,
    int size,
    float *data_ptr
) { 
    std::vector<int> indexes(size);
    for (int i = 0; i < size; ++i) {
        indexes[i] = i;
    }
    
    std::sort(
        indexes.begin(),
        indexes.end(),
        [
            source_ind_ptr
        ](int i1, int i2) {
            return source_ind_ptr[i1] < source_ind_ptr[i2];
        }
    );

    swap(source_ind_ptr, indexes.data(), size);
    swap(target_ind_ptr, indexes.data(), size);
    
    if (data_ptr) {
        swap(source_ind_ptr, indexes.data(), size);
    }
    
    return 0;
}
