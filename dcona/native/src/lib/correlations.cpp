#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cmath>
#include <thread>
#include <string>
#include <utility>
#include <queue>
#include <tuple>
#include <map>
#include <algorithm>
#include <random>
#include <vector>

#include "../correlations/correlations.h"
#include "../utils/utils.h"

namespace py = pybind11;

using NumPyFloatArray = py::array_t<float, py::array::c_style>;
using NumPyDoubleArray = py::array_t<double, py::array::c_style>;
using NumPyIntArray = py::array_t<int32_t, py::array::c_style>;


NumPyFloatArray correlation_indexed(
    const NumPyFloatArray &data,
    const NumPyIntArray &source_indexes,
    const NumPyIntArray &target_indexes,
    const std::string correlation=SPEARMAN,
    int process_num=1
) {
    py::buffer_info source_ind_buf = source_indexes.request();
    py::buffer_info target_ind_buf = target_indexes.request();
    py::buffer_info data_buf = data.request();
    
    int sample_size = data_buf.shape[1]; 
    int index_size = source_ind_buf.shape[0];
    if (source_ind_buf.size != target_ind_buf.size) {
        throw std::runtime_error("Index shapes must match");
    }

    if (process_num > index_size) {
        process_num = index_size;
    }

    if (process_num <= 0) {    
        throw std::runtime_error("Process number error");
    }
    
    NumPyFloatArray corrs = NumPyFloatArray(source_ind_buf.size);
    float *corrs_ptr = (float *) corrs.request().ptr;
    
    std::queue<std::thread> threads;
    int batch_size = index_size / process_num;
    for (int i = 0; i < process_num; ++i) {
        int left_border = i * batch_size;
        int right_border = (i + 1) * batch_size;
        if (i == process_num - 1) {
            right_border = index_size;
        }

        if (correlation == SPEARMAN) {
            std::thread thr(spearmanr,
                (float *) data_buf.ptr,
                sample_size,
                (int *) source_ind_buf.ptr,
                (int *) target_ind_buf.ptr,
                corrs_ptr,
                left_border,
                right_border,
                index_size,
                (int *) nullptr,
                -1
            );
            
            threads.push(move(thr));
        } else {
            std::thread thr(pearsonr,
                (float *) data_buf.ptr,
                sample_size,
                (int *) source_ind_buf.ptr,
                (int *) target_ind_buf.ptr,
                corrs_ptr,
                left_border,
                right_border,
                index_size,
                (int *) nullptr,
                -1
            );

            threads.push(move(thr));
        }    
    }

    while (!threads.empty()) {
        threads.front().join();
        threads.pop();
    }

    return corrs;
}

NumPyFloatArray correlation_exhaustive(
    const NumPyFloatArray &data,
    std::string correlation=SPEARMAN,
    int process_num=1
) {    
    py::buffer_info data_buf = data.request();
    
    int sample_size = data_buf.shape[1]; 
    int index_size = data_buf.shape[0];
    int pairs_num = index_size * (index_size - 1) / 2;

    if (process_num > index_size) {
        process_num = index_size;
    }

    if (process_num <= 0) {    
        throw std::runtime_error("Process number error");
    }
    
    NumPyFloatArray corrs = NumPyFloatArray(pairs_num);
    float *corrs_ptr = (float *) corrs.request().ptr;
    
    std::queue<std::thread> threads;
    int batch_size = pairs_num / process_num;
    for (int i = 0; i < process_num; ++i) {
        int left_border = i * batch_size;
        int right_border = (i + 1) * batch_size;
        if (i == process_num - 1) {
            right_border = pairs_num;
        }
        
        if (correlation == SPEARMAN) {
            std::thread thr(spearmanr,
                (float *) data_buf.ptr,
                sample_size,
                (int *) nullptr,
                (int *) nullptr,
                corrs_ptr,
                left_border,
                right_border,
                index_size,
                (int *) nullptr,
                -1
            );
            
            threads.push(move(thr));
        } else {
            std::thread thr(pearsonr,
                (float *) data_buf.ptr,
                sample_size,
                (int *) nullptr,
                (int *) nullptr,
                corrs_ptr,
                left_border,
                right_border,
                index_size,
                (int *) nullptr,
                -1
            );
            
            threads.push(move(thr));

        }
    }

    while (!threads.empty()) {
        threads.front().join();
        threads.pop();
    }

    return corrs;
}

PYBIND11_MODULE(correlations, m) {
    m.def("_correlation_indexed", &correlation_indexed);
    m.def("_correlation_exhaustive", &correlation_exhaustive);
    m.attr("UNDEFINED_CORR_VALUE") = py::float_(UNDEFINED_CORR_VALUE);
}
