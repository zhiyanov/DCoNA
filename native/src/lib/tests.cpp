#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cmath>
#include <thread>
#include <string>
#include <utility>
#include <queue>
#include <tuple>

#include <algorithm>
#include <random>
#include <vector>

#include "../tests/tests.h"

namespace py = pybind11;

using NumPyFloatArray = py::array_t<float, py::array::c_style>;
using NumPyDoubleArray = py::array_t<double, py::array::c_style>;
using NumPyIntArray = py::array_t<int32_t, py::array::c_style>;


std::pair<NumPyFloatArray, NumPyFloatArray> _ztest_sized(
    const NumPyFloatArray &first_rs,
    const NumPyIntArray &first_size,
    const NumPyFloatArray &second_rs,
    const NumPyIntArray &second_size,
    const std::string correlation=SPEARMAN,
    const std::string alternative=TWO_SIDED,
    int process_num=1
) {  
    py::buffer_info first_rs_buf = first_rs.request();
    py::buffer_info second_rs_buf = second_rs.request();
    
    int rs_number = first_rs_buf.shape[0];
    if (first_rs_buf.size != second_rs_buf.size) {
        throw std::runtime_error("Correlation shapes must match");
    }

    if (process_num > rs_number) {
        process_num = rs_number;
    }

    if (process_num <= 0) {    
        throw std::runtime_error("Process number error");
    }
    
    NumPyFloatArray stat = NumPyFloatArray(first_rs_buf.size);
    NumPyFloatArray pvalue = NumPyFloatArray(first_rs_buf.size);
    
    float *stat_ptr = (float *) stat.request().ptr;
    float *pvalue_ptr = (float *) pvalue.request().ptr;
    
    std::queue<std::thread> threads;
    int batch_size = rs_number / process_num;
    for (int i = 0; i < process_num; ++i) {
        int left_border = i * batch_size;
        int right_border = (i + 1) * batch_size;
        if (i == process_num - 1) {
            right_border = rs_number;
        }
        
        std::thread thr(ztest_sized,
            (float *) first_rs_buf.ptr,
            (int *) first_size.request().ptr,
            (float *) second_rs_buf.ptr,
            (int *) second_size.request().ptr,
            stat_ptr,
            pvalue_ptr,
            left_border,
            right_border,
            correlation,
            alternative
        );
            
        threads.push(move(thr));
    }

    while (!threads.empty()) {
        threads.front().join();
        threads.pop();
    }

    return std::pair<NumPyFloatArray,
           NumPyFloatArray>(stat, pvalue);
}

std::pair<NumPyFloatArray, NumPyFloatArray> _ztest_unsized(
    const NumPyFloatArray &first_rs,
    int first_size,
    const NumPyFloatArray &second_rs,
    int second_size,
    const std::string correlation=SPEARMAN,
    const std::string alternative=TWO_SIDED,
    int process_num=1
) {  
    py::buffer_info first_rs_buf = first_rs.request();
    py::buffer_info second_rs_buf = second_rs.request();
    
    int rs_number = first_rs_buf.shape[0];
    if (first_rs_buf.size != second_rs_buf.size) {
        throw std::runtime_error("Correlation shapes must match");
    }

    if (process_num > rs_number) {
        process_num = rs_number;
    }

    if (process_num <= 0) {    
        throw std::runtime_error("Process number error");
    }
    
    NumPyFloatArray stat = NumPyFloatArray(first_rs_buf.size);
    NumPyFloatArray pvalue = NumPyFloatArray(first_rs_buf.size);
    
    float *stat_ptr = (float *) stat.request().ptr;
    float *pvalue_ptr = (float *) pvalue.request().ptr;
    
    std::queue<std::thread> threads;
    int batch_size = rs_number / process_num;
    for (int i = 0; i < process_num; ++i) {
        int left_border = i * batch_size;
        int right_border = (i + 1) * batch_size;
        if (i == process_num - 1) {
            right_border = rs_number;
        }
        
        std::thread thr(ztest_unsized,
            (float *) first_rs_buf.ptr,
            first_size,
            (float *) second_rs_buf.ptr,
            second_size,
            stat_ptr,
            pvalue_ptr,
            left_border,
            right_border,
            correlation,
            alternative
        );
            
        threads.push(move(thr));
    }

    while (!threads.empty()) {
        threads.front().join();
        threads.pop();
    }

    return std::pair<
        NumPyFloatArray,
        NumPyFloatArray
    >(stat, pvalue);
}

PYBIND11_MODULE(tests, m) {
    m.def("_ztest_sized", &_ztest_sized);
    m.def("_ztest_unsized", &_ztest_unsized);
}
