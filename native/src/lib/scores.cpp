
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

#include "../scores/scores.h"
#include "../utils/utils.h"

namespace py = pybind11;

using NumPyFloatArray = py::array_t<float, py::array::c_style>;
using NumPyDoubleArray = py::array_t<double, py::array::c_style>;
using NumPyIntArray = py::array_t<int32_t, py::array::c_style>;

const int REPEATS_NUMBER = 1000;


std::pair<
    NumPyFloatArray,
    NumPyFloatArray
> score_indexed(
    const NumPyFloatArray &data,
    const NumPyIntArray &source_indexes,
    const std::string score=MEAN,
    int process_num=1
) {
    py::buffer_info source_ind_buf = source_indexes.request();
    py::buffer_info data_buf = data.request();

    float *data_ptr = (float *) data_buf.ptr;
    int *source_ind_ptr = (int *) source_ind_buf.ptr;
    
    int index_size = source_ind_buf.shape[0];

    if (process_num > index_size) {
        process_num = index_size;
    }

    if (process_num <= 0) {    
        throw std::runtime_error("Process number error");
    }
    
    std::vector<int> _sources;
    std::vector<int> starts;
    std::vector<int> ends;

    for (int i = 0; i < index_size;) {
        if (!_sources.size() ||
                _sources.back() < source_ind_ptr[i]) {
            int start = i;
            int end = i;
            
            while (end < index_size &&
                source_ind_ptr[start] == source_ind_ptr[end]) {
                ++end;
            }
            
            _sources.push_back(source_ind_ptr[start]);
            starts.push_back(start);
            ends.push_back(end);
            i = end;
        }
    }

    int sources_size = _sources.size(); 
    NumPyFloatArray sources = NumPyFloatArray(sources_size);
    float *sources_ptr = (float *) sources.request().ptr;
    for (int j = 0; j < sources_size; ++j) {
        sources_ptr[j] = _sources[j]; 
    }
    
    NumPyFloatArray scores = NumPyFloatArray(sources_size);
    float *scores_ptr = (float *) scores.request().ptr;
    
    if (process_num > sources_size) {
        process_num = sources_size;
    }

    if (process_num <= 0) {    
        throw std::runtime_error("Process number error");
    }

    std::queue<std::thread> threads;
    int batch_size = sources_size / process_num;
    for (int i = 0; i < process_num; ++i) {
        int left_border = i * batch_size;
        int right_border = (i + 1) * batch_size;
        if (i == process_num - 1) {
            right_border = starts.size();
        }

        if (score == MEAN) {
            std::thread thr(_mean,
                data_ptr,
                starts.data(),
                ends.data(),
                left_border,
                right_border,
                scores_ptr,
                false
            );
            
            threads.push(move(thr));
        } else if (score == MEDIAN) {
            std::thread thr(_quantile,
                data_ptr,
                starts.data(),
                ends.data(),
                left_border,
                right_border,
                scores_ptr,
                QMEDIAN,
                false
            );

            threads.push(move(thr));
        }    
    }

    while (!threads.empty()) {
        threads.front().join();
        threads.pop();
    }

    return std::pair<
        NumPyFloatArray,
        NumPyFloatArray
    >(
        sources, scores    
    );
}

NumPyFloatArray score_exhaustive(
    const NumPyFloatArray &data,
    int sources_size,
    const std::string score=MEAN,
    int process_num=1
) {
    py::buffer_info data_buf = data.request();
    float *data_ptr = (float *) data_buf.ptr;

    if (process_num > sources_size) {
        process_num = sources_size;
    }

    if (process_num <= 0) {    
        throw std::runtime_error("Process number error");
    } 
    
    NumPyFloatArray scores = NumPyFloatArray(sources_size);
    float *scores_ptr = (float *) scores.request().ptr;
    
    if (process_num > sources_size) {
        process_num = sources_size;
    }

    if (process_num <= 0) {    
        throw std::runtime_error("Process number error");
    }

    std::queue<std::thread> threads;
    int batch_size = sources_size / process_num;
    for (int i = 0; i < process_num; ++i) {
        int left_border = i * batch_size;
        int right_border = (i + 1) * batch_size;
        if (i == process_num - 1) {
            right_border = sources_size;
        }

        if (score == MEAN) {
            std::thread thr(__mean,
                data_ptr,
                sources_size,
                left_border,
                right_border,
                scores_ptr,
                false
            );
            
            threads.push(move(thr));
        } else if (score == MEDIAN) {
            std::thread thr(__quantile,
                data_ptr,
                sources_size,
                left_border,
                right_border,
                scores_ptr,
                QMEDIAN,
                false
            );

            threads.push(move(thr));
        }    
    }

    while (!threads.empty()) {
        threads.front().join();
        threads.pop();
    }
    
    return scores;
}

PYBIND11_MODULE(scores, m) {
    m.def("_score_indexed", &score_indexed);
    m.def("_score_exhaustive", &score_exhaustive);
}
