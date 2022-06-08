#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cmath>
#include <thread>
#include <string>
#include <utility>
#include <queue>
#include <cmath>
#include <algorithm>

#include "../utils/utils.h"

namespace py = pybind11;

using NumPyFloatArray = py::array_t<float, py::array::c_style>;
using NumPyDoubleArray = py::array_t<double, py::array::c_style>;
using NumPyIntArray = py::array_t<int32_t, py::array::c_style>;


NumPyIntArray unary_matrix(
    const NumPyIntArray &index,
    int base
) {
    /* This function is a vector analouge of
     * "unary_vector" */

    py::buffer_info index_buf = index.request(); 
    int index_size = index_buf.shape[0];
    int *index_ptr = (int *) index_buf.ptr;

    NumPyIntArray paired_array = NumPyIntArray(index_size * base);
    int *pa_ptr = (int *) paired_array.request().ptr;
    
    for (int i = 0; i < index_size; ++i){
        for (int j = 0; j < base; ++j) {
            pa_ptr[i * base + j] = unary_index(index_ptr[i], j, base);    
        }
    }
    
    paired_array.resize({index_size, base});
    return paired_array;
}

NumPyFloatArray quadrate(
    const NumPyFloatArray &flatten_array,
    const NumPyIntArray &index,
    int base
) {
    /* Exrtacts all elemnts of the alphabetically sorted
     * "flatten_array" with indexes in "index". After that
     * quadrate the extracted values */
    
    py::buffer_info index_buf = index.request(); 
    int index_size = index_buf.shape[0];
    int *index_ptr = (int *) index_buf.ptr;
    
    float *array_ptr = (float *) flatten_array.request().ptr;

    NumPyFloatArray matrix = NumPyFloatArray(index_size * base);
    float *matrix_ptr = (float *) matrix.request().ptr;
    
    for (int i = 0; i < index_size; ++i){
        for (int j = 0; j < base; ++j) {
            if (index_ptr[i] == j) {
                matrix_ptr[i * base + j] = 1.;
                continue;
            }
                
            matrix_ptr[i * base + j] = array_ptr[
                unary_index(index_ptr[i], j, base)
            ];    
        }
    }
    
    matrix.resize({index_size, base});
    return matrix;
}

int _reorder(
    NumPyIntArray &source_indexes,
    NumPyIntArray &target_indexes
) {
    py::buffer_info buf = 
        source_indexes.request();
    int size = buf.shape[0];

    int *source_ind_ptr =
        (int *) source_indexes.request().ptr;
    int *target_ind_ptr =
        (int *) target_indexes.request().ptr;
    
    reorder(
        source_ind_ptr,
        target_ind_ptr,
        size,
        (float *) nullptr
    );
    return 0;
}

int __reorder(
    NumPyIntArray &source_indexes,
    NumPyIntArray &target_indexes,
    NumPyFloatArray &data
) {
    py::buffer_info data_buf = data.request();
    int size = data_buf.shape[0];

    float *data_ptr = (float *) data_buf.ptr;
    int *source_ind_ptr =
        (int *) source_indexes.request().ptr;
    int *target_ind_ptr =
        (int *) target_indexes.request().ptr;
    
    reorder(
        source_ind_ptr,
        target_ind_ptr,
        size,
        data_ptr
    );

    return 0;
}

PYBIND11_MODULE(utils, m) {
    m.def("_unary_index", &unary_index);
    m.def("_paired_index", &paired_index);
    m.def("_unary_vector", &unary_vector);
    m.def("_unary_matrix", &unary_matrix);
    m.def("_quadrate", &quadrate);
    m.def("_reorder", &_reorder);
    m.def("_reorder_data", &__reorder);
    m.attr("UNDEFINED_INDEX") = py::int_(UNDEFINED_INDEX);
}
