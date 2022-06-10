#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <iostream>
#include <cmath>
#include <thread>
#include <string>
#include <utility>
#include <queue>
#include <tuple>
#include <algorithm>
#include <random>
#include <vector>

#include "../pipelines/pipelines.h"
#include "../utils/utils.h"
// #include "../tqdm/tqdm.h"
#include "../tqdm/progressbar.h"

namespace py = pybind11;

using NumPyFloatArray = py::array_t<float, py::array::c_style>;
using NumPyDoubleArray = py::array_t<double, py::array::c_style>;
using NumPyIntArray = py::array_t<int32_t, py::array::c_style>;

const int SEED = 733;

std::tuple<
    NumPyFloatArray,
    NumPyFloatArray,
    NumPyFloatArray,
    NumPyFloatArray,
    NumPyFloatArray
> ztest_pipeline_indexed(
    const NumPyFloatArray &data,
    const NumPyIntArray &source_indexes,
    const NumPyIntArray &target_indexes,
    const NumPyIntArray &reference_indexes,
    const NumPyIntArray &experimental_indexes,
    const std::string correlation=SPEARMAN,
    const std::string alternative=TWO_SIDED,
    int repeats_number=REPEATS_NUMBER,
    int process_num=1    
) {
    py::buffer_info data_buf = data.request();
    float *data_ptr = (float *) data_buf.ptr;
    int data_len    = data_buf.shape[0];
    int sample_size = data_buf.shape[1];

    py::buffer_info source_ind_buf = source_indexes.request();
    int *source_ind_ptr = (int *) source_ind_buf.ptr;
    
    py::buffer_info target_ind_buf = target_indexes.request();
    int *target_ind_ptr = (int *) target_ind_buf.ptr;
    
    int index_size = source_ind_buf.shape[0];
    
    py::buffer_info ref_ind_buf = reference_indexes.request();
    int ref_ind_size = ref_ind_buf.shape[0];
    int *ref_ind_ptr = (int *) ref_ind_buf.ptr;
    
    py::buffer_info exp_ind_buf = experimental_indexes.request();
    int exp_ind_size = exp_ind_buf.shape[0];
    int *exp_ind_ptr = (int *) exp_ind_buf.ptr;

    // Real data    
    NumPyFloatArray ref_corrs = NumPyFloatArray(index_size);
    float *ref_corrs_ptr = (float *) ref_corrs.request().ptr;
    
    NumPyFloatArray exp_corrs = NumPyFloatArray(index_size);
    float *exp_corrs_ptr = (float *) exp_corrs.request().ptr;

    NumPyFloatArray stat = NumPyFloatArray(index_size);
    float *stat_ptr = (float *) stat.request().ptr;
    
    NumPyFloatArray pvalue = NumPyFloatArray(index_size);
    float *pvalue_ptr = (float *) pvalue.request().ptr;
    
    NumPyFloatArray boot_pvalue = NumPyFloatArray(index_size);
    float *boot_pvalue_ptr = (float *) boot_pvalue.request().ptr;
    
    // Rank data    
    float *rank_ptr;
    if (correlation == SPEARMAN) {
        rank_ptr = new float[
            data_len * sample_size
        ];
    }
    
    // Bootstrapped data
    float *boot_ref_corrs_ptr =  new float[index_size];
    float *boot_exp_corrs_ptr = new float[index_size];
    float *boot_stat_ptr = new float[index_size];
     
    int *boot_ref_ind_ptr = new int[ref_ind_size];
    int *boot_exp_ind_ptr = new int[exp_ind_size];

    // Bootstrap indexes initialization
    std::vector<int> indexes(sample_size);
    for (int i = 0; i < sample_size; ++i) {
        indexes[i] = i;
    }

    // Random generator initialization
    // std::random_device random_dev;
    std::mt19937 random_gen(SEED);

    // Bootstrap pvalue computations    
    float *dpr, *rcp, *ecp, *sp, *pp;
    int *rip, *eip;
    
    std::cout << "Permutation progress: ";
    progressbar bar(repeats_number + 1);
    for (int r = 0; r < repeats_number + 1; ++r) {
        if (PyErr_CheckSignals() != 0) {
            throw py::error_already_set();
        }

        if (r == 0) {
            rcp = ref_corrs_ptr;
            ecp = exp_corrs_ptr;
            sp  = stat_ptr;
            pp  = pvalue_ptr;

            rip = ref_ind_ptr;
            eip = exp_ind_ptr;
        } else {
            std::shuffle(indexes.begin(), indexes.end(), random_gen);
            for (int i = 0; i < ref_ind_size; ++i) {
                boot_ref_ind_ptr[i] = indexes[i];
            }
            for (int i = 0; i < exp_ind_size; ++i) {
                boot_exp_ind_ptr[i] = indexes[ref_ind_size + i];
            }

            rcp = boot_ref_corrs_ptr;
            ecp = boot_exp_corrs_ptr;
            sp  = boot_stat_ptr;
            pp  = nullptr;
            
            rip = boot_ref_ind_ptr;
            eip = boot_exp_ind_ptr;
        }

        if (correlation == SPEARMAN) {
            rank_data(
                data_ptr,
                rank_ptr,
                sample_size,
                data_len,
                rip,
                ref_ind_size,
                process_num
            );

            rank_data(
                data_ptr,
                rank_ptr,
                sample_size,
                data_len,
                eip,
                exp_ind_size,
                process_num
            );

            dpr = rank_ptr;
        } else {
            dpr = data_ptr;
        }
        
        std::queue<std::thread> threads;
        int batch_size = index_size / process_num;
        for (int i = 0; i < process_num; ++i) {
            int left_border = i * batch_size;
            int right_border = (i + 1) * batch_size;
            if (i == process_num - 1) {
                right_border = index_size;
            }
            
            std::thread thr(ztest_pipeline,
                dpr,
                sample_size,
                source_ind_ptr,
                target_ind_ptr,
                left_border,
                right_border,
                index_size,
                rip,
                eip,
                ref_ind_size,
                exp_ind_size,
                rcp,
                ecp,
                sp,
                pp,
                correlation,
                alternative
            );
            
            threads.push(move(thr));
        }

        while (!threads.empty()) {
            threads.front().join();
            threads.pop();
        }

        if (r > 0) {
            for (int i = 0; i < index_size; ++i) {
                if ((alternative == TWO_SIDED) &&
                        (std::abs(stat_ptr[i]) <= std::abs(boot_stat_ptr[i]))) {
                    boot_pvalue_ptr[i] += 1;
                }

                if ((alternative == LESS) &&
                        (stat_ptr[i] <= boot_stat_ptr[i])) {
                    boot_pvalue_ptr[i] += 1;
                }
                
                if ((alternative == GREATER) &&
                        (stat_ptr[i] >= boot_stat_ptr[i])) {
                    boot_pvalue_ptr[i] += 1;
                }
            }
        } else {
            for (int i = 0; i < index_size; ++i) {
                boot_pvalue_ptr[i] = 0;
            }
        }

        bar.update();
    }
    std::cout << "\n";
    
    if (repeats_number > 0) {
        for (int i = 0; i < index_size; ++i) {
            boot_pvalue_ptr[i] /= repeats_number;
        }
    }
    
    if (correlation == SPEARMAN) {
        delete[] rank_ptr;
    }

    delete[] boot_ref_corrs_ptr;
    delete[] boot_exp_corrs_ptr;
    delete[] boot_stat_ptr;
     
    delete[] boot_ref_ind_ptr;
    delete[] boot_exp_ind_ptr;
    
    return std::tuple<
        NumPyFloatArray, NumPyFloatArray, NumPyFloatArray,
        NumPyFloatArray, NumPyFloatArray
    >(
        ref_corrs, exp_corrs, stat,
        pvalue, boot_pvalue
    );
}

std::tuple<
    NumPyIntArray, NumPyFloatArray,
    NumPyFloatArray
> _score_pipeline_indexed(
    const NumPyFloatArray &data,
    const NumPyIntArray &source_indexes,
    const NumPyIntArray &target_indexes,
    const NumPyIntArray &reference_indexes,
    const NumPyIntArray &experimental_indexes,
    const std::string correlation=SPEARMAN,
    const std::string score=MEAN,
    const std::string alternative=TWO_SIDED,
    int repeats_number=REPEATS_NUMBER,
    int process_num=1    
) {
    py::buffer_info data_buf = data.request();
    float *data_ptr = (float *) data_buf.ptr;
    int data_len    = data_buf.shape[0];
    int sample_size = data_buf.shape[1];

    py::buffer_info source_ind_buf = source_indexes.request();
    int *source_ind_ptr = (int *) source_ind_buf.ptr;
    
    py::buffer_info target_ind_buf = target_indexes.request();
    int *target_ind_ptr = (int *) target_ind_buf.ptr;
    
    int index_size = source_ind_buf.shape[0];
    
    py::buffer_info ref_ind_buf = reference_indexes.request();
    int ref_ind_size = ref_ind_buf.shape[0];
    int *ref_ind_ptr = (int *) ref_ind_buf.ptr;
    
    py::buffer_info exp_ind_buf = experimental_indexes.request();
    int exp_ind_size = exp_ind_buf.shape[0];
    int *exp_ind_ptr = (int *) exp_ind_buf.ptr;

    // Rank data    
    float *rank_ptr;
    if (correlation == SPEARMAN) {
        rank_ptr = new float[
            data_len * sample_size
        ];
    }
    
    // Bootstrapped data
    float *boot_ref_corrs_ptr =  new float[index_size];
    float *boot_exp_corrs_ptr = new float[index_size];
    float *boot_stat_ptr = new float[index_size];
     
    int *boot_ref_ind_ptr = new int[ref_ind_size];
    int *boot_exp_ind_ptr = new int[exp_ind_size];

    // Bootstrap indexes initialization
    std::vector<int> indexes(sample_size);
    for (int i = 0; i < sample_size; ++i) {
        indexes[i] = i;
    }

    // Random generator initialization
    // std::random_device random_dev;
    std::mt19937 random_gen(SEED);

    // Scores initialization
    std::vector<int> _sources;
    std::vector<int> starts;
    std::vector<int> ends;

    for (int i = 0; i < index_size;) {
        if ((!_sources.size()) ||
                (_sources.back() < source_ind_ptr[i])) {
            int start = i;
            int end = i;
            
            while ((end < index_size) &&
                    (source_ind_ptr[start] == source_ind_ptr[end])) {
                ++end;
            }

            _sources.push_back(source_ind_ptr[start]);
            starts.push_back(start);
            ends.push_back(end);
            
            i = end;
        } else {
            ++i;
        }
    }
    
    // Scores data
    int sources_size = _sources.size();

    NumPyIntArray sources = NumPyIntArray(sources_size);
    int *sources_ptr = (int *) sources.request().ptr;

    NumPyFloatArray scores = NumPyFloatArray(sources_size);
    float *scores_ptr = (float *) scores.request().ptr;

    NumPyFloatArray pvalues = NumPyFloatArray(sources_size);
    float *pvalues_ptr = (float *) pvalues.request().ptr;

    for (int j = 0; j < sources_size; ++j) {
        sources_ptr[j] = _sources[j];
        pvalues_ptr[j] = 0;
    }

    // Bootstrap scores
    float *boot_scores_ptr = new float[sources_size];
    
    if (process_num > sources_size) {
        process_num = sources_size;
    }

    if (process_num <= 0) {    
        throw std::runtime_error("Process number error");
    }
    
    // Bootstrap pvalue computations    
    float *dpr, *rcp, *ecp, *sp, *pp;
    int *rip, *eip;
    float *scp;
    
    std::cout << "Permutation progress: ";
    progressbar bar(repeats_number + 1);
    for (int r = 0; r < repeats_number + 1; ++r) {
        if (PyErr_CheckSignals() != 0) {
            throw py::error_already_set();
        }

        if (r == 0) {
            rcp = boot_ref_corrs_ptr;
            ecp = boot_exp_corrs_ptr;
            sp  = boot_stat_ptr;
            pp  = nullptr;

            rip = ref_ind_ptr;
            eip = exp_ind_ptr;

            scp = scores_ptr;
        } else {
            std::shuffle(indexes.begin(), indexes.end(), random_gen);
            for (int i = 0; i < ref_ind_size; ++i) {
                boot_ref_ind_ptr[i] = indexes[i];
            }
            for (int i = 0; i < exp_ind_size; ++i) {
                boot_exp_ind_ptr[i] = indexes[ref_ind_size + i];
            }

            rcp = boot_ref_corrs_ptr;
            ecp = boot_exp_corrs_ptr;
            sp  = boot_stat_ptr;
            pp  = nullptr;
            
            rip = boot_ref_ind_ptr;
            eip = boot_exp_ind_ptr;

            scp = boot_scores_ptr;
        }
        
        if (correlation == SPEARMAN) {
            rank_data(
                data_ptr,
                rank_ptr,
                sample_size,
                data_len,
                rip,
                ref_ind_size,
                process_num
            );

            rank_data(
                data_ptr,
                rank_ptr,
                sample_size,
                data_len,
                eip,
                exp_ind_size,
                process_num
            );

            dpr = rank_ptr;
        } else {
            dpr = data_ptr;
        }
        
        std::queue<std::thread> threads;
        int batch_size = sources_size / process_num;
        for (int i = 0; i < process_num; ++i) {
            int left_border = i * batch_size;
            int right_border = (i + 1) * batch_size;
            if (i == process_num - 1) {
                right_border = sources_size;
            }
            
            std::thread thr(score_pipeline_indexed,
                dpr,
                sample_size,
                source_ind_ptr,
                target_ind_ptr,
                index_size,
                rip,
                eip,
                ref_ind_size,
                exp_ind_size,
                rcp,
                ecp,
                sp,
                pp,
                starts.data(),
                ends.data(),
                left_border,
                right_border,
                scp,
                correlation,
                alternative,
                score
            );
            
            threads.push(move(thr));
        }

        while (!threads.empty()) {
            threads.front().join();
            threads.pop();
        }

        if (r > 0) {
            for (int i = 0; i < sources_size; ++i) {
                if (alternative == TWO_SIDED && 
                        std::abs(scores_ptr[i]) <= std::abs(boot_scores_ptr[i])) {
                    pvalues_ptr[i] += 1;
                }

                if (alternative == LESS &&
                        scores_ptr[i] >= boot_scores_ptr[i]) {
                    pvalues_ptr[i] += 1;
                }


                if (alternative == GREATER &&
                        scores_ptr[i] <= boot_scores_ptr[i]) { 
                    pvalues_ptr[i] += 1;
                }
            }
        } else {
            for (int i = 0; i < sources_size; ++i) {
                pvalues_ptr[i] = 0;
            }
        }

        bar.update();
    }
    std::cout << "\n";

    if (repeats_number > 0) {
        for (int i = 0; i < sources_size; ++i) {
            pvalues_ptr[i] /= repeats_number;
        }
    }

    if (correlation == SPEARMAN) {
        delete[] rank_ptr;
    }

    delete[] boot_ref_corrs_ptr;
    delete[] boot_exp_corrs_ptr;
    delete[] boot_stat_ptr;
 
    delete[] boot_ref_ind_ptr;
    delete[] boot_exp_ind_ptr;

    delete[] boot_scores_ptr;
    
    return std::tuple<
        NumPyIntArray,
        NumPyFloatArray,
        NumPyFloatArray
    >(
        sources,
        scores,
        pvalues
    );
}

std::tuple<
    NumPyFloatArray,
    NumPyFloatArray,
    NumPyFloatArray,
    NumPyFloatArray,
    NumPyFloatArray
> ztest_pipeline_exhaustive(
    const NumPyFloatArray &data,
    const NumPyIntArray &reference_indexes,
    const NumPyIntArray &experimental_indexes,
    const std::string correlation=SPEARMAN,
    const std::string alternative=TWO_SIDED,
    int repeats_number=REPEATS_NUMBER,
    int process_num=1    
) {
    py::buffer_info data_buf = data.request();
    float *data_ptr = (float *) data_buf.ptr;
    int index_size = data_buf.shape[0];
    int sample_size = data_buf.shape[1];
    int pairs_num = index_size * (index_size - 1) / 2;

    py::buffer_info ref_ind_buf = reference_indexes.request();
    int ref_ind_size = ref_ind_buf.shape[0];
    int *ref_ind_ptr = (int *) ref_ind_buf.ptr;
    
    py::buffer_info exp_ind_buf = experimental_indexes.request();
    int exp_ind_size = exp_ind_buf.shape[0];
    int *exp_ind_ptr = (int *) exp_ind_buf.ptr;

    // Real data    
    NumPyFloatArray ref_corrs = NumPyFloatArray(pairs_num);
    float *ref_corrs_ptr = (float *) ref_corrs.request().ptr;
    
    NumPyFloatArray exp_corrs = NumPyFloatArray(pairs_num);
    float *exp_corrs_ptr = (float *) exp_corrs.request().ptr;

    NumPyFloatArray stat = NumPyFloatArray(pairs_num);
    float *stat_ptr = (float *) stat.request().ptr;
    
    NumPyFloatArray pvalue = NumPyFloatArray(pairs_num);
    float *pvalue_ptr = (float *) pvalue.request().ptr;
    
    NumPyFloatArray boot_pvalue = NumPyFloatArray(pairs_num);
    float *boot_pvalue_ptr = (float *) boot_pvalue.request().ptr;

    // Rank data    
    float *rank_ptr;
    if (correlation == SPEARMAN) {
        rank_ptr = new float[
            index_size * sample_size
        ];
    }

    // Bootstrapped data
    float *boot_ref_corrs_ptr =  new float[pairs_num];
    float *boot_exp_corrs_ptr = new float[pairs_num];
    float *boot_stat_ptr = new float[pairs_num];
     
    int *boot_ref_ind_ptr = new int[ref_ind_size];
    int *boot_exp_ind_ptr = new int[exp_ind_size];

    // Bootstrap indexes initialization
    std::vector<int> indexes(sample_size);
    for (int i = 0; i < sample_size; ++i) {
        indexes[i] = i;
    }

    // Random generator initialization
    // std::random_device random_dev;
    std::mt19937 random_gen(SEED);

    // Bootstrap pvalue computations    
    float *dpr, *rcp, *ecp, *sp, *pp;
    int *rip, *eip;
    
    std::cout << "Permutation progress: ";
    progressbar bar(repeats_number + 1);
    for (int r = 0; r < repeats_number + 1; ++r) {
        if (PyErr_CheckSignals() != 0) {
            throw py::error_already_set();
        }

        if (r == 0) {
            rcp = ref_corrs_ptr;
            ecp = exp_corrs_ptr;
            sp  = stat_ptr;
            pp  = pvalue_ptr;

            rip = ref_ind_ptr;
            eip = exp_ind_ptr;
        } else {
            std::shuffle(indexes.begin(), indexes.end(), random_gen);
            for (int i = 0; i < ref_ind_size; ++i) {
                boot_ref_ind_ptr[i] = indexes[i];
            }
            for (int i = 0; i < exp_ind_size; ++i) {
                boot_exp_ind_ptr[i] = indexes[ref_ind_size + i];
            }

            rcp = boot_ref_corrs_ptr;
            ecp = boot_exp_corrs_ptr;
            sp  = boot_stat_ptr;
            pp  = nullptr;
            
            rip = boot_ref_ind_ptr;
            eip = boot_exp_ind_ptr;
        }
        
        if (correlation == SPEARMAN) {
            rank_data(
                data_ptr,
                rank_ptr,
                sample_size,
                index_size,
                rip,
                ref_ind_size,
                process_num
            );

            rank_data(
                data_ptr,
                rank_ptr,
                sample_size,
                index_size,
                eip,
                exp_ind_size,
                process_num
            );

            dpr = rank_ptr;
        } else {
            dpr = data_ptr;
        }

        std::queue<std::thread> threads;
        int batch_size = pairs_num / process_num;
        for (int i = 0; i < process_num; ++i) {
            int left_border = i * batch_size;
            int right_border = (i + 1) * batch_size;
            if (i == process_num - 1) {
                right_border = pairs_num;
            }
            
            std::thread thr(ztest_pipeline,
                dpr,
                sample_size,
                (int *) nullptr,
                (int *) nullptr,
                left_border,
                right_border,
                index_size,
                rip,
                eip,
                ref_ind_size,
                exp_ind_size,
                rcp,
                ecp,
                sp,
                pp,
                correlation,
                alternative
            );
            
            threads.push(move(thr));
        }

        while (!threads.empty()) {
            threads.front().join();
            threads.pop();
        }
        
        if (r > 0) {
            for (int i = 0; i < pairs_num; ++i) {
                if ((alternative == TWO_SIDED) &&
                        (std::abs(stat_ptr[i]) <= std::abs(boot_stat_ptr[i]))) {
                    boot_pvalue_ptr[i] += 1;
                }

                if ((alternative == LESS) &&
                        (stat_ptr[i] <= boot_stat_ptr[i])) {
                    boot_pvalue_ptr[i] += 1;
                }
                
                if ((alternative == GREATER) &&
                        (stat_ptr[i] >= boot_stat_ptr[i])) {
                    boot_pvalue_ptr[i] += 1;
                }
            }
        } else {
            for (int i = 0; i < pairs_num; ++i) {
                boot_pvalue_ptr[i] = 0;
            }
        }

        bar.update();
    }
    std::cout << "\n";
    
    if (repeats_number > 0) {
        for (int i = 0; i < pairs_num; ++i) {
            boot_pvalue_ptr[i] /= repeats_number;
        }
    }
    
    if (correlation == SPEARMAN) {
        delete [] rank_ptr;
    }

    delete[] boot_ref_corrs_ptr;
    delete[] boot_exp_corrs_ptr;
    delete[] boot_stat_ptr;
     
    delete[] boot_ref_ind_ptr;
    delete[] boot_exp_ind_ptr;
    
    return std::tuple<
        NumPyFloatArray,
        NumPyFloatArray,
        NumPyFloatArray,
        NumPyFloatArray,
        NumPyFloatArray
    >(
        ref_corrs,
        exp_corrs,
        stat,
        pvalue,
        boot_pvalue
    );
}

std::tuple<
    NumPyFloatArray, NumPyFloatArray
> _score_pipeline_exhaustive(
    const NumPyFloatArray &data,
    const NumPyIntArray &reference_indexes,
    const NumPyIntArray &experimental_indexes,
    const std::string correlation=SPEARMAN,
    const std::string score=MEAN,
    const std::string alternative=TWO_SIDED,
    int repeats_number=REPEATS_NUMBER,
    int process_num=1    
) {
    py::buffer_info data_buf = data.request();
    float *data_ptr = (float *) data_buf.ptr;
    int data_len    = data_buf.shape[0];
    int sample_size = data_buf.shape[1];
    int pairs_num   = data_len * (data_len - 1) / 2;

    py::buffer_info ref_ind_buf = reference_indexes.request();
    int ref_ind_size = ref_ind_buf.shape[0];
    int *ref_ind_ptr = (int *) ref_ind_buf.ptr;
    
    py::buffer_info exp_ind_buf = experimental_indexes.request();
    int exp_ind_size = exp_ind_buf.shape[0];
    int *exp_ind_ptr = (int *) exp_ind_buf.ptr;

    // Rank data    
    float *rank_ptr;
    if (correlation == SPEARMAN) {
        rank_ptr = new float[
            data_len * sample_size
        ];
    }
    
    // Bootstrapped data
    float *boot_ref_corrs_ptr =  new float[pairs_num];
    float *boot_exp_corrs_ptr = new float[pairs_num];
    float *boot_stat_ptr = new float[pairs_num];
     
    int *boot_ref_ind_ptr = new int[ref_ind_size];
    int *boot_exp_ind_ptr = new int[exp_ind_size];

    // Bootstrap indexes initialization
    std::vector<int> indexes(sample_size);
    for (int i = 0; i < sample_size; ++i) {
        indexes[i] = i;
    }

    // Random generator initialization
    // std::random_device random_dev;
    std::mt19937 random_gen(SEED);

    // Scores data
    int sources_size = data_len; 
    NumPyFloatArray scores = NumPyFloatArray(sources_size);
    float *scores_ptr = (float *) scores.request().ptr;
    
    NumPyFloatArray pvalues = NumPyFloatArray(sources_size);
    float *pvalues_ptr = (float *) pvalues.request().ptr;
    
    // Bootstrap scores
    float *boot_scores_ptr = new float[sources_size];
    
    if (process_num > sources_size) {
        process_num = sources_size;
    }

    if (process_num <= 0) {    
        throw std::runtime_error("Process number error");
    }

    // Bootstrap pvalue computations    
    float *dpr, *rcp, *ecp, *sp, *pp;
    int *rip, *eip;
    float *scp;
    
    std::cout << "Permutation progress: ";
    progressbar bar(repeats_number + 1);
    for (int r = 0; r < repeats_number + 1; ++r) {
        if (PyErr_CheckSignals() != 0) {
            throw py::error_already_set();
        }

        if (r == 0) {
            rcp = boot_ref_corrs_ptr;
            ecp = boot_exp_corrs_ptr;
            sp  = boot_stat_ptr;
            pp  = nullptr;

            rip = ref_ind_ptr;
            eip = exp_ind_ptr;

            scp = scores_ptr;
        } else {
            std::shuffle(indexes.begin(), indexes.end(), random_gen);
            for (int i = 0; i < ref_ind_size; ++i) {
                boot_ref_ind_ptr[i] = indexes[i];
            }
            for (int i = 0; i < exp_ind_size; ++i) {
                boot_exp_ind_ptr[i] = indexes[ref_ind_size + i];
            }

            rcp = boot_ref_corrs_ptr;
            ecp = boot_exp_corrs_ptr;
            sp  = boot_stat_ptr;
            pp  = nullptr;
            
            rip = boot_ref_ind_ptr;
            eip = boot_exp_ind_ptr;

            scp = boot_scores_ptr;
        }

        if (correlation == SPEARMAN) {
            rank_data(
                data_ptr,
                rank_ptr,
                sample_size,
                data_len,
                rip,
                ref_ind_size,
                process_num
            );

            rank_data(
                data_ptr,
                rank_ptr,
                sample_size,
                data_len,
                eip,
                exp_ind_size,
                process_num
            );

            dpr = rank_ptr;
        } else {
            dpr = data_ptr;
        }
        
        std::queue<std::thread> threads;
        int batch_size = pairs_num / process_num;
        for (int i = 0; i < process_num; ++i) {
            int left_border = i * batch_size;
            int right_border = (i + 1) * batch_size;
            if (i == process_num - 1) {
                right_border = pairs_num;
            }
            
            std::thread thr(ztest_pipeline,
                dpr,
                sample_size,
                (int *) nullptr,
                (int *) nullptr,
                left_border,
                right_border,
                sources_size,
                rip,
                eip,
                ref_ind_size,
                exp_ind_size,
                rcp,
                ecp,
                sp,
                pp,
                correlation,
                TWO_SIDED
            );
            
            threads.push(move(thr));
        }

        while (!threads.empty()) {
            threads.front().join();
            threads.pop();
        }

        batch_size = sources_size / process_num;
        for (int i = 0; i < process_num; ++i) {
            int left_border = i * batch_size;
            int right_border = (i + 1) * batch_size;
            if (i == process_num - 1) {
                right_border = sources_size;
            }
            
            std::thread thr(score_pipeline_exhaustive,
                sp,
                sources_size,
                left_border,
                right_border,
                scp,
                score,
                alternative
            );
            
            threads.push(move(thr));
        }

        while (!threads.empty()) {
            threads.front().join();
            threads.pop();
        }

        if (r > 0) {
            for (int i = 0; i < sources_size; ++i) {
                if (alternative == TWO_SIDED && 
                        std::abs(scores_ptr[i]) <= std::abs(boot_scores_ptr[i])) {
                    pvalues_ptr[i] += 1;
                }

                if (alternative == LESS &&
                        scores_ptr[i] >= boot_scores_ptr[i]) {
                    pvalues_ptr[i] += 1;
                }


                if (alternative == GREATER &&
                        scores_ptr[i] <= boot_scores_ptr[i]) { 
                    pvalues_ptr[i] += 1;
                }
            }
        } else {
            for (int i = 0; i < sources_size; ++i) {
                pvalues_ptr[i] = 0;
            }
        }

        bar.update();
    }
    std::cout << "\n";

    if (repeats_number > 0) {
        for (int i = 0; i < sources_size; ++i) {
            pvalues_ptr[i] /= repeats_number;
        }
    }
    
    if (correlation == SPEARMAN) {
        delete[] rank_ptr;
    }

    delete[] boot_ref_corrs_ptr;
    delete[] boot_exp_corrs_ptr;
    delete[] boot_stat_ptr;
     
    delete[] boot_ref_ind_ptr;
    delete[] boot_exp_ind_ptr;
    
    delete[] boot_scores_ptr;

    return std::tuple<
        NumPyFloatArray, NumPyFloatArray
    >(
        scores,
        pvalues
    );
}

PYBIND11_MODULE(pipelines, m) {
    m.def("_ztest_pipeline_indexed", &ztest_pipeline_indexed);
    m.def("_score_pipeline_indexed", &_score_pipeline_indexed);
    m.def("_ztest_pipeline_exhaustive", &ztest_pipeline_exhaustive);
    m.def("_score_pipeline_exhaustive", &_score_pipeline_exhaustive);
}
