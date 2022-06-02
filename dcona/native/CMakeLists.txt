cmake_minimum_required(VERSION 3.4)
project(fast_computations)

# Debug regime
# set(CMAKE_BUILD_TYPE Debug)

# Standard regime
add_compile_options(-O3)
add_compile_options(-Wall)
add_compile_options(-march=native)

add_subdirectory(extern/pybind11)

pybind11_add_module(
 	utils SHARED
	src/utils/utils.cpp
	src/lib/utils.cpp
)

pybind11_add_module(
 	correlations SHARED
	src/utils/utils.cpp
	src/correlations/correlations.cpp
 	src/lib/correlations.cpp
)

pybind11_add_module(
	tests SHARED
	src/utils/utils.cpp
	src/correlations/correlations.cpp
	src/tests/tests.cpp
	src/lib/tests.cpp
)
 
pybind11_add_module(
 	scores SHARED
	src/utils/utils.cpp
	src/scores/scores.cpp
 	src/lib/scores.cpp
)

pybind11_add_module(
  	pipelines SHARED
	src/utils/utils.cpp
	src/correlations/correlations.cpp
	src/tests/tests.cpp
    src/scores/scores.cpp
	src/pipelines/pipelines.cpp
	src/lib/pipelines.cpp
)
