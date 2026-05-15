cmake_minimum_required(VERSION 3.10)
project(SeuratObject)

# Set C++ standard to C++17 as per Makevars
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Configure Debug mode
set(CMAKE_BUILD_TYPE Debug)
add_compile_options(-g -O0)

# Discover R configuration
execute_process(
    COMMAND R RHOME
    OUTPUT_VARIABLE R_HOME
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
    COMMAND Rscript -e "cat(R.home('include'))"
    OUTPUT_VARIABLE R_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
    COMMAND Rscript -e "cat(system.file('include', package='Rcpp'))"
    OUTPUT_VARIABLE RCPP_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
    COMMAND Rscript -e "cat(system.file('include', package='RcppEigen'))"
    OUTPUT_VARIABLE RCPPEIGEN_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
    COMMAND Rscript -e "cat(system.file('include', package='RcppProgress'))"
    OUTPUT_VARIABLE RCPPPROGRESS_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Include directories
include_directories(
    ${R_INCLUDE_DIR}
    ${RCPP_INCLUDE_DIR}
    ${RCPPEIGEN_INCLUDE_DIR}
    ${RCPPPROGRESS_INCLUDE_DIR}
    .
)

# Source files
file(GLOB SOURCES "*.cpp" "*.c")

# Build the shared library
add_library(SeuratObject SHARED ${SOURCES})

# Set output name to SeuratObject.so (Linux convention for R)
set_target_properties(SeuratObject PROPERTIES 
    PREFIX ""
    OUTPUT_NAME "SeuratObject"
    SUFFIX ".so"
)

# Print configuration for verification
message(STATUS "R_HOME:                 ${R_HOME}")
message(STATUS "R_INCLUDE_DIR:          ${R_INCLUDE_DIR}")
message(STATUS "RCPP_INCLUDE_DIR:       ${RCPP_INCLUDE_DIR}")
message(STATUS "RCPPEIGEN_INCLUDE_DIR:  ${RCPPEIGEN_INCLUDE_DIR}")
message(STATUS "RCPPPROGRESS_INCLUDE_DIR: ${RCPPPROGRESS_INCLUDE_DIR}")
