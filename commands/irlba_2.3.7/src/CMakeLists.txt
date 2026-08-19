cmake_minimum_required(VERSION 3.10)
project(irlba C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Enable Position Independent Code (required for shared libraries)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Find BLAS and LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Define includes and library paths for R and Matrix package
set(R_INCLUDE_DIR "/usr/local/lib/R/include")
set(R_LIB_DIR "/usr/local/lib/R/lib")
set(Matrix_INCLUDE_DIR "/usr/local/lib/R/site-library/Matrix/include")

# Find R library
find_library(R_LIBRARY NAMES R PATHS ${R_LIB_DIR} NO_DEFAULT_PATH)
if(NOT R_LIBRARY)
    message(FATAL_ERROR "R library not found in ${R_LIB_DIR}")
endif()

# Source files
set(SOURCES
    irlb.c
    utility.c
)

# Shared library target
add_library(irlba SHARED ${SOURCES})

# Target include directories
target_include_directories(irlba PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${R_INCLUDE_DIR}
    ${Matrix_INCLUDE_DIR}
)

# Target link libraries
target_link_libraries(irlba PRIVATE
    ${R_LIBRARY}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
    m
)


