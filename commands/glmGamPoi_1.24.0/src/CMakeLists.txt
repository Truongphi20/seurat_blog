cmake_minimum_required(VERSION 3.10)
project(glmGamPoi)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Disable optimization and include debug symbols as requested ("without optimizers")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -fPIC")

# R and package include paths (absolute paths from the Docker image)
set(R_INCLUDE_DIR "/usr/local/lib/R/include")
set(RCPP_INCLUDE_DIR "/usr/local/lib/R/site-library/Rcpp/include")
set(RCPP_ARMADILLO_INCLUDE_DIR "/usr/local/lib/R/site-library/RcppArmadillo/include")
set(BEACHMAT_INCLUDE_DIR "/usr/local/lib/R/site-library/beachmat/include")
set(ASSORTHEAD_INCLUDE_DIR "/usr/local/lib/R/site-library/assorthead/include")
# The package's own headers are found in the installed location in the container
set(GLMGAMPOI_INCLUDE_DIR "/usr/local/lib/R/site-library/glmGamPoi/include")

include_directories(
    ${R_INCLUDE_DIR}
    ${RCPP_INCLUDE_DIR}
    ${RCPP_ARMADILLO_INCLUDE_DIR}
    ${BEACHMAT_INCLUDE_DIR}
    ${ASSORTHEAD_INCLUDE_DIR}
    ${GLMGAMPOI_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Source files
file(GLOB SOURCES "*.cpp")

add_library(glmGamPoi SHARED ${SOURCES})

# Libraries to link against (from R CMD config)
target_link_libraries(glmGamPoi
    blas
    lapack
    gfortran
    m
    quadmath
)

# Output name should be glmGamPoi.so
set_target_properties(glmGamPoi PROPERTIES 
    PREFIX ""
    SUFFIX ".so"
)

# Definitions common for R packages
target_compile_definitions(glmGamPoi PRIVATE R_NO_REMAP)
