cmake_minimum_required(VERSION 3.10)
project(stats LANGUAGES C Fortran)

# Define R internal include paths based on your custom folder tree layout
set(R_INTERNAL_INCLUDE_DIR 
      "${CMAKE_CURRENT_SOURCE_DIR}/../../../include"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../../gnuwin32/fixed/h"
)
set(R_INCLUDE_DIR ${R_INTERNAL_INCLUDE_DIR})

# Dependencies
find_package(OpenMP)
find_package(BLAS)
find_package(LAPACK)

# Source file definitions
set(SOURCES_C
  init.c kmeans.c ansari.c bandwidths.c chisqsim.c d2x2xk.c fexact.c kendall.c 
  ks.c line.c smooth.c prho.c swilk.c ksmooth.c loessc.c monoSpl.c isoreg.c 
  Srunmed.c dblcen.c distance.c hclust-utils.c nls.c rWishart.c HoltWinters.c 
  PPsum.c arima.c burg.c filter.c mAR.c pacf.c starma.c port.c family.c 
  sbart.c approx.c loglin.c lowess.c massdist.c splines.c lm.c complete_cases.c 
  cov.c deriv.c fft.c fourier.c model.c optim.c optimize.c integrate.c 
  random.c distn.c zeroin.c rcont.c influence.c
)

set(SOURCES_F
  bsplvd.f bvalue.f bvalus.f loessf.f ppr.f qsbart.f sgram.f sinerp.f 
  sslvrg.f stxwx.f hclust.f kmns.f eureka.f stl.f portsrc.f lminfl.f
)

# Define target
add_library(stats SHARED ${SOURCES_C} ${SOURCES_F})

set_target_properties(stats PROPERTIES
  PREFIX ""
  SUFFIX ".so"
)

# Include directories alignment
target_include_directories(stats PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${R_INCLUDE_DIR}
  ${R_INTERNAL_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../../../main  
  ${CMAKE_CURRENT_SOURCE_DIR}/../../../appl  
)

# Force-inject Linux configurations since we skipped R's ./configure
target_compile_definitions(stats PRIVATE 
  HAVE_CONFIG_H
  SUPPORT_OPENMP
  HAVE_BLAS
)

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(stats PRIVATE -Wall -O2)
endif()

set_source_files_properties(portsrc.f PROPERTIES
  COMPILE_FLAGS "-O0"
)

# Linking libraries
target_link_libraries(stats PRIVATE m)

if(OpenMP_FOUND)
  target_link_libraries(stats PRIVATE OpenMP::OpenMP_C OpenMP::OpenMP_Fortran)
endif()

if(BLAS_FOUND)
  target_link_libraries(stats PRIVATE ${BLAS_LIBRARIES})
endif()

if(LAPACK_FOUND)
  target_link_libraries(stats PRIVATE ${LAPACK_LIBRARIES})
endif()