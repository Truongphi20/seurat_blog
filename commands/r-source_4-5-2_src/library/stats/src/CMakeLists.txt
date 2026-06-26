cmake_minimum_required(VERSION 3.10)
project(stats LANGUAGES C Fortran)

# Define R internal include paths based on your custom folder tree layout
set(R_INTERNAL_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../include")

# Dependencies
find_package(OpenMP)
find_package(BLAS)
find_package(LAPACK)

# Retrieve additional BLAS/LAPACK flags from R configuration if available
if(R_COMMAND)
  execute_process(
    COMMAND ${R_COMMAND} CMD config BLAS_LIBS
    OUTPUT_VARIABLE R_BLAS_LIBS
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  execute_process(
    COMMAND ${R_COMMAND} CMD config LAPACK_LIBS
    OUTPUT_VARIABLE R_LAPACK_LIBS
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
endif()

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

# Properties match R plugin output naming conventions
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

# Compile configurations & macros
if(CONFIG_H_DIR)
  target_compile_definitions(stats PRIVATE HAVE_CONFIG_H)
  message(STATUS "HAVE_CONFIG_H enabled using: ${CONFIG_H_DIR}")
endif()

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(stats PRIVATE -Wall -O2)
endif()

# Specific low-optimization flags required by portsrc library routines
set_source_files_properties(portsrc.f PROPERTIES
  COMPILE_FLAGS "-O0"
)

# Linking libraries
target_link_libraries(stats PRIVATE m)

if(OpenMP_FOUND)
  target_link_libraries(stats PRIVATE OpenMP::OpenMP_C OpenMP::OpenMP_Fortran)
endif()

if(R_BLAS_LIBS)
  target_link_libraries(stats PRIVATE ${R_BLAS_LIBS})
elseif(BLAS_FOUND)
  target_link_libraries(stats PRIVATE ${BLAS_LIBRARIES})
endif()

if(R_LAPACK_LIBS)
  target_link_libraries(stats PRIVATE ${R_LAPACK_LIBS})
elseif(LAPACK_FOUND)
  target_link_libraries(stats PRIVATE ${LAPACK_LIBRARIES})
endif()