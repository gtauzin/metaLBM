# Project Properties
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

# Compilers
option(ICC "Use Intel compiler" OFF)

project(metaLBM)
set (CMAKE_CXX_STANDARD 11)

# Include directories
include_directories("${PROJECT_SOURCE_DIR}/include/rapidxml-1.13")
include_directories("${PROJECT_SOURCE_DIR}/src/core")
file(GLOB_RECURSE sources src/core/*.cpp)
file(GLOB_RECURSE headers ${PROJECT_SOURCE_DIR}/src/core/*.h)

file(GLOB_RECURSE test_sources test/*.cpp)
set(test_sources ${test_sources} ${sources})

# Project Output
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin")

# Options
option(SERIAL "Enabling serial implementation" OFF)
option(OMP "Enabling OpenMP implementation" OFF)
option(MPI "Enabling parallel MPI implementation" OFF)
option(MPI_OMP "Enabling parallel MPI+OMP implementation" OFF)
option(CUDA "Enabling parallel CUDA implementation" OFF)
option(MPI_CUDA "Enabling parallel MPI+CUDA implementation" OFF)

link_libraries(-lpthread)

# Definitions
add_definitions(-D_SOA)
add_definitions(${DATA_STRUCT})
add_definitions(-DDATA_TYPE=double)
add_definitions(${DATA_TYPE})
add_definitions(-DLOG=)
add_definitions(${LOG})
message(STATUS "${NPROCS}")

if(SERIAL)
  add_definitions(-D_SERIAL)
  add_definitions(-DNPROCS=1)
  add_definitions(-DNTHREADS=1)
endif()

if(OMP)
  add_definitions(-D_OMP)
  add_definitions(-DNPROCS=1)
  add_definitions(-DNTHREADS=${NTHREADS})
endif()

if(MPI)
  add_definitions(-D_MPI)
  add_definitions(-DNPROCS=${NPROCS})
  add_definitions(-DNTHREADS=1)
endif()

if(MPI_OMP)
  add_definitions(-D_MPI_OMP)
  add_definitions(-DNPROCS=${NPROCS})
  add_definitions(-DNTHREADS=${NTHREADS})
endif()

if(CUDA)
  add_definitions(-D_CUDA)
  add_definitions(-DNPROCS=1)
  add_definitions(-DNTHREADS=1)
endif()

if(MPI_CUDA)
  add_definitions(-D_MPI_CUDA)
  add_definitions(-DNPROCS=${NPROCS})
  add_definitions(-DNTHREADS=1)
endif()

# Flags
if(ICC)
  set(ENV{CC} "icc")
  set(ENV{CXX} "icpc")
  set(ICC_COVERAGE_COMPILE_FLAGS "-std=c++11 -O3 -Wall -Wno-unknown-pragmas -pedantic -Wfatal-errors")
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${ICC_COVERAGE_COMPILE_FLAGS}")
  set(ICC_COVERAGE_LINK_FLAGS    "")
  set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${ICC_COVERAGE_LINK_FLAGS}")
else()
  set(GCC_COVERAGE_COMPILE_FLAGS "-std=c++11 -O3 -Wall -Wno-unknown-pragmas -pedantic")# -Wfatal-errors")
  #set(GCC_COVERAGE_COMPILE_FLAGS "-std=c++11 -O0 -g -Wall")
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
  set(GCC_COVERAGE_LINK_FLAGS    " ")
  set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -pg")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g -pg")
endif()

# Libraries
## Boost
set(BOOST_ROOT $ENV{BOOST_ROOT})
set(BOOST_INCLUDEDIR $ENV{BOOST_INCLUDEDIR})
set(BOOST_LIBRARYDIR $ENV{BOOST_LIBRARYDIR})

add_definitions(-DBOOST_LOG_DYN_LINK)
set(Boost_USE_STATIC_LIBS OFF) # enable static linking
set(Boost_USE_MULTITHREAD ON) # enable multithreading
find_package(Boost REQUIRED)
find_package(Boost QUIET COMPONENTS log log_setup system chrono thread filesystem unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

## FFTW
set(FFTW_INCLUDE_DIR $ENV{FFTW_INCLUDE_DIR})
set(FFTW_LIBRARY_DIR $ENV{FFTW_LIBRARY_DIR})

if (FFTW_INCLUDE_DIR)
  set (FFTW_FIND_QUIETLY FALSE)
endif (FFTW_INCLUDE_DIR)

find_path (FFTW_INCLUDES NAMES fftw3.h PATHS ${FFTW_INCLUDE_DIR})
find_library (FFTW_LIBRARIES NAMES fftw3 PATHS ${FFTW_LIBRARY_DIR})

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)

mark_as_advanced (FFTW_LIBRARIES FFTW_INCLUDES)
include_directories(${FFTW_INCLUDES})

## OpenMP
if (OMP OR MPI_OMP)
  find_package(OpenMP)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

## MPI
if(MPI OR MPI_CUDA OR MPI_OMP)
  find_package(MPI QUIET REQUIRED)
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
endif()

## CUDA
if(CUDA OR MPI_CUDA)
  message(STATUS "CUDA OR MPI_CUDA ON")
  find_package(CUDA QUIET REQUIRED)
  set(GENCODE_SM13 -gencode arch=compute_13,code=sm_13)
  set(GENCODE_SM20 -gencode arch=compute_20,code=sm_20)
  set(GENCODE_SM30 -gencode arch=compute_30,code=sm_30)
  set(GENCODE_SM35 -gencode arch=compute_35,code=sm_35)
  set(GENCODE_SM35 -gencode arch=compute_37,code=sm_37)
  set(GENCODE_SM50 -gencode arch=compute_50,code=sm_50)
  set(GENCODE_SM52 -gencode arch=compute_52,code=\"sm_52,compute_52\")
  set(GENCODE_SM60 -gencode arch=compute_60,code=\"sm_60,compute_60\")
  set(GENCODE_FLAGS ${GENCODE_SM20} ${GENCODE_SM30} ${GENCODE_SM35}
    ${GENCODE_SM37} ${GENCODE_SM50} ${GENCODE_SM52}
    ${GENCODE_SM60})
  set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};
    ${GENCODE_FLAGS};
    -O3;
    -std=c++11)
  include_directories("${CUDA_INCLUDE_DIRS}")
endif()

# Build Targets
enable_testing()
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(doc)

# CMake debug
# get_cmake_property(_variableNames VARIABLES)
# foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()
