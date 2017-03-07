# Project Properties
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

# Compilers
option(ICC "Use Intel compiler" OFF)

project(LBM)
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
add_definitions(-DNTHREADS=${NTHREADS})
add_definitions(-DNPROCS=${NPROCS})
add_definitions(${DATA_STRUCT})
add_definitions(${LOG})
message(STATUS "${NPROCS}")

if(SERIAL)
  add_definitions(-D_SERIAL)
endif()

if(OMP)
  add_definitions(-D_OMP)
endif()

if(MPI)
  add_definitions(-D_MPI)
endif()

if(MPI_OMP)
  add_definitions(-D_MPI_OMP)
endif()

if(CUDA)
  add_definitions(-D_CUDA)
endif()

if(MPI_CUDA)
  add_definitions(-D_MPI_CUDA)
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

# Hostname
execute_process(COMMAND hostname OUTPUT_VARIABLE HOST_NAME OUTPUT_STRIP_TRAILING_WHITESPACE)

# Libraries
## Boost
if( ${HOST_NAME} MATCHES "hpc" )
  set(BOOST_ROOT "/home/rookstar/software/")
  set(BOOST_INCLUDEDIR "/home/rookstar/software/include")
  set(BOOST_LIBRARYDIR "/home/rookstar/software/lib")
  add_definitions(-DBOOST_LOG_DYN_LINK)
  set(Boost_USE_STATIC_LIBS OFF) # enable dynamic linking
elseif( ${HOST_NAME} MATCHES "frontend" )
  set(BOOST_ROOT "/home/tauzin/software/boost/")
  set(BOOST_INCLUDEDIR "/home/tauzin/software/boost/include")
  set(BOOST_LIBRARYDIR "/home/tauzin/software/boost/lib")
  add_definitions(-DBOOST_LOG_DYN_LINK)
  set(Boost_USE_STATIC_LIBS OFF) # enable dynamic linking
elseif( ${HOST_NAME} MATCHES "newturb" )
  set(BOOST_ROOT "/scratch/tauzin/software/")
  set(BOOST_INCLUDEDIR "/scratch/tauzin/software/include")
  set(BOOST_LIBRARYDIR "/scratch/tauzin/software/lib")
  set(Boost_USE_STATIC_LIBS ON) # enable static linking
elseif( ${HOST_NAME} MATCHES "node*" )
#  set(BOOST_ROOT "/galileo/home/userexternal/gtauzin0/software/")
#  set(BOOST_INCLUDEDIR "/galileo/home/userexternal/gtauzin0/software/include")
#  set(BOOST_LIBRARYDIR "/galileo/home/userexternal/gtauzin0/software/lib")
  set(Boost_USE_STATIC_LIBS ON) # enable static linking
elseif( ${HOST_NAME} MATCHES "hydra-08" )
  set(BOOST_ROOT "/home/tauzin/software/")
  set(BOOST_INCLUDEDIR "/home/tauzin/software/include")
  set(BOOST_LIBRARYDIR "/home/tauzin/software/lib")
  add_definitions(-DBOOST_LOG_DYN_LINK)
  set(Boost_USE_STATIC_LIBS OFF) # enable dynamic linking
else()
  add_definitions(-DBOOST_LOG_DYN_LINK)
  set(Boost_USE_STATIC_LIBS OFF) # enable dynamic linking
endif()

set(Boost_USE_MULTITHREAD ON) # enable multithreading
find_package(Boost REQUIRED)
find_package(Boost QUIET COMPONENTS log log_setup system chrono thread filesystem unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

## FFTW
if( ${HOST_NAME} MATCHES "hpc" )
  set(FFTW_INCLUDE_DIR "/home/rookstar/software/fftw/include")
  set(FFTW_LIBRARY_DIR "/home/rookstar/software/fftw/lib")
elseif( ${HOST_NAME} MATCHES "frontend" )
  set(FFTW_INCLUDE_DIR "/home/tauzin/software/fftw/include")
  set(FFTW_LIBRARY_DIR "/home/tauzin/software/fftw/lib")
elseif( ${HOST_NAME} MATCHES "newturb" )
  set(FFTW_INCLUDE_DIR "/scratch/software/fftw/include")
  set(FFTW_LIBRARY_DIR "/scratch/software/fftw/lib")
elseif( ${HOST_NAME} MATCHES "node*" )
  set(FFTW_INCLUDE_DIR "/galileo/home/userexternal/gtauzin0/software/fftw/include")
  set(FFTW_LIBRARY_DIR "/galileo/home/userexternal/gtauzin0/software/fftw/lib")
elseif( ${HOST_NAME} MATCHES "hydra-08" )
else()
endif()

if (FFTW_INCLUDE_DIR)
  set (FFTW_FIND_QUIETLY FALSE)
endif (FFTW_INCLUDE_DIR)

find_path (FFTW_INCLUDES NAMES fftw3.h PATHS ${FFTW_INCLUDE_DIR})
find_library (FFTW_LIBRARIES NAMES fftw3 PATHS ${FFTW_LIBRARY_DIR})

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)

mark_as_advanced (FFTW_LIBRARIES FFTW_INCLUDES)
include_directories(${FFTW_INCLUDES})

## OPENMP
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
  #set(CUDA_SEPARABLE_COMPILATION ON)
  #set(CUDA_PROPAGATE_HOST_FLAGS OFF)
  list(APPEND CUDA_NVCC_FLAGS -std=c++11 -gencode arch=compute_20,code=sm_20) #-std=c++11 -O3 -Xcompiler "-funroll-loops"
endif()

# Build Targets
add_subdirectory(src)
add_subdirectory(test)

# CMake debug
# get_cmake_property(_variableNames VARIABLES)
# foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()
