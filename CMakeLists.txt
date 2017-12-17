# Project Properties
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(metaLBM VERSION 0.0.1 LANGUAGES CXX)
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)

# Check for CUDA support
include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
  enable_language(CUDA)
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_35 --expt-extended-lambda")
  set(MPI_CXX_COMPILE_OPTIONS "")
  add_definitions(-O0 -lineinfo -G)
else()
    message(STATUS "No CUDA support")
endif()

# Project Output
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/bin")

# Setup the `metalbm` target
add_library(metalbm INTERFACE)
target_compile_features(metalbm INTERFACE cxx_std_11)
target_include_directories(metalbm INTERFACE
  "$<BUILD_INTERFACE:${metaLBM_SOURCE_DIR}/include>"
  $<INSTALL_INTERFACE:include>)

find_package(RapidXML REQUIRED)
target_link_libraries(metalbm INTERFACE RapidXML::RapidXML)

find_package(NVSHMEM)
if(NVSHMEM_FOUND)
  target_link_libraries(metalbm INTERFACE NVSHMEM::NVSHMEM)
endif()

#find_package(HDF5 REQUIRED)
#target_link_libraries(metalbm INTERFACE HDF5::HDF5)

find_package(MPI REQUIRED)
target_link_libraries(metalbm INTERFACE MPI::MPI_CXX)


#find_package(OpenMP)
#if (OpenMP_CXX_FOUND)
#    target_link_libraries(metalbm INTERFACE OpenMP::OpenMP_CXX)
#endif()

option(USE_SCOREP "Enabling manual instrumentation with Score-P" OFF)
if(USE_SCOREP)
  find_package(Scorep)
  if(SCOREP_FOUND)
        include_directories(${SCOREP_INCLUDE_DIRS})
        add_definitions(-DUSE_SCOREP)
    else()
        message("Score-P was not found!")
    endif()
endif()

# if CUDA presente?
option(USE_NVTX "Enabling manual instrumentation with NVTX" OFF)
if(USE_NVTX)
  add_definitions(-DUSE_NVTX)
endif()

#find_package(FFTW) # TODO:

# Build Targets
enable_testing()
add_subdirectory(src)
add_subdirectory(example)
#add_subdirectory(test)
#add_subdirectory(doc)
