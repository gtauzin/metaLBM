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
target_link_libraries(metalbm INTERFACE NVSHMEM::NVSHMEM)


#find_package(HDF5 REQUIRED)
#target_link_libraries(metalbm INTERFACE HDF5::HDF5)

find_package(MPI REQUIRED)
target_link_libraries(metalbm INTERFACE MPI::MPI_CXX)

find_package(OpenMP)
if (NOT OpenMP_CXX_FOUND)
    target_link_libraries(metalbm INTERFACE OpenMP::OpenMP_CXX)
endif()

option(PROFILE_SCOREP "Enabling manual instrumentation with Score-P" OFF)
if(PROFILE_SCOREP)
  find_package(Scorep)
  if(SCOREP_FOUND)
        include_directories(${SCOREP_INCLUDE_DIRS})
        add_definitions(-DPROFILE_SCOREP)
    else()
        message("Score-P was not found!")
    endif()
endif()

#find_package(FFTW) # TODO:

# Build Targets
enable_testing()
add_subdirectory(src)
add_subdirectory(example)
#add_subdirectory(test)
#add_subdirectory(doc)
