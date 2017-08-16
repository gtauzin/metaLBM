# Project Properties
cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(metaLBM VERSION 0.0.1 LANGUAGES CXX)
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)

# Check for CUDA support
include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
else()
    message(STATUS "No CUDA support")
endif()

# Setup the `metalbm` target
add_library(metalbm INTERFACE)
target_compile_features(metalbm INTERFACE cxx_std_11)
target_include_directories(metalbm INTERFACE
        "$<BUILD_INTERFACE:${metaLBM_SOURCE_DIR}/include>"
        $<INSTALL_INTERFACE:include>)

find_package(RapidXML REQUIRED)
target_link_libraries(metalbm INTERFACE RapidXML::RapidXML)

find_package(MPI)  # Imported library support in cmake 3.9; we need the module locally
if(MPI_CXX_FOUND)
    target_link_libraries(metalbm INTERFACE MPI::MPI_CXX)
endif()

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(metalbm INTERFACE OpenMP::OpenMP_CXX)
endif()

#find_package(FFTW) # TODO: FFTW no FOUND_LIBRARY
#if(FFTW_FOUND) 
    # TODO
#endif()

# Build Targets
enable_testing()
add_subdirectory(src/mpi)
add_subdirectory(src/serial)
add_subdirectory(test)
add_subdirectory(doc)
