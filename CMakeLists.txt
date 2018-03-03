# Project Properties
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(metaLBM VERSION 0.0.1 LANGUAGES C CXX)
#add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)

# Check for CUDA support
include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
  enable_language(CUDA)
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_70 --expt-extended-lambda")
  set(MPI_C_COMPILE_OPTIONS "")
  set(MPI_CXX_COMPILE_OPTIONS "")
  set(MPI_C_COMPILE_FLAGS "")
  set(MPI_CXX_COMPILE_FLAGS "")
  #  add_definitions(-O0 -lineinfo -g)
  add_definitions(-O0 -g)
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

find_package(MPI REQUIRED)
target_link_libraries(metalbm INTERFACE MPI::MPI_CXX)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP has been found but is not linked.")
  #target_link_libraries(metalbm INTERFACE OpenMP::OpenMP_CXX)
endif()

find_package(HDF5 REQUIRED)
if(HDF5_FOUND)
  target_include_directories(metalbm INTERFACE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(metalbm INTERFACE ${HDF5_LIBRARIES})
  # target_link_libraries(metalbm INTERFACE HDF5::HDF5)
  add_definitions(-DUSE_HDF5 -lz -ldl)
endif()

find_package(RapidXML REQUIRED)
target_link_libraries(metalbm INTERFACE RapidXML::RapidXML)

option(USE_FFTW "Enabling FFTW" OFF)
if(USE_FFTW)
  find_package(FFTW)
  if(FFTW_FOUND)
    target_link_libraries(metalbm INTERFACE FFTW::FFTW)
    add_definitions(-DUSE_FFTW)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lm -pthread" )
    else()
      message(STATUS "FFTW was not found!")
    endif()
endif()

option(USE_NVSHMEM "Enabling NVSHMEM" OFF)
if(USE_NVSHMEM)
  find_package(NVSHMEM)
  if(NVSHMEM_FOUND)
    target_link_libraries(metalbm INTERFACE NVSHMEM::NVSHMEM)
    add_definitions(-DUSE_NVSHMEM)
    else()
      message(STATUS "NVSHMEM was not found!")
    endif()
endif()

option(USE_SCOREP "Enabling manual instrumentation with Score-P" OFF)
if(USE_SCOREP)
  find_package(Scorep)
  if(SCOREP_FOUND)
        include_directories(${SCOREP_INCLUDE_DIRS})
        add_definitions(-DUSE_SCOREP)
    else()
        message(STATUS "Score-P was not found!")
    endif()
endif()

# if CUDA present?
option(USE_NVTX "Enabling manual instrumentation with NVTX" OFF)
if(USE_NVTX)
  set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lnvToolsExt" )
  add_definitions(-DUSE_NVTX)
endif()


# Build Targets
enable_testing()
add_subdirectory(src)
add_subdirectory(example)
#add_subdirectory(test)
#add_subdirectory(doc)
