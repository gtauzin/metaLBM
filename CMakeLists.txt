# Project Properties
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(metaLBM VERSION 0.0.1 LANGUAGES C CXX)

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

find_package(HDF5 REQUIRED)
target_include_directories(metalbm INTERFACE ${HDF5_INCLUDE_DIRS})
target_link_libraries(metalbm INTERFACE ${HDF5_LIBRARIES} -ldl -lz)

# Check for CUDA support
option(USE_CUDA "Enabling CUDA" ON)
if(USE_CUDA)
  enable_language(CUDA)
endif()

option(USE_FFTW "Enabling FFTW" ON)
if(USE_FFTW)
  find_package(FFTW REQUIRED)
  target_link_libraries(metalbm INTERFACE FFTW::FFTW)
  target_compile_definitions(metalbm INTERFACE USE_FFTW)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lm -pthread" )
endif()

option(USE_NVSHMEM "Enabling NVSHMEM" OFF)
if(USE_NVSHMEM)
  find_package(NVSHMEM REQUIRED)
  target_link_libraries(metalbm INTERFACE NVSHMEM::NVSHMEM)
  target_compile_definitions(metalbm INTERFACE USE_NVSHMEM)
endif()

option(USE_NVTX "Enabling manual instrumentation with NVTX" OFF)

option(USE_SCOREP "Enabling manual instrumentation with Score-P" OFF)
if(USE_SCOREP)
  find_package(Scorep REQUIRED)
  include_directories(${SCOREP_INCLUDE_DIRS})
  target_compile_definitions(metalbm INTERFACE USE_SCOREP)
endif()


# Build Targets
enable_testing()
add_subdirectory(src)
add_subdirectory(example)
#add_subdirectory(test)
#add_subdirectory(doc)
