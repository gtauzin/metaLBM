file(GLOB_RECURSE MPI_CUDA_sources_cpp ${PROJECT_SOURCE_DIR}/src/mpi_cuda/*.cpp)
file(GLOB_RECURSE MPI_CUDA_sources_cu ${PROJECT_SOURCE_DIR}/src/mpi_cuda/*.cu)
file(GLOB_RECURSE MPI_CUDA_headers ${PROJECT_SOURCE_DIR}/src/mpi_cuda/*.h)

include_directories("${PROJECT_SOURCE_DIR}/src/mpi_cuda")

set_source_files_properties(${sources} ${MPI_CUDA_sources_cpp} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)

cuda_add_executable(lbm ${headers} ${MPI_CUDA_headers} ${sources} ${MPI_CUDA_sources_cpp} ${MPI_CUDA_sources_cu})


target_link_libraries(lbm ${MPI_CXX_LIBRARIES}
                          ${CUDA_LIBRARIES}
                          ${Boost_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(lbm PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(lbm PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
