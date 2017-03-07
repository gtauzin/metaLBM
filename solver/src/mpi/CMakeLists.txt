file(GLOB_RECURSE MPI_sources ${PROJECT_SOURCE_DIR}/src/mpi/*.cpp)
add_executable(lbm ${sources} ${MPI_sources})

include_directories("${PROJECT_SOURCE_DIR}/src/mpi")

target_link_libraries(lbm
  ${MPI_CXX_LIBRARIES}
  ${Boost_LIBRARIES}
  ${FFTW_LIBRARIES}
  )

if(MPI_OMP)
  target_link_libraries(lbm ${OpenMP_CXX_LIBRARIES})
endif()

target_compile_features(lbm PRIVATE cxx_range_for)

if(MPI_COMPILE_FLAGS)
  set_target_properties(lbm PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(lbm PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
