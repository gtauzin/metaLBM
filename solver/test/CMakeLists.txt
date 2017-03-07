add_executable(all_tests ${test_sources})

target_link_libraries(all_tests ${MPI_CXX_LIBRARIES}
                                ${Boost_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(all_tests PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(all_tests PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

add_test(NAME all_tests
         WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
         ${MPIEXEC_PREFLAGS} all_tests ${MPIEXEC_POSTFLAGS})
