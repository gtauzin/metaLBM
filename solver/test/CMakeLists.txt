add_custom_target(tests COMMENT "Builds all unit tests.")

# Custom test functions

function(metaLBM_add_test targetName targetFile)
    add_executable(${targetName} ${targetFile})
    target_link_libraries(${targetName} ${Boost_LIBRARIES})
    add_test(NAME ${targetName}
            WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
            COMMAND ${targetName})
    add_dependencies(tests ${targetName})
endfunction()

function(metaLBM_add_cuda_test targetName targetFile)
    cuda_add_executable(${targetName} ${targetFile})
    target_link_libraries(${targetName} ${Boost_LIBRARIES})
    add_test(NAME ${targetName}
            WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
            COMMAND ${targetName})
    add_dependencies(tests ${targetName})
endfunction()

function(metaLBM_add_mpi_test targetName targetFile numProcs)
    add_executable(${targetName} ${targetFile})
    target_link_libraries(${targetName} ${MPI_CXX_LIBRARIES}
      ${Boost_LIBRARIES})
    set_target_properties(${targetName} PROPERTIES
            COMPILE_FLAGS "${MPI_CXX_COMPILER_FLAGS}"
            LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    add_test(NAME ${targetName}
            WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
            COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${numProcs}
            ${MPIEXEC_PREFLAGS} ${targetName} ${MPIEXEC_POSTFLAGS})
    add_dependencies(tests ${targetName})
endfunction()

metaLBM_add_test(test_commons test_commons.cpp)
metaLBM_add_test(test_calculate test_calculate.cpp)
