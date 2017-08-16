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

# Core headers tests
metaLBM_add_test(TestInitialize TestInitialize.cpp)
metaLBM_add_test(TestEquilibrium TestEquilibrium.cpp)
metaLBM_add_test(TestCollision TestCollision.cpp)
metaLBM_add_test(TestForce TestForce.cpp)
metaLBM_add_test(TestWriter TestWriter.cpp)
metaLBM_add_test(TestDomain TestDomain.cpp)
metaLBM_add_test(TestCommons TestCommons.cpp)
metaLBM_add_test(TestMoment TestMoment.cpp)
metaLBM_add_test(TestForcingScheme TestForcingScheme.cpp)
metaLBM_add_test(TestField TestField.cpp)
metaLBM_add_test(TestDistribution TestDistribution.cpp)
 metaLBM_add_test(TestBoundary TestBoundary.cpp)
metaLBM_add_test(TestHelpers TestHelpers.cpp)
metaLBM_add_test(TestReader TestReader.cpp)
metaLBM_add_test(TestLattice TestLattice.cpp)
metaLBM_add_test(TestMathVector TestMathVector.cpp)
metaLBM_add_test(TestStaticArray TestStaticArray.cpp)
metaLBM_add_test(TestDynamicArray TestDynamicArray.cpp)
