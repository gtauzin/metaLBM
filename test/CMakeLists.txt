## Boost
find_package(Boost REQUIRED COMPONENTS system filesystem unit_test_framework)
if(NOT Boost_FOUND)
    message(WARNING
            "Boost was not found; tests will not be available!")
    return()
endif()

add_custom_target(tests
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
        COMMENT "Builds all unit tests.")

# Custom test functions
function(metaLBM_add_test targetName targetFile numProcs)
    add_executable(${targetName} ${targetFile})
    target_link_libraries(${targetName}
            metalbm
            Boost::system
            Boost::filesystem
            Boost::unit_test_framework)

    if (${numProcs} EQUAL 1)
        add_test(NAME ${targetName}
                WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
                COMMAND ${targetName})
    else ()
        add_test(NAME ${targetName}
                WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
                COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${numProcs}
                ${MPIEXEC_PREFLAGS} ./${targetName} ${MPIEXEC_POSTFLAGS})
    endif ()

    add_dependencies(tests ${targetName})
endfunction()


# Core headers tests
# metaLBM_add_test(TestInitialize TestInitialize.cpp 1)
# metaLBM_add_test(TestEquilibrium TestEquilibrium.cpp 1)
# metaLBM_add_test(TestCollision TestCollision.cpp 1)
metaLBM_add_test(TestForce TestForce.cpp 1)
# metaLBM_add_test(TestWriter TestWriter.cpp 1)
# metaLBM_add_test(TestDomain TestDomain.cpp 1)
# metaLBM_add_test(TestCommons TestCommons.cpp 1)
# metaLBM_add_test(TestMoment TestMoment.cpp 1)
# metaLBM_add_test(TestForcingScheme TestForcingScheme.cpp 1)
# metaLBM_add_test(TestField TestField.cpp 1)
# metaLBM_add_test(TestDistribution TestDistribution.cpp 1)
# metaLBM_add_test(TestBoundary TestBoundary.cpp 1)
# metaLBM_add_test(TestHelpers TestHelpers.cpp 1)
# metaLBM_add_test(TestReader TestReader.cpp 1)
# metaLBM_add_test(TestLattice TestLattice.cpp 1)
# metaLBM_add_test(TestMathVector TestMathVector.cpp 1)
# metaLBM_add_test(TestStaticArray TestStaticArray.cpp 1)
# metaLBM_add_test(TestDynamicArray TestDynamicArray.cpp 1)
