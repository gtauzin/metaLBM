add_custom_target(examples COMMENT "Builds all examples.")

function(metalbm_add_exec targetName targetFile)
    add_executable(${targetName} ${targetFile})
    target_link_libraries(${targetName}
        PRIVATE metaLBM::metaLBM)
    add_dependencies(examples ${targetName})
endfunction()

metalbm_add_exec(lbm_mpi main.cpp)
