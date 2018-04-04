add_custom_target(gpulbm COMMENT "Builds GPU targets.")
function(gpu_target NPROCS NTHREADS X_SIZE Y_SIZE Z_SIZE)
  set(target_name gpulbm_${NPROCS}_${NTHREADS}_${X_SIZE}_${Y_SIZE}_${Z_SIZE})
  add_executable(${target_name} main.cu)
  target_link_libraries(${target_name} PRIVATE metalbm)
  target_compile_options(${target_name} PRIVATE --expt-extended-lambda)
  target_compile_definitions(${target_name} PRIVATE
    NPROCS=${NPROCS}
    NTHREADS=${NTHREADS}
    GLOBAL_LENGTH_X=${X_SIZE}
    GLOBAL_LENGTH_Y=${Y_SIZE}
    GLOBAL_LENGTH_Z=${Z_SIZE})
  add_dependencies(gpulbm ${target_name})
endfunction()

add_custom_target(cpulbm COMMENT "Builds CPU targets.")
function(cpu_target NPROCS NTHREADS X_SIZE Y_SIZE Z_SIZE)
  set(target_name cpulbm_${NPROCS}_${NTHREADS}_${X_SIZE}_${Y_SIZE}_${Z_SIZE})
  add_executable(${target_name} main.cpp)
  target_link_libraries(${target_name} PRIVATE metalbm)
  target_compile_definitions(${target_name} PRIVATE
    NPROCS=${NPROCS}
    NTHREADS=${NTHREADS}
    GLOBAL_LENGTH_X=${X_SIZE}
    GLOBAL_LENGTH_Y=${Y_SIZE}
    GLOBAL_LENGTH_Z=${Z_SIZE})
  add_dependencies(cpulbm ${target_name})
endfunction()

message(STATUS "les params sont ${PARAMS}")
foreach(params_list ${PARAMS})
  message(STATUS "la liste est ${params_list}")
  separate_arguments(params_list) # convert " " to ";" (i.e. lists)
  list(GET params_list 0 NPROCS)
  list(GET params_list 1 NTHREADS)
  list(GET params_list 2 X_SIZE)
  list(GET params_list 3 Y_SIZE)
  list(GET params_list 4 Z_SIZE)
  if(USE_CUDA)
    gpu_target(${NPROCS} ${NTHREADS} ${X_SIZE} ${Y_SIZE} ${Z_SIZE})
  endif()
  cpu_target(${NPROCS} ${NTHREADS} ${X_SIZE} ${Y_SIZE} ${Z_SIZE})
endforeach()
