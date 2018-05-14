add_custom_target(gpulbm COMMENT "Builds GPU targets.")
function(gpu_target_helper TARGET_STR NPROCS NTHREADS X_SIZE Y_SIZE Z_SIZE)
  add_executable(${TARGET_STR} main.cu)
  target_link_libraries(${TARGET_STR} PRIVATE metalbm)
  target_compile_options(${TARGET_STR} PRIVATE --expt-extended-lambda)
  target_compile_definitions(${TARGET_STR} PRIVATE
    NPROCS=${NPROCS}
    NTHREADS=${NTHREADS}
    GLOBAL_LENGTH_X=${X_SIZE}
    GLOBAL_LENGTH_Y=${Y_SIZE}
    GLOBAL_LENGTH_Z=${Z_SIZE})
  add_dependencies(gpulbm ${TARGET_STR})
endfunction()

function(gpu_target NPROCS NTHREADS X_SIZE Y_SIZE Z_SIZE)
  set(TARGET_STR gpulbm_${NPROCS}_${NTHREADS}_${X_SIZE}_${Y_SIZE}_${Z_SIZE})
  gpu_target_helper(${TARGET_STR} ${NPROCS} ${NTHREADS} ${X_SIZE} ${Y_SIZE} ${Z_SIZE})
endfunction()

function(gpu_target_nvtx NPROCS NTHREADS X_SIZE Y_SIZE Z_SIZE)
  set(TARGET_STR gpulbm_${NPROCS}_${NTHREADS}_${X_SIZE}_${Y_SIZE}_${Z_SIZE}_nvtx)
  gpu_target_helper(${TARGET_STR} ${NPROCS} ${NTHREADS} ${X_SIZE} ${Y_SIZE} ${Z_SIZE})
  target_compile_options(${TARGET_STR} PRIVATE -lnvToolsExt)
  target_compile_definitions(${TARGET_STR} PRIVATE USE_NVTX)
  #set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lnvToolsExt" )
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



foreach(params_list ${PARAMS})
  separate_arguments(params_list) # convert " " to ";" (i.e. lists)
  list(GET params_list 0 NPROCS)
  list(GET params_list 1 NTHREADS)
  list(GET params_list 2 X_SIZE)
  list(GET params_list 3 Y_SIZE)
  list(GET params_list 4 Z_SIZE)
  if(USE_CUDA)
    gpu_target(${NPROCS} ${NTHREADS} ${X_SIZE} ${Y_SIZE} ${Z_SIZE})
    if(USE_NVTX)
      gpu_target_nvtx(${NPROCS} ${NTHREADS} ${X_SIZE} ${Y_SIZE} ${Z_SIZE})
    endif()
  endif()
  cpu_target(${NPROCS} ${NTHREADS} ${X_SIZE} ${Y_SIZE} ${Z_SIZE})
endforeach()
