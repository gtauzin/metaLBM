add_custom_target(gpulbm COMMENT "Builds GPU targets.")
function(gpu_target_helper TARGET_STR NPROCS NTHREADS X_SIZE Y_SIZE Z_SIZE LBM_POSTFIX)
  add_executable(${TARGET_STR} main.cu)
  target_link_libraries(${TARGET_STR} PRIVATE metalbm)
  target_compile_options(${TARGET_STR} PRIVATE --expt-extended-lambda)
  target_compile_definitions(${TARGET_STR} PRIVATE
    NPROCS=${NPROCS}
    NTHREADS=${NTHREADS}
    GLOBAL_LENGTH_X=${X_SIZE}
    GLOBAL_LENGTH_Y=${Y_SIZE}
    GLOBAL_LENGTH_Z=${Z_SIZE}
    LBM_POSTFIX=${LBM_POSTFIX})
  set_target_properties(${TARGET_STR} PROPERTIES CUDA_SEPARABLE_COMPILATION OFF)
  add_dependencies(gpulbm ${TARGET_STR})
endfunction()

function(gpu_target NPROCS NTHREADS X_SIZE Y_SIZE Z_SIZE LBM_POSTFIX)
  set(TARGET_STR gpulbm_${NPROCS}_${NTHREADS}_${X_SIZE}_${Y_SIZE}_${Z_SIZE}_${LBM_POSTFIX})
  gpu_target_helper(${TARGET_STR} ${NPROCS} ${NTHREADS} ${X_SIZE} ${Y_SIZE} ${Z_SIZE} ${LBM_POSTFIX})
endfunction()

function(gpu_target_nvtx NPROCS NTHREADS X_SIZE Y_SIZE Z_SIZE LBM_POSTFIX)
  set(TARGET_STR gpulbm_${NPROCS}_${NTHREADS}_${X_SIZE}_${Y_SIZE}_${Z_SIZE}_${LBM_POSTFIX}_nvtx)
  gpu_target_helper(${TARGET_STR} ${NPROCS} ${NTHREADS} ${X_SIZE} ${Y_SIZE} ${Z_SIZE} ${LBM_POSTFIX})
  target_link_libraries(${TARGET_STR} PRIVATE -lnvToolsExt)
  target_compile_definitions(${TARGET_STR} PRIVATE USE_NVTX)
endfunction()

add_custom_target(cpulbm COMMENT "Builds CPU targets.")
function(cpu_target NPROCS NTHREADS X_SIZE Y_SIZE Z_SIZE LBM_POSTFIX)
  set(target_name cpulbm_${NPROCS}_${NTHREADS}_${X_SIZE}_${Y_SIZE}_${Z_SIZE}_${LBM_POSTFIX})
  add_executable(${target_name} main.cpp)
  target_link_libraries(${target_name} PRIVATE metalbm)
  target_compile_definitions(${target_name} PRIVATE
    NPROCS=${NPROCS}
    NTHREADS=${NTHREADS}
    GLOBAL_LENGTH_X=${X_SIZE}
    GLOBAL_LENGTH_Y=${Y_SIZE}
    GLOBAL_LENGTH_Z=${Z_SIZE}
    LBM_POSTFIX=${LBM_POSTFIX})
  add_dependencies(cpulbm ${target_name})
endfunction()


foreach(params_list ${PARAMS})
  separate_arguments(params_list) # convert " " to ";" (i.e. lists)
  list(GET params_list 0 NPROCS)
  list(GET params_list 1 NTHREADS)
  list(GET params_list 2 X_SIZE)
  list(GET params_list 3 Y_SIZE)
  list(GET params_list 4 Z_SIZE)
  list(GET params_list 5 LBM_POSTFIX)
  if(USE_CUDA)
    gpu_target(${NPROCS} ${NTHREADS} ${X_SIZE} ${Y_SIZE} ${Z_SIZE} ${LBM_POSTFIX})
    if(USE_NVTX)
      gpu_target_nvtx(${NPROCS} ${NTHREADS} ${X_SIZE} ${Y_SIZE} ${Z_SIZE} ${LBM_POSTFIX})
    endif()
  endif()
  cpu_target(${NPROCS} ${NTHREADS} ${X_SIZE} ${Y_SIZE} ${Z_SIZE} ${LBM_POSTFIX})
endforeach()
