add_executable(cpulbm main.cpp)
target_link_libraries(cpulbm PRIVATE metalbm)

if(USE_CUDA)
  add_executable(gpulbm main.cu)
  target_link_libraries(gpulbm PRIVATE metalbm)
  target_compile_options(gpulbm PRIVATE --expt-extended-lambda)
endif()

