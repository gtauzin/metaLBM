add_executable(lbm_mpi main.cpp)
target_link_libraries(lbm_mpi PRIVATE metalbm)
target_include_directories(lbm_mpi PRIVATE "./")
