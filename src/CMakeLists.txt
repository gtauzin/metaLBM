include_directories("${PROJECT_SOURCE_DIR}/src")
add_executable(cpulbm main.cpp)
target_link_libraries(cpulbm PRIVATE metalbm)

add_executable(gpulbm main.cu)
target_link_libraries(gpulbm PRIVATE metalbm)
set_target_properties(gpulbm PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
