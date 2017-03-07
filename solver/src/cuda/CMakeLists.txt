file(GLOB_RECURSE CUDA_sources ${PROJECT_SOURCE_DIR}/src/cuda/*.cpp)
add_executable(lbm ${sources} ${CUDA_sources})

include_directories("${PROJECT_SOURCE_DIR}/src/cuda")

target_compile_features(lbm PRIVATE cxx_range_for)

target_link_libraries(lbm ${Boost_LIBRARIES}
                          ${CUDA_LIBRARIES})
