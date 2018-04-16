#include_directories("${PROJECT_SOURCE_DIR}/example")
#add_executable(example_constexpr GlobalConstexpr.cu)
#target_link_libraries(example_constexpr PRIVATE metalbm)
#set_target_properties(example_constexpr PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#
#add_executable(templated_kernel TemplatedKernel.cu)
#target_link_libraries(templated_kernel PRIVATE metalbm)
#set_target_properties(templated_kernel PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
