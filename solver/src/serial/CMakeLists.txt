file(GLOB_RECURSE SERIAL_sources ${PROJECT_SOURCE_DIR}/src/serial/*.cpp)
add_executable(lbm ${sources} ${SERIAL_sources})

include_directories("${PROJECT_SOURCE_DIR}/src/serial")

target_compile_features(lbm PRIVATE cxx_range_for)

target_link_libraries(lbm ${Boost_LIBRARIES})
