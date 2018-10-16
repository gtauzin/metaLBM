find_package(Doxygen)
if (NOT DOXYGEN_FOUND)
    message(WARNING "Doxygen was not found; the 'doc' and 'docset' targets "
        "will be unavailable.")
    return()
endif()

configure_file(Doxyfile.in Doxyfile @ONLY)
add_custom_target(doc
    COMMAND ${CMAKE_COMMAND} -E remove_directory html
    COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)
