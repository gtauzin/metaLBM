# Copyright Teodor Nikolov 2017
# Distributed under the MIT License (MIT)

# .rst:
# FindRapidXML
# ------------
#
# Locate the XML parsing library RapidXml
#
#
# Imported targets
# ^^^^^^^^^^^^^^^^
# 
# ``RapidXML::RapidXML``
#
#
# Result variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project:
#
# ``RAPIDXML_FOUND``            - Found the RapidXml library
# ``RAPIDXML_INCLUDE_DIRS``     - RapidXML include directories
#
#
# Cache variables
# ^^^^^^^^^^^^^^^
# 
# The following Cache variables may also be set:
# 
# ``RAPIDXML_ROOT``             - The root directory of RapidXml installation
#                                 (may also be set as an environment variable)
#

set(RAPIDXML_INCLUDE_SEARCH_DIRS "")
if(RAPIDXML_ROOT)
    list(APPEND RAPIDXML_INCLUDE_SEARCH_DIRS 
                ${RAPIDXML_ROOT}/include 
                ${RAPIDXML_ROOT})
  elseif( $ENV{RAPIDXML_ROOT} )
      list(APPEND RAPIDXML_INCLUDE_SEARCH_DIRS 
                  $ENV{RAPIDXML_ROOT}/include 
                  $ENV{RAPIDXML_ROOT})
endif()

set(RAPIDXML_KNOWN_VERSIONS "1.0" "1.1" "1.11" "1.12" "1.13")

set(RAPIDXML_PATH_SUFFIXES) 
foreach(RAPIDXML_VERSION ${RAPIDXML_KNOWN_VERSIONS})
    list(APPEND RAPIDXML_PATH_SUFFIXES "rapidxml-${RAPIDXML_VERSION}")
endforeach()

find_path(
    RAPIDXML_INCLUDE_DIRS
    NAMES         rapidxml.hpp 
    HINTS         ${RAPIDXML_INCLUDE_SEARCH_DIRS}
    PATH_SUFFIXES ${RAPIDXML_PATH_SUFFIXES})

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(RapidXML DEFAULT_MSG RAPIDXML_INCLUDE_DIRS)

mark_as_advanced(RAPIDXML_INCLUDE_DIRS)

if(RAPIDXML_FOUND)
	 add_library(RapidXML::RapidXML INTERFACE IMPORTED)
     set_target_properties(RapidXML::RapidXML PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${RAPIDXML_INCLUDE_DIRS})
endif()
