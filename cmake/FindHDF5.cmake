# Copyright Guillaume Tauzin 2018
# Distributed under the MIT License (MIT)

# .rst:
# FindHDF5
# ------------
#
# Locate the FFT library HDF5
#
#
# Imported targets
# ^^^^^^^^^^^^^^^^
#
# ``HDF5::HDF5``
#
#
# Result variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project:
#
# ``HDF5_FOUND``            - Found the HDF5 library
# ``HDF5_INCLUDE_DIRS``     - HDF5 include directories
#
#
# Cache variables
# ^^^^^^^^^^^^^^^
#
# The following Cache variables may also be set:
#
# ``HDF5_ROOT``             - The root directory of HDF5 installation
#                                 (may also be set as an environment variable)
#

if(NOT HDF5_ROOT AND (DEFINED ENV{HDF5_ROOT}))
  set(HDF5_ROOT $ENV{HDF5_ROOT})
endif()

set(HDF5_INCLUDE_SEARCH_DIRS ${HDF5_ROOT}/include)
set(HDF5_LIBRARY_SEARCH_DIRS ${HDF5_ROOT}/lib)

set(HDF5_KNOWN_VERSIONS "1.8.17")

set(HDF5_PATH_SUFFIXES)
foreach(HDF5_VERSION ${HDF5_KNOWN_VERSIONS})
  list(APPEND HDF5_PATH_SUFFIXES "hdf5-${HDF5_VERSION}")
endforeach()

find_path(HDF5_INCLUDE_DIRS
  NAMES         hdf5.h
  PATHS         ${HDF5_INCLUDE_SEARCH_DIRS}
  PATH_SUFFIXES ${HDF5_PATH_SUFFIXES}
  NO_DEFAULT_PATH)

find_library(HDF5_LIBRARY_DIRS
  NAMES         hdf5
  PATHS         ${HDF5_LIBRARY_SEARCH_DIRS}
  PATH_SUFFIXES ${HDF5_PATH_SUFFIXES}
  NO_DEFAULT_PATH)

message(STATUS "Hey... ${HDF5_ROOT}")
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HDF5 DEFAULT_MSG HDF5_INCLUDE_DIRS HDF5_LIBRARY_DIRS)

mark_as_advanced(HDF5_INCLUDE_DIRS HDF5_LIRARY_DIRS)

if(HDF5_FOUND)
  add_library(HDF5::HDF5 INTERFACE IMPORTED)
  set_target_properties(HDF5::HDF5 PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIRS})
  set_target_properties(HDF5::HDF5 PROPERTIES INTERFACE_LINK_LIBRARIES ${HDF5_LIBRARY_DIRS})
endif()
