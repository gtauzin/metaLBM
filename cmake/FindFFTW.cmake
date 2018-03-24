# Copyright Guillaume Tauzin 2018
# Distributed under the MIT License (MIT)

# .rst:
# FindFFTW
# ------------
#
# Locate the FFT library FFTW
#
#
# Imported targets
# ^^^^^^^^^^^^^^^^
#
# ``FFTW::FFTW``
#
#
# Result variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project:
#
# ``FFTW_FOUND``            - Found the FFTW library
# ``FFTW_INCLUDE_DIRS``     - FFTW include directories
#
#
# Cache variables
# ^^^^^^^^^^^^^^^
#
# The following Cache variables may also be set:
#
# ``FFTW_ROOT``             - The root directory of FFTW installation
#                                 (may also be set as an environment variable)
#

if(NOT FFTW_ROOT AND (DEFINED ENV{FFTW_ROOT}))
  set(FFTW_ROOT $ENV{FFTW_ROOT})
endif()

set(FFTW_INCLUDE_SEARCH_DIRS ${FFTW_ROOT}/include)
set(FFTW_LIBRARY_SEARCH_DIRS ${FFTW_ROOT}/lib)

set(FFTW_KNOWN_VERSIONS "3.3.7")

set(FFTW_PATH_SUFFIXES)
foreach(FFTW_VERSION ${FFTW_KNOWN_VERSIONS})
  list(APPEND FFTW_PATH_SUFFIXES "fftw-${FFTW_VERSION}")
endforeach()

find_path(FFTW_INCLUDE_DIRS
  NAMES         fftw3.h
  PATHS         ${FFTW_INCLUDE_SEARCH_DIRS}
  PATH_SUFFIXES ${FFTW_PATH_SUFFIXES}
  NO_DEFAULT_PATH)


find_library(FFTW_LIBRARY_DIRS
  NAMES         fftw3
  PATHS         ${FFTW_LIBRARY_SEARCH_DIRS}
  PATH_SUFFIXES ${FFTW_PATH_SUFFIXES}
  NO_DEFAULT_PATH)


find_path(FFTW_MPI_INCLUDE_DIRS
  NAMES         fftw3-mpi.h
  PATHS         ${FFTW_INCLUDE_SEARCH_DIRS}
  PATH_SUFFIXES ${FFTW_PATH_SUFFIXES}
  NO_DEFAULT_PATH)

find_library(FFTW_MPI_LIBRARY_DIRS
  NAMES         fftw3_mpi
  PATHS         ${FFTW_LIBRARY_SEARCH_DIRS}
  PATH_SUFFIXES ${FFTW_PATH_SUFFIXES}
  NO_DEFAULT_PATH)

if(FFTW_WITH_THREADS)
  find_library(FFTW_THREAD_LIBRARY_DIRS
    NAMES         fftw3_threads
    PATHS         ${FFTW_LIBRARY_SEARCH_DIRS}
    PATH_SUFFIXES ${FFTW_PATH_SUFFIXES}
    NO_DEFAULT_PATH)
  set(FFTW_LIBRARY_DIRS ${FFTW_LIBRARY_DIRS} ${FFTW_THREAD_LIBRARY_DIRS})
endif()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_INCLUDE_DIRS FFTW_LIBRARY_DIRS)

mark_as_advanced(FFTW_INCLUDE_DIRS FFTW_LIRARY_DIRS)

if(FFTW_FOUND)
  add_library(FFTW::FFTW INTERFACE IMPORTED)
  set_property(TARGET FFTW::FFTW PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${FFTW_MPI_INCLUDE_DIRS} ${FFTW_INCLUDE_DIRS})
  set_property(TARGET FFTW::FFTW PROPERTY INTERFACE_LINK_LIBRARIES ${FFTW_MPI_LIBRARY_DIRS} ${FFTW_LIBRARY_DIRS})
endif()
