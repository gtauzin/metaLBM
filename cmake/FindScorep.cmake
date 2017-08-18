# Copyright (c) 2016, Technische Universität Dresden, Germany
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted
# provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions
#    and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
#    and the following disclaimer in the documentation and/or other materials provided with the
#    distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
#    or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
# THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

if(SCOREP_CONFIG_PATH)
  find_program(SCOREP_CONFIG NAMES scorep-config
    PATHS /opt/scorep/bin)
else(SCOREP_CONFIG_PATH)
  find_program(SCOREP_CONFIG NAMES scorep-config
    PATHS /opt/scorep/bin)
endif(SCOREP_CONFIG_PATH)

if(NOT SCOREP_CONFIG)
  message(STATUS "no scorep-config found")
  set(SCOREP_FOUND false)
else(NOT SCOREP_CONFIG)

    message(STATUS "SCOREP library found. (using ${SCOREP_CONFIG})")

    execute_process(COMMAND ${SCOREP_CONFIG} "--user" "--cppflags" OUTPUT_VARIABLE SCOREP_INCLUDE_DIRS)
    string(REPLACE "-I" ";" SCOREP_INCLUDE_DIRS ${SCOREP_INCLUDE_DIRS})

    execute_process(COMMAND ${SCOREP_CONFIG} "--user" "--ldflags" OUTPUT_VARIABLE _LINK_LD_ARGS)
    string( REPLACE " " ";" _LINK_LD_ARGS ${_LINK_LD_ARGS} )
    foreach( _ARG ${_LINK_LD_ARGS} )
        if(${_ARG} MATCHES "^-L")
          string(REGEX REPLACE "^-L" "" _ARG ${_ARG})
          set(SCOREP_LINK_DIRS ${SCOREP_LINK_DIRS} ${_ARG})
        endif(${_ARG} MATCHES "^-L")
    endforeach(_ARG)

    execute_process(COMMAND ${SCOREP_CONFIG} "--user" "--libs" OUTPUT_VARIABLE _LINK_LD_ARGS)
    string( REPLACE " " ";" _LINK_LD_ARGS ${_LINK_LD_ARGS} )
    foreach( _ARG ${_LINK_LD_ARGS} )
      if(${_ARG} MATCHES "^-l")
        string(REGEX REPLACE "^-l" "" _ARG ${_ARG})
        find_library(_SCOREP_LIB_FROM_ARG NAMES ${_ARG}
          PATHS ${SCOREP_LINK_DIRS})
        if(_SCOREP_LIB_FROM_ARG)
          set(SCOREP_LIBRARIES ${SCOREP_LIBRARIES} ${_SCOREP_LIB_FROM_ARG})
        endif(_SCOREP_LIB_FROM_ARG)
        unset(_SCOREP_LIB_FROM_ARG CACHE)
        endif(${_ARG} MATCHES "^-l")
      endforeach(_ARG)

    set(SCOREP_FOUND true)
endif(NOT SCOREP_CONFIG)

mark_as_advanced(SCOREP_CONFIG)
