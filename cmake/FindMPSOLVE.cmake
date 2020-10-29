# Try to find the MPSolve library
# https://github.com/robol/MPSolve
#
#
# Once done this will define
#
#  MPSOLVE_FOUND - system has MPSOLVE lib with correct version
#  MPSOLVE_INCLUDE_DIRS - the MPSOLVE include directory
#  MPSOLVE_LIBRARIES - the MPSOLVE library
#
# Copyright (c) 2019 Alex Brandt, <abrandt5@uwo.ca>
# Copyright (c) 2020 Pavel Holoborodko, <pavel@advanpix.com>
#
# Redistribution and use is allowed according to the terms of the BSD license.

find_path(MPSOLVE_INCLUDE_DIRS NAMES mps/mps.h PATHS $ENV{MPSOLVEDIR} ${INCLUDE_INSTALL_DIR})

if(NOT MPSOLVE_FIND_VERSION)
  if(NOT MPSOLVE_FIND_VERSION_MAJOR)
    set(MPSOLVE_FIND_VERSION_MAJOR 3)
  endif()
  if(NOT MPSOLVE_FIND_VERSION_MINOR)
    set(MPSOLVE_FIND_VERSION_MINOR 2)
  endif()
  if(NOT MPSOLVE_FIND_VERSION_PATCH)
    set(MPSOLVE_FIND_VERSION_PATCH 1)
  endif()
  set(MPSOLVE_FIND_VERSION
    "${MPSOLVE_FIND_VERSION_MAJOR}.${MPSOLVE_FIND_VERSION_MINOR}.${MPSOLVE_FIND_VERSION_PATCH}")
endif()

if(MPSOLVE_INCLUDE_DIRS)
  file(GLOB MPSOLVE_HEADERS "${MPSOLVE_INCLUDE_DIRS}/version.h")
  foreach(mps_header_filename ${MPSOLVE_HEADERS})
    file(READ "${mps_header_filename}" _mps_version_header)
    string(REGEX MATCH
      "This file is part of MPSolve ([0-9]+\\.[0-9]+\\.?[0-9]*)" _mps_version_match
      "${_mps_version_header}")
    if(_mps_version_match)
      set(MPSOLVE_VERSION "${CMAKE_MATCH_1}")
      #message("MPSOLVE_VERSION=${MPSOLVE_VERSION}")
      break()
    endif()
    string(REGEX MATCH
      "\\*\\*.*Version ([0-9]+\\.[0-9]+)" _mps_version_match
      "${_mps_version_header}")
    if(_mps_version_match)
      set(MPSOLVE_VERSION "${CMAKE_MATCH_1}")
      #message("MPSOLVE_VERSION=${MPSOLVE_VERSION}")
      break()
    endif()
  endforeach()

  # Check whether found version exists and exceeds the minimum requirement
  if(NOT MPSOLVE_VERSION)
    set(MPSOLVE_VERSION_OK FALSE)
    message(STATUS "MPSOLVE version was not detected")
  elseif(${MPSOLVE_VERSION} VERSION_LESS ${MPSOLVE_FIND_VERSION})
    set(MPSOLVE_VERSION_OK FALSE)
    message(STATUS "MPSOLVE version ${MPSOLVE_VERSION} found in ${MPSOLVE_INCLUDE_DIRS}, "
                   "but at least version ${MPSOLVE_FIND_VERSION} is required")
  else()
    set(MPSOLVE_VERSION_OK TRUE)
  endif()
endif()


find_library(MPSOLVE_LIBRARIES mps PATHS $ENV{MPSOLVEDIR} ${LIB_INSTALL_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPSOLVE DEFAULT_MSG
                                  MPSOLVE_INCLUDE_DIRS MPSOLVE_LIBRARIES MPSOLVE_VERSION_OK)

mark_as_advanced(MPSOLVE_INCLUDE_DIRS MPSOLVE_LIBRARIES)
