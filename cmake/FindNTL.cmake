# Try to find the NTL library
# https://www.shoup.net/ntl/
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(NTL 11.4.3)
# to require version 11.4.3 to newer of GMP.
#
# Once done this will define
#
#  NTL_FOUND - system has NTL lib with correct version
#  NTL_INCLUDE_DIRS - the NTL include directory
#  NTL_LIBRARIES - the NTL library
#  NTL_VERSION - NTL version
#
# Copyright (c) 2019 Alex Brandt, <abrandt5@uwo.ca>
# Copyright (c) 2020 Pavel Holoborodko, <pavel@advanpix.com>
#
# Redistribution and use is allowed according to the terms of the BSD license.

find_path(NTL_INCLUDE_DIRS NAMES NTL/version.h PATHS $ENV{NTLDIR} ${INCLUDE_INSTALL_DIR})

# Set NTL_FIND_VERSION to 11.4.0 if no minimum version is specified
if(NOT NTL_FIND_VERSION)
  if(NOT NTL_FIND_VERSION_MAJOR)
    set(NTL_FIND_VERSION_MAJOR 11)
  endif()
  if(NOT NTL_FIND_VERSION_MINOR)
    set(NTL_FIND_VERSION_MINOR 4)
  endif()
  if(NOT NTL_FIND_VERSION_PATCH)
    set(NTL_FIND_VERSION_PATCH 0)
  endif()
  set(NTL_FIND_VERSION
    "${NTL_FIND_VERSION_MAJOR}.${NTL_FIND_VERSION_MINOR}.${NTL_FIND_VERSION_PATCH}")
endif()

if(NTL_INCLUDE_DIRS)
  file(GLOB NTL_HEADERS "${NTL_INCLUDE_DIRS}/NTL/version.h")
  foreach(ntl_header_filename ${NTL_HEADERS})
    file(READ "${ntl_header_filename}" _ntl_version_header)
    string(REGEX MATCH
      "define[ \t]*NTL_VERSION[ \t]*\"*([0-9]+\\.[0-9]+\\.[0-9]+)\"*" _ntl_version_match
      "${_ntl_version_header}")
    if(_ntl_version_match)
      set(NTL_VERSION "${CMAKE_MATCH_1}")
      break()
    endif()
  endforeach()

  # Check whether found version exists and exceeds the minimum requirement
  if(NOT NTL_VERSION)
    set(NTL_VERSION_OK FALSE)
    message(STATUS "NTL version was not detected")
  elseif(${NTL_VERSION} VERSION_LESS ${NTL_FIND_VERSION})
    set(NTL_VERSION_OK FALSE)
    message(STATUS "NTL version ${NTL_VERSION} found in ${NTL_INCLUDE_DIRS}, "
                   "but at least version ${NTL_FIND_VERSION} is required")
  else()
    set(NTL_VERSION_OK TRUE)
  endif()
endif()


find_library(NTL_LIBRARIES ntl PATHS $ENV{NTLDIR} ${LIB_INSTALL_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NTL DEFAULT_MSG NTL_INCLUDE_DIRS NTL_LIBRARIES NTL_VERSION_OK)

mark_as_advanced(NTL_INCLUDE_DIRS NTL_LIBRARIES)
