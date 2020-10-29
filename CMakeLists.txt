/*
 * This file is part of "Polynomials with Small Mahler Measure" (PSMM) project.
 *
 * Copyright (C) 2020, Advanpix LLC.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Pavel Holoborodko <pavel@advanpix.com>
 */

cmake_minimum_required(VERSION 3.1)

project(psmm LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_VERBOSE_MAKEFILE on)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build" FORCE)
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(SOURCES psmm.cpp)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

if(NOT GMP_INCLUDE_DIRS OR NOT GMP_LIBRARIES)
    find_package(GMP REQUIRED)
endif()

if(NOT GMPXX_INCLUDE_DIRS OR NOT GMPXX_LIBRARIES)
    find_package(GMPXX REQUIRED)
endif()

if(NOT NTL_INCLUDE_DIRS OR NOT NTL_LIBRARIES)
    find_package(NTL REQUIRED)
endif()

if(NOT MPSOLVE_INCLUDE_DIRS OR NOT MPSOLVE_LIBRARIES)
    find_package(MPSOLVE 3.2.1 REQUIRED)
endif()

add_executable(psmm ${SOURCES})

target_link_libraries(psmm PUBLIC
    ${NTL_LIBRARIES}
    ${MPSOLVE_LIBRARIES}
    ${GMP_LIBRARIES}
    ${GMPXX_LIBRARIES}
    Threads::Threads
)

target_include_directories(psmm PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${GMP_INCLUDE_DIRS}
    ${GMPXX_INCLUDE_DIRS}
    ${NTL_INCLUDE_DIRS}
    ${MPSOLVE_INCLUDE_DIRS}
)
