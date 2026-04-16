#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM) project.
#
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Authors:
#   Pavel Holoborodko <pavel@advanpix.com>
#

cmake_minimum_required(VERSION 3.20)

project(psmm LANGUAGES C CXX VERSION 0.2)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_CXX_STANDARD          17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

option(PSMM_FAST_ESTIMATOR "Use the Mahler-measure early-exit estimator in the main search loop" ON)
option(PSMM_BUILD_TESTS    "Build PSMM unit and integration tests"                                ON)

# ---- Dependencies ---------------------------------------------------------

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

find_package(GMP    REQUIRED)
find_package(GMPXX  REQUIRED)
find_package(NTL    REQUIRED)

# MPSolve is built from our patched upstream via ExternalProject. The
# subdirectory defines the imported target `mpsolve::mps`.
add_subdirectory(mpsolve)

# ---- psmm executable ------------------------------------------------------

add_executable(psmm psmm.cpp)

target_link_libraries(psmm PRIVATE
    mpsolve::mps
    ${NTL_LIBRARIES}
    ${GMPXX_LIBRARIES}
    ${GMP_LIBRARIES}
    Threads::Threads
)

target_include_directories(psmm PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${GMP_INCLUDE_DIRS}
    ${GMPXX_INCLUDE_DIRS}
    ${NTL_INCLUDE_DIRS}
)

if(PSMM_FAST_ESTIMATOR)
    target_compile_definitions(psmm PRIVATE USE_FAST_MAHLER_ESTIMATOR)
endif()

target_compile_options(psmm PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wno-sign-compare -Wno-unused-parameter>
)

install(TARGETS psmm RUNTIME DESTINATION bin)

# ---- Tests ----------------------------------------------------------------

if(PSMM_BUILD_TESTS)
    enable_testing()
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt)
        add_subdirectory(tests)
    endif()
endif()
