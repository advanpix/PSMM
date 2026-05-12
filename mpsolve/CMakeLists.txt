#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
# Build a patched MPSolve as an ExternalProject and expose it to PSMM
# as an imported target `mpsolve::mps`.
#
# The MPSolve project ships an autotools build, not CMake. We therefore:
#   1. Download the upstream source tarball pinned by upstream-version.cmake.
#   2. Apply the patches in mpsolve/patches/ (in sorted order).
#   3. Run autogen.sh / configure / make / make install into a build-local
#      prefix under ${CMAKE_BINARY_DIR}/mpsolve-install.
#   4. Export an IMPORTED static library target that PSMM links against.
#
# Build-host prerequisites: autoconf, automake, libtool, bison, flex,
# pkg-config, gmp headers.
#

include(ExternalProject)
include(${CMAKE_CURRENT_SOURCE_DIR}/upstream-version.cmake)

set(_mpsolve_install  ${CMAKE_BINARY_DIR}/mpsolve-install)
set(_mpsolve_include  ${_mpsolve_install}/include)
set(_mpsolve_libdir   ${_mpsolve_install}/lib)
set(_mpsolve_libname  ${CMAKE_STATIC_LIBRARY_PREFIX}mps${CMAKE_STATIC_LIBRARY_SUFFIX})
set(_mpsolve_libfile  ${_mpsolve_libdir}/${_mpsolve_libname})

file(GLOB _mpsolve_patches CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/patches/*.patch")
list(SORT _mpsolve_patches)

# Build a shell snippet that applies every patch in order from inside
# <SOURCE_DIR>. Using `patch -p1 -N` makes re-builds idempotent.
set(_apply_patches_cmd "true")
foreach(_p ${_mpsolve_patches})
    set(_apply_patches_cmd "${_apply_patches_cmd} && patch -p1 -N -r - < ${_p}")
endforeach()

# MPSolve's configure doesn't accept --disable-ui cleanly in all versions,
# but we can at least strip the GUI/examples/docs and pick static output.
set(_mpsolve_configure_args
    --prefix=${_mpsolve_install}
    --disable-shared
    --enable-static
    --disable-examples
    --disable-ui
    --disable-documentation
    --without-octave
    --without-matlab
    --without-ui
)

ExternalProject_Add(mpsolve_ext
    PREFIX            ${CMAKE_BINARY_DIR}/mpsolve-ext
    URL               ${PSMM_MPSOLVE_URL}
    URL_HASH          ${PSMM_MPSOLVE_URL_HASH}
    DOWNLOAD_NO_PROGRESS      ON
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE

    PATCH_COMMAND     sh -c "${_apply_patches_cmd}"

    # autogen.sh regenerates configure from configure.ac — needed because
    # the tarball ships without a prebuilt configure script.
    CONFIGURE_COMMAND sh -c "cd <SOURCE_DIR> && ./autogen.sh"
              COMMAND <SOURCE_DIR>/configure ${_mpsolve_configure_args}

    BUILD_COMMAND     make -j
    INSTALL_COMMAND   make install

    BUILD_BYPRODUCTS  ${_mpsolve_libfile}
    LOG_DOWNLOAD      ON
    LOG_PATCH         ON
    LOG_CONFIGURE     ON
    LOG_BUILD         ON
    LOG_INSTALL       ON
)

file(MAKE_DIRECTORY ${_mpsolve_include})

add_library(mpsolve::mps STATIC IMPORTED GLOBAL)
set_target_properties(mpsolve::mps PROPERTIES
    IMPORTED_LOCATION             ${_mpsolve_libfile}
    INTERFACE_INCLUDE_DIRECTORIES ${_mpsolve_include}
    INTERFACE_LINK_LIBRARIES      "Threads::Threads;${GMP_LIBRARIES};m"
)
add_dependencies(mpsolve::mps mpsolve_ext)

set(MPSOLVE_VERSION      ${PSMM_MPSOLVE_VERSION} PARENT_SCOPE)
set(MPSOLVE_INCLUDE_DIRS ${_mpsolve_include}     PARENT_SCOPE)
set(MPSOLVE_LIBRARIES    mpsolve::mps            PARENT_SCOPE)
