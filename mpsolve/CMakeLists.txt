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
# <SOURCE_DIR>. We use `patch -p1 -N -r -` to forward-apply each patch
# while ignoring hunks that are already in place. GNU patch returns:
#   0 — every hunk applied,
#   1 — at least one hunk was skipped (already-applied or rejected),
#   2 — fatal error (can't find file to patch, malformed patch, etc.).
# Re-runs of the patch step legitimately hit case 1 (every hunk is
# already in the source tree from the previous run), so we treat exit
# code 1 as success and only fail on 2+. Without this, retriggering the
# patch step via ExternalProject_Add_StepDependencies (e.g. after
# editing a patch file) would fail spuriously the second time around.
set(_apply_patches_cmd "true")
foreach(_p ${_mpsolve_patches})
    set(_apply_patches_cmd "${_apply_patches_cmd} && ( patch -p1 -N -r - < ${_p} || [ $? -eq 1 ] )")
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

# Invalidate the download step (which performs extract) whenever any
# patch file changes. Without this, ExternalProject caches the patch
# stamp once and never re-applies — so editing or adding a patch
# would silently be ignored. CONFIGURE_DEPENDS on the GLOB above
# only triggers cmake reconfigure; it does not by itself invalidate
# the patch stamp.
#
# Why depend on `download`, not `patch`? GNU `patch -p1 -N` is not
# truly idempotent: when re-run against an already-patched tree, its
# fuzz-match heuristic can find spurious second locations to apply
# (we observed patch 0001 silently producing duplicate function
# definitions on re-run). Forcing the download step to re-run gives
# us a freshly-extracted pristine source tree, on which all patches
# apply cleanly exactly once. Re-extract is essentially free —
# URL_HASH ensures the cached tarball is reused without a network
# fetch.
ExternalProject_Add_StepDependencies(mpsolve_ext download ${_mpsolve_patches})

# Convenience target: force a full MPSolve rebuild from the cached
# tarball. Use this if the source/install tree gets into a bad state,
# or as a definitive "the patches really did re-apply" reset. Running
# `make mpsolve_force_rebuild` is equivalent to deleting
# build/mpsolve-ext + build/mpsolve-install and rebuilding.
#
# IMPORTANT (CMake gotcha): ExternalProject_Add bundles extract into
# the download step, so deleting only the extracted source dir
# (mpsolve-ext/src/mpsolve_ext) is not enough -- the mpsolve_ext-download
# stamp must also be invalidated, otherwise the next build will jump
# straight to the patch step against an empty source tree and fail
# with "can't find file to patch". This target invalidates both.
add_custom_target(mpsolve_force_rebuild
    COMMAND ${CMAKE_COMMAND} -E rm -rf
            ${CMAKE_BINARY_DIR}/mpsolve-ext
            ${CMAKE_BINARY_DIR}/mpsolve-install
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target mpsolve_ext
    COMMENT "Force full MPSolve re-download + re-extract + re-patch + rebuild"
    USES_TERMINAL
    VERBATIM
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
