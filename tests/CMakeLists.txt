#
# PSMM test suite.
#
# Uses doctest (header-only) pulled via FetchContent. Tests compile
# against the PSMM headers directly and the patched libmps built by
# ../mpsolve.
#

include(FetchContent)

FetchContent_Declare(
    doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG        v2.4.11
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(doctest)

#
# One test binary per thematic group; each gets registered as a CTest
# case so `ctest` / `ctest -V` runs them all.
#

function(psmm_add_test name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE
        doctest::doctest
        mpsolve::mps
        ${NTL_LIBRARIES}
        ${GMPXX_LIBRARIES}
        ${GMP_LIBRARIES}
        Threads::Threads
    )
    target_include_directories(${name} PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${GMP_INCLUDE_DIRS}
        ${GMPXX_INCLUDE_DIRS}
        ${NTL_INCLUDE_DIRS}
    )
    target_compile_definitions(${name} PRIVATE
        PSMM_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data"
        USE_FAST_MAHLER_ESTIMATOR
    )
    add_test(NAME ${name} COMMAND ${name})
endfunction()

psmm_add_test(test_mahler_reference)
psmm_add_test(test_mahler_early_exit)
psmm_add_test(test_iterator)
