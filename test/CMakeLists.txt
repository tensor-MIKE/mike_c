# SPDX-License-Identifier: Apache-2.0

FOREACH(SVARIANT ${SVARIANT_S})

    string(TOLOWER ${SVARIANT} SVARIANT_LOWER)
    string(TOUPPER ${SVARIANT} SVARIANT_UPPER)

    add_executable(mike_test_scheme_${SVARIANT} test_mike.c )
    target_link_libraries(mike_test_scheme_${SVARIANT} mike_${SVARIANT_LOWER} )
    target_include_directories(mike_test_scheme_${SVARIANT} PUBLIC ${PROJECT_SOURCE_DIR}/src/common ${INC_PUBLIC})
    add_test(mike_${SVARIANT}_SELFTEST mike_test_scheme_${SVARIANT})

ENDFOREACH()



FOREACH(SVARIANT ${SVARIANT_S})

    string(TOLOWER ${SVARIANT} SVARIANT_LOWER)
    string(TOUPPER ${SVARIANT} SVARIANT_UPPER)

    # Benchmarking tool
    add_executable(benchmark_${SVARIANT} benchmark.c)
    target_link_libraries(benchmark_${SVARIANT} mike_${SVARIANT_LOWER} )
    target_include_directories(benchmark_${SVARIANT} PUBLIC ${PROJECT_SOURCE_DIR}/src/common ${INC_PUBLIC})
    target_compile_definitions(benchmark_${SVARIANT_LOWER} PUBLIC VARIANT=${SVARIANT_LOWER})

ENDFOREACH()


