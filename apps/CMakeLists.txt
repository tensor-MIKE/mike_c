# SPDX-License-Identifier: Apache-2.0

FOREACH(SVARIANT ${SVARIANT_S})

    string(TOLOWER ${SVARIANT} SVARIANT_LOWER)
    string(TOUPPER ${SVARIANT} SVARIANT_UPPER)

    # Benchmarking tool
    add_executable(mike_cli_${SVARIANT} cli.c)
    target_link_libraries(mike_cli_${SVARIANT} mike_${SVARIANT_LOWER} )
    target_include_directories(mike_cli_${SVARIANT} PUBLIC ${PROJECT_SOURCE_DIR}/src/common ${INC_PUBLIC})
    target_compile_definitions(mike_cli_${SVARIANT_LOWER} PUBLIC VARIANT=${SVARIANT_LOWER})

ENDFOREACH()