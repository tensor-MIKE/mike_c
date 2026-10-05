add_executable(mike_test_gf_${SVARIANT_LOWER}_fp ${${CCSD_NAME_UPPER}_GENERIC_DIR}/test/test_fp.c ${${CCSD_NAME_UPPER}_GENERIC_DIR}/test/test_utils.c)
target_link_libraries(mike_test_gf_${SVARIANT_LOWER}_fp ${LIB_GF_${SVARIANT_UPPER}} ${LIB_PRECOMP_${SVARIANT_UPPER}} mike_common_test)
target_include_directories(mike_test_gf_${SVARIANT_LOWER}_fp PRIVATE ${INC_GF} ${INC_COMMON} ${PROJECT_SOURCE_DIR}/src/precomp/ref/${SVARIANT_LOWER}/include ${INC_PUBLIC})

add_executable(mike_test_gf_${SVARIANT_LOWER}_fp2 ${${CCSD_NAME_UPPER}_GENERIC_DIR}/test/test_fp2.c ${${CCSD_NAME_UPPER}_GENERIC_DIR}/test/test_utils.c)
target_link_libraries(mike_test_gf_${SVARIANT_LOWER}_fp2 ${LIB_GF_${SVARIANT_UPPER}} ${LIB_PRECOMP_${SVARIANT_UPPER}} mike_common_test)
target_include_directories(mike_test_gf_${SVARIANT_LOWER}_fp2 PRIVATE ${INC_GF} ${INC_COMMON} ${PROJECT_SOURCE_DIR}/src/precomp/ref/${SVARIANT_LOWER}/include ${INC_PUBLIC})

add_test(mike_test_gf_${SVARIANT_LOWER}_fp mike_test_gf_${SVARIANT_LOWER}_fp test)
add_test(mike_test_gf_${SVARIANT_LOWER}_fp2 mike_test_gf_${SVARIANT_LOWER}_fp2 test)

add_executable(mike_bench_gf_${SVARIANT_LOWER}_fp ${${CCSD_NAME_UPPER}_GENERIC_DIR}/test/bench_fp.c ${${CCSD_NAME_UPPER}_GENERIC_DIR}/test/test_utils.c)
target_link_libraries(mike_bench_gf_${SVARIANT_LOWER}_fp ${LIB_GF_${SVARIANT_UPPER}} ${LIB_PRECOMP_${SVARIANT_UPPER}} mike_common_sys)
target_include_directories(mike_bench_gf_${SVARIANT_LOWER}_fp PRIVATE ${INC_GF} ${INC_COMMON} ${PROJECT_SOURCE_DIR}/src/precomp/ref/${SVARIANT_LOWER}/include ${INC_PUBLIC})

add_executable(mike_bench_gf_${SVARIANT_LOWER}_fp2 ${${CCSD_NAME_UPPER}_GENERIC_DIR}/test/bench_fp2.c ${${CCSD_NAME_UPPER}_GENERIC_DIR}/test/test_utils.c)
target_link_libraries(mike_bench_gf_${SVARIANT_LOWER}_fp2 ${LIB_GF_${SVARIANT_UPPER}} ${LIB_PRECOMP_${SVARIANT_UPPER}} mike_common_sys)
target_include_directories(mike_bench_gf_${SVARIANT_LOWER}_fp2 PRIVATE ${INC_GF} ${INC_COMMON} ${PROJECT_SOURCE_DIR}/src/precomp/ref/${SVARIANT_LOWER}/include ${INC_PUBLIC})

set(BM_BINS ${BM_BINS}
    mike_bench_gf_${SVARIANT_LOWER}_fp mike_bench_gf_${SVARIANT_LOWER}_fp2
    CACHE INTERNAL "List of benchmark executables")

