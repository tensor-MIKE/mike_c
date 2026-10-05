add_executable(mike_test_nike_${SVARIANT_LOWER} ${LVLX_DIR}/test/test_nike.c)
target_link_libraries(mike_test_nike_${SVARIANT_LOWER} ${LIB_EC_${SVARIANT_UPPER}} ${LIB_WEIL_${SVARIANT_UPPER}} ${LIB_THETA_${SVARIANT_UPPER}} ${LIB_GF_${SVARIANT_UPPER}} ${LIB_NIKE_${SVARIANT_UPPER}}  mike_common_test)
target_include_directories(mike_test_nike_${SVARIANT_LOWER} PRIVATE ${INC_PUBLIC} ${INC_COMMON} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_GF} ${INC_WEIL} ${INC_THETA} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_NIKE})

add_test(mike_test_nike_${SVARIANT_LOWER} mike_test_nike_${SVARIANT_LOWER} 3)

add_custom_command(
  TARGET mike_test_nike_${SVARIANT_LOWER}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy $<TARGET_FILE:mike_test_nike_${SVARIANT_LOWER}> ${CMAKE_BINARY_DIR}/test/mike_${SVARIANT_LOWER}
)

add_executable(mike_bench_nike_${SVARIANT_LOWER} ${LVLX_DIR}/test/bench_nike.c)
target_link_libraries(mike_bench_nike_${SVARIANT_LOWER} ${LIB_EC_${SVARIANT_UPPER}} ${LIB_WEIL_${SVARIANT_UPPER}} ${LIB_THETA_${SVARIANT_UPPER}} ${LIB_GF_${SVARIANT_UPPER}} ${LIB_NIKE_${SVARIANT_UPPER}}  mike_common_sys)
target_include_directories(mike_bench_nike_${SVARIANT_LOWER} PRIVATE ${INC_PUBLIC} ${INC_COMMON} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_GF} ${INC_WEIL} ${INC_THETA} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_NIKE})

set(BM_BINS ${BM_BINS} mike_bench_nike_${SVARIANT_LOWER} CACHE INTERNAL "List of benchmark executables")
