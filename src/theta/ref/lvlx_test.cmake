add_executable(mike_test_theta_${SVARIANT_LOWER} ${LVLX_DIR}/test/test_theta.c)
target_link_libraries(mike_test_theta_${SVARIANT_LOWER} ${LIB_THETA_${SVARIANT_UPPER}}  mike_common_test)
target_include_directories(mike_test_theta_${SVARIANT_LOWER} PRIVATE ${INC_PUBLIC} ${INC_COMMON} ${INC_THETA} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}})

add_test(mike_test_theta_${SVARIANT_LOWER} mike_test_theta_${SVARIANT_LOWER} 3)


