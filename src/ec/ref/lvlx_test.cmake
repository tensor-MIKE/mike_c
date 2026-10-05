add_executable(mike_test_ec_${SVARIANT_LOWER} ${LVLX_DIR}/test/test_extras.c ${LVLX_DIR}/test/ec_test.c)
target_link_libraries(mike_test_ec_${SVARIANT_LOWER} ${LIB_PRECOMP_${SVARIANT_UPPER}} ${LIB_EC_${SVARIANT_UPPER}}  mike_common_test)
target_include_directories(mike_test_ec_${SVARIANT_LOWER} PRIVATE ${INC_PUBLIC} ${INC_COMMON} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC})

add_test(mike_test_ec_${SVARIANT_LOWER} mike_test_ec_${SVARIANT_LOWER} 3)

add_executable(mike_test_ec_isog_${SVARIANT_LOWER} ${LVLX_DIR}/test/test_extras.c ${LVLX_DIR}/test/isog_test.c)
target_link_libraries(mike_test_ec_isog_${SVARIANT_LOWER} ${LIB_PRECOMP_${SVARIANT_UPPER}} ${LIB_EC_${SVARIANT_UPPER}}  mike_common_test)
target_include_directories(mike_test_ec_isog_${SVARIANT_LOWER} PRIVATE ${INC_PUBLIC} ${INC_COMMON} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC})

add_test(mike_test_ec_isog_${SVARIANT_LOWER} mike_test_ec_isog_${SVARIANT_LOWER} 3)
