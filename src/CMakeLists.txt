add_subdirectory(common)
add_subdirectory(gf)
add_subdirectory(precomp)
add_subdirectory(ec)
add_subdirectory(theta)
add_subdirectory(weil)
add_subdirectory(nike)

FOREACH(SVARIANT ${SVARIANT_S})
    string(TOLOWER ${SVARIANT} SVARIANT_LOWER)
    string(TOUPPER ${SVARIANT} SVARIANT_UPPER)
    set(SOURCE_FILES_VARIANT mike.c)

    # Library for MIKE variant
    add_library(mike_${SVARIANT_LOWER} ${SOURCE_FILES_VARIANT})
    target_link_libraries(mike_${SVARIANT_LOWER} PUBLIC
        ${LIB_NIKE_${SVARIANT_UPPER}}
        ${LIB_GF_${SVARIANT_UPPER}}
        ${LIB_EC_${SVARIANT_UPPER}}
        ${LIB_THETA_${SVARIANT_UPPER}}
        ${LIB_WEIL_${SVARIANT_UPPER}}
        ${LIB_PRECOMP_${SVARIANT_UPPER}}
        mike_common_sys
    )

    target_include_directories(mike_${SVARIANT_LOWER} PUBLIC ${INC_NIKE} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_WEIL} ${INC_THETA} ${INC_EC} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_COMMON} ${INC_HD} ../include PRIVATE common/generic internal)
    target_compile_definitions(mike_${SVARIANT_LOWER} PUBLIC MIKE_VARIANT=${SVARIANT})

    # Library for MIKE variant (test)
    # STATIC: shared test libraries would duplicate the hidden-visibility KAT DRBG state
    add_library(mike_${SVARIANT_LOWER}_test STATIC ${SOURCE_FILES_VARIANT})
    target_link_libraries(mike_${SVARIANT_LOWER}_test PUBLIC
        ${LIB_NIKE_${SVARIANT_UPPER}}
        ${LIB_GF_${SVARIANT_UPPER}}
        ${LIB_EC_${SVARIANT_UPPER}}
        ${LIB_THETA_${SVARIANT_UPPER}}
        ${LIB_WEIL_${SVARIANT_UPPER}}
        ${LIB_PRECOMP_${SVARIANT_UPPER}}
        mike_common_test
    )

    target_include_directories(mike_${SVARIANT_LOWER}_test PUBLIC  ${INC_NIKE} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_WEIL} ${INC_THETA} ${INC_EC} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_COMMON} ${INC_HD} ../include PRIVATE common/generic internal)
    target_compile_definitions(mike_${SVARIANT_LOWER}_test PUBLIC MIKE_VARIANT=${SVARIANT})
ENDFOREACH()