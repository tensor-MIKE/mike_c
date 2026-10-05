set(SOURCE_FILES_EC_GENERIC_REF
    ${LVLX_DIR}/basis.c
    ${LVLX_DIR}/ec.c
    ${LVLX_DIR}/normalize.c
    ${LVLX_DIR}/isog_chains.c
)

add_library(${LIB_EC_${SVARIANT_UPPER}} STATIC ${SOURCE_FILES_EC_GENERIC_REF})
target_link_libraries(${LIB_EC_${SVARIANT_UPPER}} ${LIB_PRECOMP_${SVARIANT_UPPER}} ${LIB_GF_${SVARIANT_UPPER}})
target_include_directories(${LIB_EC_${SVARIANT_UPPER}} PRIVATE ${INC_COMMON} ${INC_EC} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_PUBLIC} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}})
target_compile_options(${LIB_EC_${SVARIANT_UPPER}} PRIVATE ${C_OPT_FLAGS})
target_compile_definitions(${LIB_EC_${SVARIANT_UPPER}} PUBLIC MIKE_VARIANT=${SVARIANT_LOWER})

add_subdirectory(test)