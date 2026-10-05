set(SOURCE_FILES_THETA_GENERIC_REF
    ${LVLX_DIR}/isogeny_dim2.c
    ${LVLX_DIR}/isogeny_dim4.c
    ${LVLX_DIR}/theta_struct_dim2.c
    ${LVLX_DIR}/theta_struct_dim4.c
    ${LVLX_DIR}/isogeny_chain_weil.c
)

add_library(${LIB_THETA_${SVARIANT_UPPER}} STATIC ${SOURCE_FILES_THETA_GENERIC_REF})
target_link_libraries(${LIB_THETA_${SVARIANT_UPPER}} ${LIB_PRECOMP_${SVARIANT_UPPER}} ${LIB_GF_${SVARIANT_UPPER}})
target_include_directories(${LIB_THETA_${SVARIANT_UPPER}} PRIVATE ${INC_COMMON} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_PUBLIC} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_THETA})
target_compile_options(${LIB_THETA_${SVARIANT_UPPER}} PRIVATE ${C_OPT_FLAGS})
target_compile_definitions(${LIB_THETA_${SVARIANT_UPPER}} PUBLIC MIKE_VARIANT=${SVARIANT_LOWER})

add_subdirectory(test)