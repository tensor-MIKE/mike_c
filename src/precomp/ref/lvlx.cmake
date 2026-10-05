set(SOURCE_FILES_PRECOMP_${SVARIANT_UPPER}_REF
    ec_params.c
    e0_basis.c
    constants.c
)

add_library(${LIB_PRECOMP_${SVARIANT_UPPER}} STATIC ${SOURCE_FILES_PRECOMP_${SVARIANT_UPPER}_REF})
target_link_libraries(${LIB_PRECOMP_${SVARIANT_UPPER}})
target_include_directories(${LIB_PRECOMP_${SVARIANT_UPPER}} PRIVATE ${INC_PUBLIC} ${INC_EC} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_COMMON})
#target_include_directories(${INC_EC} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_COMMON})
target_compile_options(${LIB_PRECOMP_${SVARIANT_UPPER}} PRIVATE ${C_OPT_FLAGS})
target_compile_definitions(${LIB_PRECOMP_${SVARIANT_UPPER}} PUBLIC MIKE_VARIANT=${SVARIANT_LOWER})

set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM true)