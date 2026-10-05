set(SOURCE_FILES_COMMON_GENERIC
    randombytes_system.c 
    fips202.c 
    mem.c
    tools.c
)

foreach (MIKE_COMMON_TARGET mike_common_test mike_common_sys)
    add_library(${MIKE_COMMON_TARGET} STATIC ${SOURCE_FILES_COMMON_GENERIC})
    target_include_directories(${MIKE_COMMON_TARGET} PRIVATE include ${INC_PUBLIC})
    target_compile_options(${MIKE_COMMON_TARGET} PUBLIC ${C_OPT_FLAGS})
    if (ENABLE_CT_TESTING)
        target_compile_definitions(${MIKE_COMMON_TARGET} PUBLIC ENABLE_CT_TESTING)
    endif()
endforeach()
