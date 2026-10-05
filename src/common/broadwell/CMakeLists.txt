set(SOURCE_FILES_COMMON_AESNI
    aes_ni.c
    ctr_drbg.c
    randombytes_ctrdrbg_aesni.c
    vaes256_key_expansion.S
)

foreach(MIKE_COMMON_TARGET mike_common_test mike_common_sys)
    target_sources(${MIKE_COMMON_TARGET} PRIVATE ${SOURCE_FILES_COMMON_AESNI})
    target_include_directories(${MIKE_COMMON_TARGET} PRIVATE include)
    target_compile_definitions(${MIKE_COMMON_TARGET} PRIVATE RANDOMBYTES_AES_NI)
    target_compile_options(${MIKE_COMMON_TARGET} PRIVATE -maes -mavx2)
endforeach()

set(SOURCE_FILES_CTRDRBG_TEST_BENCHMARK
    ../ref/aes_c.c
    aes_ni.c
    ctr_drbg.c
    randombytes_ctrdrbg_aesni.c
    ../ref/randombytes_ctrdrbg.c
    ../generic/randombytes_system.c
    vaes256_key_expansion.S
)

add_executable(mike_test_ctrdrbg_intel ${SOURCE_FILES_CTRDRBG_TEST_BENCHMARK} ../generic/test/test_ctrdrbg.c)
target_include_directories(mike_test_ctrdrbg_intel PRIVATE ${INC_PUBLIC} ${INC_COMMON} include ../ref/include)
target_compile_definitions(mike_test_ctrdrbg_intel PRIVATE
    CTRDRBG_TEST_BENCH
    RANDOMBYTES_INIT_PLATFORM=randombytes_init_aes_ni
    RANDOMBYTES_PLATFORM=randombytes_aes_ni)
target_compile_options(mike_test_ctrdrbg_intel PRIVATE -maes -mavx2)

add_test(mike_test_ctrdrbg_intel mike_test_ctrdrbg_intel)

add_executable(mike_bench_ctrdrbg_intel ${SOURCE_FILES_CTRDRBG_TEST_BENCHMARK} ../generic/test/bench_ctrdrbg.c)
target_include_directories(mike_bench_ctrdrbg_intel PRIVATE ${INC_PUBLIC} ${INC_COMMON} include ../ref/include)
target_compile_definitions(mike_bench_ctrdrbg_intel PRIVATE
    CTRDRBG_TEST_BENCH
    RANDOMBYTES_INIT_PLATFORM=randombytes_init_aes_ni
    RANDOMBYTES_PLATFORM=randombytes_aes_ni)
target_compile_options(mike_bench_ctrdrbg_intel PRIVATE -maes -mavx2)

set(BM_BINS ${BM_BINS} mike_bench_ctrdrbg_intel CACHE INTERNAL "List of benchmark executables")
