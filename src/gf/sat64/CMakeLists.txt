# The saturated 64-bit backend shares its C layer (headers, lvlx/ and the
# per-prime generated fp.c) across architectures; only <prime>/<arch>/fp_asm.S
# differs.
set(LVLX_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lvlx)

if (MIKE_BUILD_TYPE MATCHES "broadwell")
    set(GF_SAT64_ARCH broadwell)
    set(GF_SAT64_FLAGS -mbmi2 -madx -fno-tree-vectorize -fno-tree-slp-vectorize -fno-stack-protector)
elseif (MIKE_BUILD_TYPE STREQUAL "arm64")
    set(GF_SAT64_ARCH arm64)
    set(GF_SAT64_FLAGS "")
else()
    message(FATAL_ERROR "No sat64 assembly for build type ${MIKE_BUILD_TYPE}")
endif()

include(${SELECT_MIKE_VARIANT})
