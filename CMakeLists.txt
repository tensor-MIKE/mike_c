# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.13)
project(MIKE VERSION 0.1 LANGUAGES C ASM)

set(MIKE_SO_VERSION "0")
set(CMAKE_C_STANDARD 11)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CTest)

option(ENABLE_STRICT "Build with strict compile options." ON)
option(ENABLE_TESTS  "Enable compilation of tests." ON)
option(ENABLE_CT_TESTING  "Enable compilation for constant time testing." OFF)
set(GF_RADIX "AUTO" CACHE STRING "Set the radix for the gf module (currently supported values: 32 or 64), or AUTO.")

if (NOT DEFINED MIKE_BUILD_TYPE)
  SET(MIKE_BUILD_TYPE "ref")
endif()

if (${MIKE_BUILD_TYPE} MATCHES "ref")
  add_compile_definitions(MIKE_BUILD_TYPE_REF MIKE_GF_IMPL_REF)
elseif (${MIKE_BUILD_TYPE} MATCHES "opt")
  add_compile_definitions(MIKE_BUILD_TYPE_OPT MIKE_GF_IMPL_REF)
elseif (${MIKE_BUILD_TYPE} MATCHES "broadwell")
  add_compile_definitions(MIKE_BUILD_TYPE_BROADWELL MIKE_GF_IMPL_SAT64)
elseif (${MIKE_BUILD_TYPE} STREQUAL "arm64")
  add_compile_definitions(MIKE_BUILD_TYPE_ARM64 MIKE_GF_IMPL_SAT64)
elseif (${MIKE_BUILD_TYPE} MATCHES "arm64crypto")
  add_compile_definitions(MIKE_BUILD_TYPE_ARM64CRYPTO MIKE_GF_IMPL_REF)
endif()

if (NOT DEFINED MIKE_PRIME_CHOICE)
  SET(MIKE_PRIME_CHOICE "fast")
endif()

SET(VARIANT_ONE "p308_633;p474_593;p628_317")
SET(VARIANT_TWO "p374_117;p566_77;p758_41")

if (${MIKE_PRIME_CHOICE} MATCHES "fast")
  SET(SVARIANT_S ${VARIANT_ONE})
elseif (${MIKE_PRIME_CHOICE} MATCHES "conservative")
  SET(SVARIANT_S ${VARIANT_TWO})
endif()

include(.cmake/flags.cmake)
include(.cmake/sanitizers.cmake)
include(.cmake/target.cmake)

set(BM_BINS "" CACHE INTERNAL "List of benchmark binaries")

set(SELECT_IMPL_TYPE ${PROJECT_SOURCE_DIR}/.cmake/impl_type.cmake)
set(SELECT_MIKE_VARIANT ${PROJECT_SOURCE_DIR}/.cmake/mike_variant.cmake)

set(INC_PUBLIC ${PROJECT_SOURCE_DIR}/include)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(apps)

include(.cmake/bm.cmake)

#endif()
