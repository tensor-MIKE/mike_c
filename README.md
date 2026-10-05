# MIKE

This library is a C implementation of MIKE.

## Requirements

- CMake (version 3.13 or later)
- C11-compatible compiler


## Build

For a generic build
```
$ mkdir -p build
$ cd build
$ cmake -DMIKE_BUILD_TYPE=ref -DMIKE_PRIME_CHOICE=fast ..
$ make
$ make test
```

An optimized executable with debug code and assertions disabled can be built
replacing the `cmake` command above by
```
cmake -DMIKE_BUILD_TYPE=ref -DMIKE_PRIME_CHOICE=<fast/conservative> -DCMAKE_BUILD_TYPE=Release ..
```

## Build options

CMake build options can be specified with `-D<BUILD_OPTION>=<VALUE>`.

### MIKE_BUILD_TYPE

Specifies the build type for which MIKE is built. The currently supported values are:
- `ref`: builds the plain reference implementation (pure C).
- `broadwell`: builds an additional assembly optimized implementation targeting the Intel Broadwell architecture (and later). The optimizations are applied to the finite field arithmetic.
- `arm64`: builds an additional assembly optimized implementation targeting the ARM 64 architecture (and later). The optimizations are applied to the finite field arithmetic.

### MIKE_PRIME_CHOICE

Specifies if fast or conservative primes should be chosen to reach NIST security levels I, III, V. The currently supported values are:
- `fast`: for primes `p308_633`, `p474_593` and `p628_317` for NIST levels I, III and V respectively.
- `conservative`: for primes `p374_117`, `p566_77` and `p758_41` for NIST levels I, III and V respectively.

where `pf_c` denotes a prime of the form $c\cdot 2^f-1$. 

### CMAKE_BUILD_TYPE

Can be used to specify special build types. The options are:

- `Release`: Builds with optimizations enabled and assertions disabled.
- `Debug`: Builds with debug symbols.
- `ASAN`: Builds with AddressSanitizer memory error detector.
- `UBSAN`: Builds with UndefinedBehaviorSanitizer for undefined behavior detection.

The default build type uses the flags `-O2 -Wstrict-prototypes -Wno-error=strict-prototypes -fvisibility=hidden -Wno-error=implicit-function-declaration -Wno-error=attributes`. (Notice that assertions remain enabled in this configuration, which harms performance.)


## Pre-computation

Finite field arithmetic in pure C (`ref`) and in ASM (`broadwell`, `arm64`), constants from `src/precomp` and MIKE namespace can be generated from the `scripts` directory. We refer to the file `script/README.md` for detailed instructions. 

For each prime parameter, a synthetic description of the parameter set is provided in `src/precomp/ref/.../mike_parameters.txt`.

## Test

In the build directory, run `make test` or `ctest`.

The test harness consists of the following units:

- Self-tests: `MIKE_test_scheme_<pe_c>` - runs random self-tests (key generation & key exchange).
- Sub-library specific unit-tests.

Note that, `ctest` has a default timeout of 1500s, which is applied to all tests. To override the default timeout, run
`ctest --timeout <seconds>`.


## Benchmarks

A benchmarking suite is built and can be executed with the following command:
```
test/benchmark_<pe_c> [--iterations=<iterations>]
```
where `<pe_c>` specifies the MIKE parameter set and `<iterations>` is the
number of iterations used for benchmarking; if the `--iterations` option is omitted, a default of 100 iterations is used.

The benchmarks profile the key generation and key exchange functions. The results are reported in CPU cycles if available on the host platform, and timing in nanoseconds otherwise.

## Project Structure

The source code consists of a number of sub-libraries used to implement the
final MIKE library:
- `common`: common code for hash function, seed expansion, PRNG, memory handling.
- `gf`: GF(p^2) and GF(p) arithmetic.
- `ec`: elliptic curves, dim1 isogenies and pairings. Everything that is purely defined using GF(p^2) arithmetic.
- `precomp`: constants and precomputed values.
- `theta`: dimension 2 and 4 abelian varieties and their generic isogenies in the theta model. Everything that is purely defined using GF(p) arithmetic.
- `weil`: gluing, invariant, mike's basis generation. Everything requiring both GF(p) and GF(p^2) and highly specific to the mike key exchange.  
- `nike`: code for the key generation and exchange protocols.


## Acknowledgements

- The general infrastructure of this project is inherited from [SQIsign](https://github.com/SQIsign/the-sqisign).

- The ec module implementation of the (`src/ec`) is a modification from the ec module of [qt-Pegasis-Fp](https://github.com/Pierrick-Dartois/qt-pegasis-Fp/tree/main) and the isogeny code and normalization of the latest version of [SQIsign](https://github.com/SQIsign/the-sqisign). 


- The reference implementation for finite field arithemtic (i.e., `src/gf/ref`)
was generated using [modarith](https://github.com/mcarrickscott/modarith) by
Michael Scott.


## License

MIKE is licensed under Apache-2.0. See [LICENSE](LICENSE) and [NOTICE](NOTICE).

Third party code is used in some files:

- `src/common/aes_c.c`; MIT: "Copyright (c) 2016 Thomas Pornin <pornin@bolet.org>"
- `src/common/fips202.c`: CC0: Copyright (c) 2023, the PQClean team
- `src/common/randombytes_system.c`: MIT: Copyright (c) 2017 Daan Sprenkels <hello@dsprenkels.com>
- `src/common/broadwell/{aes_ni.c, vaes256_key_expansion.S}`: Apache-2.0: Copyright 2019 Amazon.com, Inc.
- `src/common/broadwell/ctr_drbg.c`: ISC: Copyright (c) 2017, Google Inc.
