// SPDX-License-Identifier: Apache-2.0

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>

#include <ex.h>
#include <rng.h>
#include <bench.h>
#include <bench_test_arguments.h>
#if defined(TARGET_BIG_ENDIAN)
#include <tutil.h>
#endif

void
bench(size_t runs)
{
    const size_t sharedkey_len = 32;

    unsigned char *pk = calloc(runs, PUBLICKEY_BYTES);
    unsigned char *sk = calloc(runs, SECRETKEY_BYTES);
    unsigned char *sharedkey = calloc(runs, sharedkey_len);

    printf(" MIKE (%zu iterations)\n", runs);


    BENCH_CODE_1(runs);
    mike_keypair(pk, sk);
    BENCH_CODE_2("KeyGen");

    BENCH_CODE_1(runs);
    mike_exchange(sharedkey, pk, sk);
    BENCH_CODE_2("Exchange");

    free(pk);
    free(sk);
    free(sharedkey);

}

int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int iterations = MIKE_TEST_REPS;
    int help = 0;
    int seed_set = 0;

#ifndef NDEBUG
    fprintf(stderr,
            "\x1b[31mIt looks like MIKE was compiled with assertions enabled.\n"
            "This will severely impact performance measurements.\x1b[0m\n");
#endif

    for (int i = 1; i < argc; i++) {
        if (!help && strcmp(argv[i], "--help") == 0) {
            help = 1;
            continue;
        }

        if (!seed_set && !parse_seed(argv[i], seed)) {
            seed_set = 1;
            continue;
        }

        if (sscanf(argv[i], "--iterations=%d", &iterations) == 1) {
            continue;
        }
    }

    if (help || iterations <= 0) {
        printf("Usage: %s [--iterations=<iterations>] [--seed=<seed>]\n", argv[0]);
        printf("Where <iterations> is the number of iterations used for benchmarking; if not "
               "present, uses the default: %d)\n",
               iterations);
        printf("Where <seed> is the random seed to be used; if not present, a random seed is "
               "generated\n");
        return 1;
    }

    if (!seed_set) {
        randombytes_select((unsigned char *)seed, sizeof(seed));
    }

    print_seed(seed);

#if defined(TARGET_BIG_ENDIAN)
    for (int i = 0; i < 12; i++) {
        seed[i] = BSWAP32(seed[i]);
    }
#endif

    randombytes_init((unsigned char *)seed, NULL, 256);
    cpucycles_init();

    bench(iterations);

    return 0;
}
