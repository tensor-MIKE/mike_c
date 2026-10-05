#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <locale.h>
#include <time.h>

#include <nike.h>

#include <tools.h>
#include <rng.h>
#include <bench.h>
#include <bench_test_arguments.h>

#define STRINGIFY2(x) #x
#define STRINGIFY(x) STRINGIFY2(x)

static __inline__ uint64_t
rdtsc(void)
{
    return (uint64_t)cpucycles();
}

int
bench_mike(uint64_t bench)
{
    int res = 0;
    setlocale(LC_NUMERIC, "");
    uint64_t t0, t1;
    clock_t t;
    double ms;

    public_key_t *pksA = calloc(bench, sizeof(public_key_t));
    public_key_t *pksB = calloc(bench, sizeof(public_key_t));

    secret_key_t *sksA = calloc(bench, sizeof(secret_key_t));
    secret_key_t *sksB = calloc(bench, sizeof(secret_key_t));

    shared_t *sharesA = calloc(bench, sizeof(shared_t));
    shared_t *sharesB = calloc(bench, sizeof(shared_t));

    printf("\n\nBenchmarking MIKE for " STRINGIFY(MIKE_VARIANT) ":\n\n");
    t = tic();
    t0 = rdtsc();
    for (uint64_t i = 0; i < bench; i++) {
        if (!(protocols_keygen(&pksA[i], &sksA[i]) && protocols_keygen(&pksB[i], &sksB[i]))) {
            res = 1;
            printf("keygen failed\n");
            goto fin;
        }
    }
    t1 = rdtsc();
    ms = (1000. * (double)(clock() - t) / CLOCKS_PER_SEC);
    printf("Average keygen time [%.2f ms]\n", (double)(ms / (2 * bench)));
    printf("\x1b[34mAvg keygen: %'" PRIu64 " cycles\x1b[0m\n", (t1 - t0) / (2 * bench));

    t = tic();
    t0 = rdtsc();
    for (uint64_t i = 0; i < bench; ++i) {
        protocols_exchange(&sharesA[i], &pksB[i], &sksA[i]);
        protocols_exchange(&sharesB[i], &pksA[i], &sksB[i]);
    }
    t1 = rdtsc();
    ms = (1000. * (double)(clock() - t) / CLOCKS_PER_SEC);
    printf("Average key exchange time [%.2f ms]\n", (double)(ms / (2 * bench)));
    printf("\x1b[34mAvg key exchange: %'" PRIu64 " cycles\x1b[0m\n", (t1 - t0) / (2 * bench));

    for (uint64_t i = 0; i < bench; ++i) {
        for (uint8_t j = 0; j < 4; j++) {
            if (!fp_is_equal(&(sharesA[i].invariant[j]), &(sharesB[i].invariant[j]))) {
                res = 1;
                printf("mike does not compute the same key\n");
                goto fin;
            }
        }
    }

fin:;

    free(pksA);
    free(pksB);
    free(sksA);
    free(sksB);
    free(sharesA);
    free(sharesB);
    return (res);
}

// run all tests in module
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

    int res = bench_mike(iterations);

    return res;
}
