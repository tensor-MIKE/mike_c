#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <rng.h>
#include <encoded_sizes.h>
#include <inttypes.h>
#include <bench_test_arguments.h>

#include <theta.h>

void
fp_random_test(fp_t *a)
{
    uint8_t tmp[FP_ENCODED_BYTES];

    randombytes(tmp, FP_ENCODED_BYTES);

    fp_decode_reduce(a, tmp, FP_ENCODED_BYTES);
}

int
test_hadamard(int iterations)
{
    fp_t tmp, tmp2, tmp3, tmp4;
    theta_point_dim2_t before_2, after_2;
    theta_point_dim4_t before_4, after_4;
    theta_point_weil_t before_weil, after_weil;

    for (int g = 0; g < iterations; g++) {

        // TEST DIM 2

        for (uint8_t i = 0; i < 4; i++) {
            fp_random_test(&before_2[i]);
        }

        theta_hadamard_dim2(&after_2, &before_2);
        theta_hadamard_dim2(&after_2, &after_2);

        fp_copy(&tmp, &before_2[0]);
        fp_copy(&tmp3, &after_2[0]);
        fp_inv(&tmp);
        fp_inv(&tmp3);

        for (uint8_t i = 0; i < 4; i++) {
            fp_mul(&tmp2, &tmp, &before_2[i]);
            fp_mul(&tmp4, &tmp3, &after_2[i]);
            if (!fp_is_equal(&tmp2, &tmp4)) {
                printf("dim2 hadamard not self inverse");
                return 1;
            }
        }

        // TEST DIM 4

        for (uint8_t i = 0; i < 16; i++) {
            fp_random_test(&before_4[i]);
        }

        theta_hadamard_dim4(&after_4, &before_4);
        theta_hadamard_dim4(&after_4, &after_4);

        fp_copy(&tmp, &before_4[0]);
        fp_copy(&tmp3, &after_4[0]);
        fp_inv(&tmp);
        fp_inv(&tmp3);

        for (uint8_t i = 0; i < 16; i++) {
            fp_mul(&tmp2, &tmp, &before_4[i]);
            fp_mul(&tmp4, &tmp3, &after_4[i]);
            if (!fp_is_equal(&tmp2, &tmp4)) {
                printf("dim4 hadamard not self inverse");
                return 1;
            }
        }

        // TEST WEIL

        for (uint8_t i = 0; i < 10; i++) {
            fp_random_test(&before_weil[i]);
        }

        theta_hadamard_weil(&after_weil, &before_weil);
        theta_hadamard_weil(&after_weil, &after_weil);

        fp_copy(&tmp, &before_weil[0]);
        fp_copy(&tmp3, &after_weil[0]);
        fp_inv(&tmp);
        fp_inv(&tmp3);

        for (uint8_t i = 0; i < 10; i++) {
            fp_mul(&tmp2, &tmp, &before_weil[i]);
            fp_mul(&tmp4, &tmp3, &after_weil[i]);
            if (!fp_is_equal(&tmp2, &tmp4)) {
                printf("weil-hadamard not self inverse");
                return 1;
            }
        }
    }

    printf("Hadamard ......................................................PASSED\n");

    return 0;
}

int
test_proj_inverse(int iterations)
{
    fp_t tmp, tmp2, tmp3, tmp4;
    theta_point_dim2_t before_2, after_2;
    theta_point_dim4_t before_4, after_4;
    theta_point_weil_t before_weil, after_weil;

    for (int g = 0; g < iterations; g++) {

        // TEST DIM 2

        for (uint8_t i = 0; i < 4; i++) {
            fp_random_test(&before_2[i]);
        }

        theta_invert_dim2(&after_2, &before_2);
        theta_invert_dim2(&after_2, &after_2);

        fp_copy(&tmp, &before_2[0]);
        fp_copy(&tmp3, &after_2[0]);
        fp_inv(&tmp);
        fp_inv(&tmp3);

        for (uint8_t i = 0; i < 4; i++) {
            fp_mul(&tmp2, &tmp, &before_2[i]);
            fp_mul(&tmp4, &tmp3, &after_2[i]);
            if (!fp_is_equal(&tmp2, &tmp4)) {
                printf("dim2 inverse not self inverse\n");
                return 1;
            }
        }

        // TEST DIM 4

        for (uint8_t i = 0; i < 16; i++) {
            fp_random_test(&before_4[i]);
        }

        theta_invert_theta_dim4(&after_4, &before_4);
        theta_invert_theta_dim4(&after_4, &after_4);

        fp_copy(&tmp, &before_4[0]);
        fp_copy(&tmp3, &after_4[0]);
        fp_inv(&tmp);
        fp_inv(&tmp3);

        for (uint8_t i = 0; i < 16; i++) {
            fp_mul(&tmp2, &tmp, &before_4[i]);
            fp_mul(&tmp4, &tmp3, &after_4[i]);
            if (!fp_is_equal(&tmp2, &tmp4)) {
                printf("dim2 inverse not self inverse\n");
                return 1;
            }
        }

        // TEST WEIL

        for (uint8_t i = 0; i < 10; i++) {
            fp_random_test(&before_weil[i]);
        }

        theta_invert_weil(&after_weil, &before_weil);
        theta_invert_weil(&after_weil, &after_weil);

        fp_copy(&tmp, &before_weil[0]);
        fp_copy(&tmp3, &after_weil[0]);
        fp_inv(&tmp);
        fp_inv(&tmp3);

        for (uint8_t i = 0; i < 10; i++) {
            fp_mul(&tmp2, &tmp, &before_weil[i]);
            fp_mul(&tmp4, &tmp3, &after_weil[i]);
            if (!fp_is_equal(&tmp2, &tmp4)) {
                printf("Weil inverse not self inverse\n");
                return 1;
            }
        }
    }

    printf("Inverse .......................................................PASSED\n");

    return 0;
}

int
test_weil_compress(int iterations)
{
    fp_t tmp, tmp2, tmp3, tmp4;
    theta_point_dim4_t after_theta, before_theta;
    theta_point_dim4_t point;
    theta_point_weil_t before_weil, after_weil;

    for (int g = 0; g < iterations; g++) {

        // Test conversion between dim 2 and dim 4
        for (uint8_t i = 0; i < 10; i++) {
            fp_random_test(&before_weil[i]);
        }

        theta_weil_to_dim4_point(&before_theta, &before_weil);
        theta_dim4_to_weil(&after_weil, &before_theta);

        fp_copy(&tmp, &before_weil[0]);
        fp_copy(&tmp3, &after_weil[0]);
        fp_inv(&tmp);
        fp_inv(&tmp3);

        for (uint8_t i = 0; i < 10; i++) {
            fp_mul(&tmp2, &tmp, &before_weil[i]);
            fp_mul(&tmp4, &tmp3, &after_weil[i]);
            if (!fp_is_equal(&tmp2, &tmp4)) {
                printf("Weil inverse not self inverse\n");
                return 1;
            }
        }

        // TEST weil dot product

        for (uint8_t i = 0; i < 16; i++) {
            fp_random_test(&point[i]);
        }

        theta_dot_prod_weil(&after_theta, &point, &before_weil);
        theta_weil_to_dim4_point(&before_theta, &before_weil);
        theta_dot_prod_dim4(&point, &point, &before_theta);

        fp_copy(&tmp, &after_theta[0]);
        fp_copy(&tmp3, &point[0]);
        fp_inv(&tmp);
        fp_inv(&tmp3);

        for (uint8_t i = 0; i < 16; i++) {
            fp_mul(&tmp2, &tmp, &after_theta[i]);
            fp_mul(&tmp4, &tmp3, &point[i]);
            if (!fp_is_equal(&tmp2, &tmp4)) {
                printf("Weil dot product not correct\n");
                return 1;
            }
        }
    }

    printf("Weil compression ..............................................PASSED\n");

    return 0;
}

int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int iterations = 100 * MIKE_TEST_REPS;
    int help = 0;
    int seed_set = 0;
    int res = 0;

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
        printf("Where <iterations> is the number of iterations used for testing; if not "
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

    res |= test_hadamard(iterations);     // Replace by unit tests (output 0 on success) here
    res |= test_proj_inverse(iterations); // Replace by unit tests (output 0 on success) here
    res |= test_weil_compress(iterations);

    if (res) {
        printf("Tests failed!\n");
    } else {
        printf("All Theta module tests passed.\n");
    }

    return res;
}
