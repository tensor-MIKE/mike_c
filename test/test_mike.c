// SPDX-License-Identifier: Apache-2.0

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <rng.h>
#include <ex.h>
#include <bench_test_arguments.h>
#ifdef TARGET_BIG_ENDIAN
#include <tutil.h>
#endif

static void
print_hex(const unsigned char *hex, int len)
{
    for (int i = 0; i < len; ++i) {
        printf("%02x", hex[i]);
    }
    printf("\n");
}

// Tests that crypto_keygen() and crypto_exchange() reject a wrong pk length, and that they do so before
// reading the pk buffer at all.
//
// Two details ensure that these tests notice if those length checks ever go missing:
//
//  - Every buffer is allocated to exactly the length that is then passed to the function, so a sanitizer flags any read
//    past that length. Writing the tests with one full-size buffer, called with a smaller length, would read from a
//    memory region that is still validly allocated, where nothing would complain.
//
//  - key exchange works
static int
test_invalid_lengths(int interations)
{

    unsigned char *sk_good = calloc(1, SECRETKEY_BYTES);
    unsigned char *pk_good = calloc(1, PUBLICKEY_BYTES);

    unsigned char *pk_bad = calloc(1, PUBLICKEY_BYTES + 1);
    unsigned char *shared_key = calloc(1, 32);

    fp2_t A;
    int res = -1;
    int ret = 0;

    // sample good keys
    mike_keypair(pk_good, sk_good);

    for (int g = 0; g < interations; g++) {
        // test with curve defined over Fp.
        randombytes(pk_bad, FP_ENCODED_BYTES);
        fp_decode_reduce(&A.re, pk_bad, FP_ENCODED_BYTES);
        fp_set_zero(&A.im);
        fp2_encode(pk_bad, &A);

        res = mike_exchange(shared_key, pk_bad, sk_good);

        if (res) {
            printf("mike exchange accepts curves defined over Fp \n");
            ret = 1;
            goto fin;
        }
        res = -1;

        // test with random Fp2 A ( = regular curves with proba 1/p)
        randombytes(pk_bad, FP2_ENCODED_BYTES);
        fp_decode_reduce(&A.re, pk_bad, FP_ENCODED_BYTES);
        fp_decode_reduce(&A.im, pk_bad + FP_ENCODED_BYTES, FP_ENCODED_BYTES);
        fp2_encode(pk_bad, &A);
        res = mike_exchange(shared_key, pk_bad, sk_good);

        if (res) {
            printf("mike exchange worked on a regular curve \n");
            ret = 1;
            goto fin;
        }
        res = -1;

        // test with encoding that is too small
        randombytes(pk_bad, FP2_ENCODED_BYTES - 1);
        res = mike_exchange(shared_key, pk_bad, sk_good);

        if (res) {
            printf("mike exchange worked on bad inputs \n");
            ret = 1;
            goto fin;
        }
        res = -1;

        // test with encoding that is too big
        randombytes(pk_bad, FP2_ENCODED_BYTES + 1);
        res = mike_exchange(shared_key, pk_bad, sk_good);

        if (res) {
            printf("mike exchange worked on bad inputs \n");
            ret = 1;
            goto fin;
        }
        res = -1;

        // test zero vector
        for (int i = 0; i < FP2_ENCODED_BYTES; i++) {
            pk_bad[i] = 0;
        }
        res = mike_exchange(shared_key, pk_bad, sk_good);

        if (res) {
            printf("mike exchange worked on bad inputs \n");
            ret = 1;
            goto fin;
        }
        res = -1;
    }
fin:;
    free(sk_good);
    free(pk_good);
    free(pk_bad);
    free(shared_key);
    if (!res)
        printf("Reject wrong public keys............................... PASSED\n");
    return ret;
}

static int
test_mike(int interations)
{

    unsigned char *skA = calloc(1, SECRETKEY_BYTES);
    unsigned char *skB = calloc(1, SECRETKEY_BYTES);

    unsigned char *pkA = calloc(1, PUBLICKEY_BYTES);
    unsigned char *pkB = calloc(1, PUBLICKEY_BYTES);

    unsigned char *sharedA = calloc(1, 32);
    unsigned char *sharedB = calloc(1, 32);

    int resA = 0, resB = 0;
    int res = 0;

    for (int g = 0; g < interations; g++) {

        mike_keypair(pkA, skA);
        mike_keypair(pkB, skB);

        resA = mike_exchange(sharedA, pkB, skA);
        resB = mike_exchange(sharedB, pkA, skB);

        if (!resA) {
            printf("Pb during Key Agreement Alice \n");
            return 1;
        }
        if (!resB) {
            printf("Pb during Key Agreement Bob \n");
            return 1;
        }

        for (uint8_t i = 0; i < 32; i++) {
            if (sharedA[i] != sharedB[i]) {
                res = 1;
                goto end;
            }
        }
    }

end:;
    free(skA);
    free(skB);
    free(pkA);
    free(pkB);
    free(sharedA);
    free(sharedB);

    if (res)
        printf("Alice and Bob did not compute the same shared secret \n");
    else
        printf("Same shared key........................................ PASSED\n");
    return res;
}

int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int help = 0;
    int seed_set = 0;
    int interations = 50;
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
    }

    if (help) {
        printf("Usage: %s [--seed=<seed>]\n", argv[0]);
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

    res = test_invalid_lengths(interations);
    res |= test_mike(interations);

    if (res) {
        printf("Tests failed!\n");
    } else {
        printf("All MIKE tests passed.\n");
    }
}
