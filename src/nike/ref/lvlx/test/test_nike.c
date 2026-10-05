#include <nike.h>
#include <rng.h>
#include <constants.h>
#include <bench_test_arguments.h>


// Make n random-ish field elements (for tests only!).
static void
test_nike_fp_random_test(fp_t *a)
{
    uint8_t tmp[FP_ENCODED_BYTES];

    randombytes(tmp, sizeof(tmp));

    fp_decode_reduce(a, tmp, sizeof(tmp));
}

static void
nike_test_fp2_random_test(fp2_t *a)
{
    test_nike_fp_random_test(&(a->re));
    test_nike_fp_random_test(&(a->im));
}

int
nike_test_sk_encoding(int iterations)
{

    int res = 0;
    int index = 0;
    secret_key_t sk;
    digit_t mask, val;
    unsigned char xbytes[SECRETKEY_BYTES];
    unsigned char ybytes[SECRETKEY_BYTES];

    for (int i = 0; i < iterations; i++) {
        index = i % SECRETKEY_BYTES;
        for (int k = 0; k < SECRETKEY_BYTES; k++) {
            xbytes[k] = 0;
        }
        randombytes(&xbytes[index], 1);
        secret_key_from_bytes(&sk, xbytes);
        for (int k = 0; k < SECRETKEY_WORDS; k++) {
            if (k != (int)(index / sizeof(digit_t)))
                res = res | (sk.x[k] != 0);
        }
        mask = 0xff;
        mask = mask << (8 * (index % sizeof(digit_t)));
        val = sk.x[index / (sizeof(digit_t))] & mask;
        val = val >> (8 * (index % sizeof(digit_t)));

        res = res | (((uint8_t)val) != xbytes[index]);

        if (res) {
            printf("NIKE sk decoding failed\n");
            break;
        }

        secret_key_to_bytes(ybytes, &sk);

        for (int k = 0; k < SECRETKEY_BYTES; k++) {
            res = res | (xbytes[k] != ybytes[k]);
        }

        if (res) {
            printf("NIKE sk encoding failed\n");
            break;
        }
    }

    // try of fully random sk keys.
    for (int i = 0; i < iterations; i++) {
        index = i % SECRETKEY_BYTES;
        for (int k = 0; k < SECRETKEY_BYTES; k++) {
            xbytes[k] = 0;
        }
        randombytes(xbytes, SECRETKEY_BYTES);
        secret_key_from_bytes(&sk, xbytes);
        secret_key_to_bytes(ybytes, &sk);

        for (int k = 0; k < SECRETKEY_BYTES; k++) {
            res = res | (xbytes[k] != ybytes[k]);
        }

        if (res) {
            printf("NIKE sk encoding failed\n");
            break;
        }
    }

    if (!res)
        printf("NIKE sk encoding test................................................ PASSED\n");
    return (res);
}

int
nike_test_pk_encoding(int iterations)
{

    int res = 0;
    public_key_t pk, pk2;
    unsigned char xbytes[PUBLICKEY_BYTES];

    // try of fully random sk keys.
    for (int i = 0; i < iterations; i++) {
        for (int k = 0; k < PUBLICKEY_BYTES; k++) {
            xbytes[k] = 0;
        }
        randombytes(xbytes, PUBLICKEY_BYTES);
        nike_test_fp2_random_test(&pk.curveA); 
        public_key_to_bytes(xbytes, &pk);
        public_key_from_bytes(&pk2, xbytes);

        res = !fp2_is_equal(&pk.curveA, &pk2.curveA); 

        if (res) {
            printf("NIKE pk encoding failed\n");
            break;
        }
    }

    if (!res)
        printf("NIKE pk encoding test................................................ PASSED\n");
    return (res);
}

int
nike_test_keygen(int iterations)
{
    int res = 0;
    secret_key_t sk;
    public_key_t pk;
    ec_jac_point_t P, Q, pc, qc;
    for (int iter = 0; iter < iterations; iter++) {
        res = res | !protocols_keygen(&pk, &sk);
        if (res) {
            printf("NIKE keygen test failed: returned error\n");
            break;
        }
        res = res | !((sk.x[0] & 7) == 3);
        if (res) {
            printf("NIKE keygen test failed: sk not 3 mod 8\n");
            break;
        }
        int nonzeros = (SECRETKEY_BITS % (sizeof(digit_t) * 8));
        digit_t mask = -1;
        mask = (mask << (nonzeros));
        res = res | !((sk.x[SECRETKEY_WORDS - 1] & mask) == 0);
        if (res) {
            printf("NIKE keygen test failed: sk too long\n");
            break;
        }
        ec_curve_t curve;
        ec_curve_init_from_A(&curve, &pk.curveA);
        ec_jac_basis_2f(&P, &Q, &curve, TORSION_EVEN_POWER);
        ec_jac_dbl_iter(&pc, &P, TORSION_EVEN_POWER - 1, &curve, 1);
        ec_jac_dbl_iter(&qc, &Q, TORSION_EVEN_POWER - 1, &curve, 1);
        res = res | ec_jac_is_zero(&pc);
        res = res | ec_jac_is_zero(&qc);
        ec_DBL(&P, &pc, &curve);
        ec_DBL(&Q, &qc, &curve);
        res = res | !ec_jac_is_zero(&P);
        res = res | !ec_jac_is_zero(&Q);
        if (res) {
            printf("NIKE keygen test failed: pk not supersingular\n");
            break;
        }
    }
    if (!res)
        printf("NIKE keygen test..................................................... PASSED\n");
    return (res);
}

int
nike_test_nike(int iterations)
{
    secret_key_t sk_A, sk_B;
    public_key_t pk_A, pk_B;
    shared_t share_A, share_B;
    int resA, resB;
    for (int iter = 0; iter < iterations; iter++) {
        protocols_keygen(&pk_A, &sk_A);
        protocols_keygen(&pk_B, &sk_B);

        resA = protocols_exchange(&share_A, &pk_B, &sk_A);
        resB = protocols_exchange(&share_B, &pk_A, &sk_B);

        if (!resA) {
            printf("Pb during Key Agreement Alice \n");
            return 1;
        }
        if (!resB) {
            printf("Pb during Key Agreement Bob \n");
            return 1;
        }

        for (uint8_t i = 0; i < 4; i++) {
            if (!fp_is_equal(&share_A.invariant[i], &share_B.invariant[i])) {
                printf("Alice and Bob do not compute the same shared secret \n");
                return 1;
            }
        }
    }
    printf("NIKE same invariant.................................................. PASSED\n");
    return 0;
}

int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int iterations = MIKE_TEST_REPS;
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

    printf("--------------------------------------------------------------------------------\n\n");
    printf("Testing NIKE functions and encodings:\n\n");
    print_seed(seed);

    // #if defined(TARGET_BIG_ENDIAN)
    // for (int i = 0; i < 12; i++) {
    // seed[i] = BSWAP32(seed[i]);
    //}
    // #endif

    randombytes_init((unsigned char *)seed, NULL, 256);
    res |= nike_test_sk_encoding(iterations);
    res |= nike_test_pk_encoding(iterations);
    res |= nike_test_keygen(iterations);
    res |= nike_test_nike(iterations);

    if (res) {
        printf("Tests failed!\n");
    } else {
        printf("All NIKE tests passed.\n");
    }

    return res;
}
