#include <ec.h>
#include <ec_isog.h>
#include <stdio.h>
#include <bench_test_arguments.h>
#include <rng.h>
#include "test_extras.h"

// test helper functions
int
ec_test_A24_conversion(int iterations)
{
    int res = 0;
    fp2_t rnd;
    ec_xz_point_t A24;
    ec_curve_t curve, test;
    ec_curve_init_precomputed(&curve);
    fp2_set_small(&curve.A, 6);
    ec_compute_A24(&curve);

    for (int i = 0; i < iterations; i++) {
        ec_curve_to_A24(&A24, &curve);
        // normalize
        if (fp2_is_zero(&A24.z)) {
            fp2_set_one(&A24.x);
        } else {
            fp2_inv(&A24.z);
            fp2_mul(&A24.x, &A24.x, &A24.z);
            fp2_set_one(&A24.z);
        }
        // test
        ec_A24_to_curve_A(&test.A, &A24);
        ec_curve_init_from_A(&test, &test.A);
        // check here if both curves are the same
        res = res || !fp2_is_equal(&test.A, &test.A);
        if (res) {
            printf("Curve to A24 conversion tests failed A equality at iteration %d\n", i);
            break;
        }
        // sample next curve
        ec_fp2_random_test(&(rnd));
        fp2_mul(&(curve.A), &(curve.A), &(rnd));
        fp2_inv(&rnd);
        fp2_mul(&curve.A, &curve.A, &rnd);
    }
    if (!res) {
        printf("Curve to A24 conversion tests........................................ PASSED\n");
    } else {
        printf("Curve to A24 conversion tests........................................ FAILED\n");
    }
    return (res);
}

int
ec_test_A24_xDBL(int iterations)
{
    int res = 0;
    fp2_t rnd;
    ec_xz_point_t A24, P, Q, T;
    ec_curve_t curve;
    ec_curve_init_precomputed(&curve);
    fp2_set_small(&curve.A, 6);
    ec_compute_A24(&curve);

    for (int i = 0; i < iterations; i++) {
        ec_xz_random_test(&P, &curve);
        ec_curve_to_A24(&A24, &curve);
        ec_xDBL_A24(&Q, &P, &A24);
        ec_xDBL(&T, &P, &curve);
        res = res || !ec_xz_is_equal(&T, &Q);
        if (res) {
            printf("A24 point doubling test failed in iteration %d\n", i);
            break;
        }
        // sample next curve
        ec_fp2_random_test(&(rnd));
        fp2_mul(&(curve.A), &(curve.A), &(rnd));
        fp2_inv(&rnd);
        fp2_mul(&curve.A, &curve.A, &rnd);
    }
    if (!res) {
        printf("A24 point doubling................................................... PASSED\n");
    } else {
        printf("A24 point doubling................................................... FAILED\n");
    }
    return (res);
}

int
ec_xz_order_even(const ec_xz_point_t *P, ec_curve_t *curve)
{
    ec_xz_point_t Q, R;
    ec_compute_A24(curve);
    fp2_copy(&Q.x, &P->x);
    fp2_copy(&Q.z, &P->z);
    if (fp2_is_zero(&Q.x) && !fp2_is_zero(&Q.z))
        return (1);
    ec_xMUL(&Q, &R, &Q, TORSION_ODD, ec_n_bits(TORSION_ODD, NWORDS_ORDER), curve);
    int order = 1;
    if (ec_xz_is_zero(&Q)) {
        return (0);
    }
    while (!ec_xz_is_two_torsion(&Q, curve)) {
        ec_xDBL(&Q, &Q, curve);
        order = order + 1;
    }
    return (order);
}

int
ec_test_generated_basis()
{
    int res = 0;
    int f = 6;
    ec_xz_point_t p, q;
    ec_jac_point_t jp, jq;
    ec_curve_t curve;
    fp2_t A;
    fp_set_zero(&A.im);
    fp_set_small(&A.re, 6);
    ec_curve_init_from_A(&curve, &A);
    ec_jac_basis_2f(&jp, &jq, &curve, f);
    ec_jac_to_xz(&p, &jp);
    ec_jac_to_xz(&q, &jq);
    int Porder = ec_xz_order_even(&p, &curve);
    int Qorder = ec_xz_order_even(&q, &curve);
    res = res || !(Porder == Qorder);
    res = res || !(Porder == f);
    for (int i = 1; i < Qorder; i++) {
        ec_xDBL(&q, &q, &curve);
        ec_xDBL(&p, &p, &curve);
    }
    res = res || (ec_xz_is_equal(&p, &q));
    ec_xDBL(&q, &q, &curve);
    ec_xDBL(&p, &p, &curve);
    res = res || !ec_xz_is_zero(&p);
    res = res || !ec_xz_is_zero(&q);
    if (!res) {
        printf("Generated basis ..................................................... PASSED\n");
    } else {
        printf("Generated basis test failed\n");
    }

    return (res);
}

int
ec_test_precomputed_basis()
{
    int res = 0;
    ec_xz_point_t p, q;
    ec_curve_t curve;
    fp2_t A;
    fp_set_zero(&A.im);
    fp_copy(&A.re, &A0);
    ec_curve_init_from_A(&curve, &A);
    fp2_copy(&p.x, &BASIS_EM_PX);
    fp2_copy(&q.x, &BASIS_EM_QX);
    fp2_set_one(&p.z);
    fp2_set_one(&q.z);
    int Porder = ec_xz_order_even(&p, &curve);
    int Qorder = ec_xz_order_even(&q, &curve);
    res = res || !(Porder == Qorder);
    res = res || !(Porder == TORSION_EVEN_POWER);
    for (int i = 1; i < Qorder; i++) {
        ec_xDBL(&q, &q, &curve);
        ec_xDBL(&p, &p, &curve);
    }
    res = res || (ec_xz_is_equal(&p, &q));
    res = res || !(fp2_is_zero(&q.x));
    ec_xDBL(&q, &q, &curve);
    ec_xDBL(&p, &p, &curve);
    res = res || !ec_xz_is_zero(&p);
    res = res || !ec_xz_is_zero(&q);
    if (!res) {
        printf("Precomputed basis ................................................... PASSED\n");
    } else {
        printf("Precomputed basis ................................................... FAILED\n");
    }
    return (res);
}

int
ec_test_chain(int iterations)
{
    int res = 0;
    int length;
    int found = 0;
    digit_t diff[NWORDS_ORDER];
    fp2_t imageA;
    ec_xz_point_t P, Q, R, S;
    ec_xz_point_t ps[5];
    ec_xz_point_t psc[5];
    ec_jac_point_t psj[5];
    ec_jac_point_t Pj;
    ec_curve_t curve, save_curve;
    ec_isog_even_t phi;

    ec_curve_init_precomputed(&curve);
    fp2_set_small(&curve.A, 6);
    for (int iter = 0; iter < iterations; iter++) {
        found = 0;
        length = 2 * iter + 10;
        for (long unsigned int k = 0; k < NWORDS_ORDER; k++) {
            diff[k] = 0;
            if (k == ((TORSION_EVEN_POWER - length) / (sizeof(digit_t) * 8))) {
                int shift = TORSION_EVEN_POWER - length -
                            (TORSION_EVEN_POWER - length) / (sizeof(digit_t) * 8) * (sizeof(digit_t) * 8);
                diff[k] = ((digit_t)1) << shift;
            }
        }

        ec_compute_A24(&curve);
        while (!found) {
            // sample random point
            ec_xz_random_normalized_test(&P, &curve);
            assert(ec_x_is_on_curve(&P.x, &curve));
            // multiply with odd*2^(two_torsion_power-length)
            ec_xMUL(&R, &Q, &P, TORSION_ODD, ec_n_bits(TORSION_ODD, NWORDS_ORDER), &curve);
            ec_xMUL(&P, &Q, &R, diff, ec_n_bits(diff, NWORDS_ORDER), &curve);
            // alternative point sampling using a basis (commented)
            //  ec_basis_2f(&P, &Q, &curve, TORSION_EVEN_POWER);
            //  printf("order P basis %d\n",ec_xz_order_even(&P,&curve));
            //  ec_xMUL(&R, &Q, &P, diff, ec_n_bits(diff, NWORDS_ORDER), &curve);
            //  ec_copy_xz_point(&P,&R);
            //  check if order is right
            if (ec_xz_order_even(&P, &curve) != length)
                continue;
            // check if above (0,0)
            ec_copy_xz_point(&Q, &P);
            assert(ec_xz_order_even(&Q, &curve) == length);
            for (int k = length; k > 1; k--) {
                ec_xDBL(&Q, &Q, &curve);
                ec_xz_normalize(&Q);
                assert(ec_x_is_on_curve(&Q.x, &curve));
            }
            assert(ec_xz_order_even(&Q, &curve) == 1);
            if (fp2_is_zero(&Q.x))
                continue;
            found = 1;
        }
        assert(ec_xz_order_even(&P, &curve) == length);

        // compute 4P  (used to test if it gets mapped to 0)
        ec_xDBL(&Q, &P, &curve);
        ec_xDBL(&ps[0], &Q, &curve);
        if (!fp2_is_zero(&ps[0].z)) {
            fp2_inv(&ps[0].z);
            fp2_mul(&ps[0].x, &ps[0].x, &ps[0].z);
            fp2_set_one(&ps[0].z);
            assert(ec_x_is_on_curve(&ps[0].x, &curve));
            assert(ec_xz_order_even(&ps[0], &curve) == length - 2);
        }
        // save data defining isogeny for 2nd isogeny computation (with normalization)
        ec_copy_xz_point(&S, &P);
        ec_copy_curve(&save_curve, &curve);
        // length = length - 2;

        // Create inputs of form Q,R,R+Q,Q- to test if group morphism
        // where ps[1]=Q, ps[2]=R, ps[3]=Q-R, ps[4]=R+Q
        ec_xz_random_normalized_test(&ps[1], &curve);
        ec_xz_random_normalized_test(&ps[2], &curve);
        ec_xMUL(&ps[1], &Q, &ps[1], TORSION_EVEN, ec_n_bits(TORSION_EVEN, NWORDS_ORDER), &curve);
        ec_xMUL(&ps[2], &Q, &ps[2], TORSION_EVEN, ec_n_bits(TORSION_EVEN, NWORDS_ORDER), &curve);

        ec_xz_to_jac(&psj[1], &ps[1], &curve);
        ec_xz_to_jac(&psj[2], &ps[2], &curve);
        ec_ADD(&psj[4], &psj[1], &psj[2], &curve);
        ec_jac_neg(&Pj, &psj[2]);
        ec_ADD(&psj[3], &psj[1], &Pj, &curve);
        ec_jac_to_xz(&ps[1], &psj[1]);
        ec_jac_to_xz(&ps[2], &psj[2]);
        // switch since xADD seems to prefer doing R,Q,R+Q->R-Q
        ec_jac_to_xz(&ps[3], &psj[4]);
        ec_jac_to_xz(&ps[4], &psj[3]);
        for (int i = 0; i < 5; i++)
            ec_copy_xz_point(&psc[i], &ps[i]);
        ec_xADD(&Q, &ps[1], &ps[2], &ps[3]);
        ec_xADD(&R, &ps[1], &ps[2], &ps[4]);
        assert(ec_xz_is_equal(&Q, &ps[4]));

        // compute chain on that kernel
        // map points though
        ec_copy_curve(&phi.curve, &curve);
        ec_copy_xz_point(&phi.kernel, &P);
        phi.length = length;
        res = res || (int)ec_iso_isogeny_2chain_with_strategy(&imageA, &curve, &P, length, &ps[0], 5, 0);
        ec_curve_init_from_A(&curve, &imageA);

        // check group morphism
        ec_compute_A24(&curve);
        ec_xADD(&Q, &ps[1], &ps[2], &ps[3]);
        ec_xz_normalize(&ps[4]);
        res = res | !ec_x_is_on_curve(&ps[4].x, &curve);
        ec_xz_normalize(&Q);
        res = res | !ec_x_is_on_curve(&Q.x, &curve);
        res = res | !ec_xz_is_equal(&Q, &ps[4]);
        // check kernel point
        ec_xz_normalize(&ps[0]);
        res = res | !ec_x_is_on_curve(&ps[0].x, &curve);
        res = res | !ec_xz_is_zero(&ps[0]);

        ec_copy_curve(&phi.curve, &save_curve);
        ec_copy_xz_point(&phi.kernel, &S);
        phi.length = length;

        // test normalisation

        res = res || (int)ec_iso_isogeny_2chain(&imageA, &phi, 1);
        ec_curve_init_from_A(&curve, &imageA);

        if (res) {
            printf("chain test failed in iteration %d\n", iter);
            break;
        }
        // use curve as next curve
    }
    if (!res) {
        printf("isog chain .......................................................... PASSED\n");
    } else {
        printf("isog chain .......................................................... FAILED\n");
    }
    return (res);
}

int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int iterations = 100;
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
    printf("Testing elliptic curve arithmetic over GF(p):\n\n");
    print_seed(seed);

    // #if defined(TARGET_BIG_ENDIAN)
    // for (int i = 0; i < 12; i++) {
    // seed[i] = BSWAP32(seed[i]);
    //}
    // #endif

    randombytes_init((unsigned char *)seed, NULL, 256);

    res |= ec_test_A24_conversion(iterations);
    res |= ec_test_A24_xDBL(iterations);
    res |= ec_test_generated_basis();
    res |= ec_test_precomputed_basis();
    res |= ec_test_chain(iterations);

    if (res) {
        printf("Tests failed!\n");
    } else {
        printf("All ec isogenies tests passed.\n");
    }

    return res;
}
