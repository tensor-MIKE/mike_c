#include "test_extras.h"
#include "rng.h"
#include <ec.h>
#include <ec_params.h>
#include <encoded_sizes.h>
#include <constants.h>

// Make n random-ish field elements (for tests only!).
void
ec_fp_random_test(fp_t *a)
{
    uint8_t tmp[FP_ENCODED_BYTES];

    randombytes(tmp, sizeof(tmp));

    fp_decode_reduce(a, tmp, sizeof(tmp));
}

void
ec_fp2_random_test(fp2_t *a)
{
    ec_fp_random_test(&a->im);
    ec_fp_random_test(&a->re);
}

uint32_t
ec_ec_point_is_on_curve(fp2_t *y2, const fp2_t *x, const ec_curve_t *E)
{
    fp2_t t0, t1;
    fp2_sqr(&t0, x);        // t0 = x^2
    fp2_mul(&t1, &E->A, x); // t1 = A*x
    fp2_add(&t0, &t0, &t1); // t0 = x^2+A*x
    fp2_set_one(&t1);       // t1 =1
    fp2_add(&t0, &t0, &t1); // t0 = x^2+A*x+1
    fp2_mul(y2, &t0, x);    // y2 = x*(x^2+A*x+1)
    return fp2_is_square(y2);
}

void
ec_xz_random_normalized_test(ec_xz_point_t *P, const ec_curve_t *curve)
{
    fp2_set_one(&P->z);
    fp2_t y2;
    while (1) {
        ec_fp2_random_test(&P->x);
        if (ec_ec_point_is_on_curve(&y2, &P->x, curve)) {
            break;
        }
    }
}

void
ec_xz_random_test(ec_xz_point_t *P, const ec_curve_t *curve)
{
    ec_xz_random_normalized_test(P, curve);
    ec_fp2_random_test(&P->z);
    fp2_mul(&P->x, &P->x, &P->z);
}