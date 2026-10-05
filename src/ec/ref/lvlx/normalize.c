#include <ec.h>
#include <ec_isog.h>
#include <encoded_sizes.h>

// Returns UINT32_MAX if a < b, 0 otherwise
static inline uint32_t
ct_lt_u8(uint8_t a, uint8_t b)
{
    uint32_t borrow_bit = (((uint32_t)a - (uint32_t)b) >> 8) & 1u;
    return (uint32_t)(-(int32_t)borrow_bit);
}

// Returns UINT32_MAX if a == b, 0 otherwise
static inline uint32_t
ct_eq_u8(uint8_t a, uint8_t b)
{
    uint32_t x = (uint32_t)(a ^ b);
    uint32_t nonzero = (x | (uint32_t)(-(int32_t)x)) >> 31;
    return ~((uint32_t)(-(int32_t)nonzero));
}

// Returns UINT32_MAX if x1 < x2, 0 otherwise where x are represented as little endian integers of the form
// encode(x.re) || encode(x.im) where the real and imaginary components themselves are encoded as little endian
// integers.
// In other words, x1 < x2 if im(x1) < im(x2) or im(x1) == im(x2) and re(x1) < re(x2)
int
ec_fp2_less_than(const fp2_t *x1, const fp2_t *x2)
{
    uint8_t buf1[FP2_ENCODED_BYTES];
    uint8_t buf2[FP2_ENCODED_BYTES];

    fp2_encode(buf1, x1);
    fp2_encode(buf2, x2);

    uint32_t result = 0;
    uint32_t all_equal_so_far = -1;

    for (size_t idx = FP2_ENCODED_BYTES; idx-- > 0;) {
        uint32_t less = ct_lt_u8(buf1[idx], buf2[idx]);
        uint32_t equal = ct_eq_u8(buf1[idx], buf2[idx]);
        result |= all_equal_so_far & less;
        all_equal_so_far &= equal;
    }

    return result;
}

void
ec_fp2_mul_by_i(fp2_t *x, const fp2_t *y)
{
    fp_copy(&x->im, &y->re);
    fp_neg(&x->re, &y->im);
}

void
ec_fp2_frob(fp2_t *x, const fp2_t *y)
{
    fp_copy(&x->re, &y->re);
    fp_neg(&x->im, &y->im);
}

void
ec_compute_montgomery_coefficient(ec_xz_point_t *mont, const ec_xz_point_t *th)
{
    fp2_t xx, zz;

    // xx = x^2, zz = z^2
    fp2_sqr(&xx, &th->x);
    fp2_sqr(&zz, &th->z);
    // xx = x^4, zz = y^4
    fp2_sqr(&xx, &xx);
    fp2_sqr(&zz, &zz);

    // A = 2(x^4+z^4)/(x^4-z^4)
    fp2_add(&mont->x, &xx, &zz);
    fp2_sub(&mont->z, &xx, &zz);
    fp2_add(&mont->x, &mont->x, &mont->x);
}

// Given a list of Fp2 element of length len, find max the biggest element and return the index of the maximum in the
// list
void
ec_find_max_coefficient_in_list(fp2_t *max, const fp2_t *list, const uint8_t len)
{
    uint32_t ctl;
    fp2_copy(max, &list[0]);
    for (uint8_t i = 1; i < len; i++) {
        ctl = ec_fp2_less_than(max, &list[i]);
        fp2_select(max, max, &list[i], ctl);
    }
}

bool
ec_theta_to_montgomery(fp2_t *curveA, const ec_xz_point_t *th, const int weil_rest)
{
    ec_xz_point_t mont1, mont2, mont3;
    ec_xz_point_t th2, th3;
    fp2_t ma, mb, ia, ib;
    fp2_t Montgomery_coefs[12];

    fp2_neg(&ma, &th->x);
    fp2_neg(&mb, &th->z);
    ec_fp2_mul_by_i(&ia, &th->x);
    ec_fp2_mul_by_i(&ib, &th->z);

    ec_compute_montgomery_coefficient(&mont1, th);

    // a2,b2 = a+b, a-b
    fp2_add(&th2.x, &th->x, &th->z);
    fp2_sub(&th2.z, &th->x, &th->z);
    ec_compute_montgomery_coefficient(&mont2, &th2);

    // a3,b3 = ia+b, a+ib
    fp2_add(&th3.x, &ia, &th->z);
    fp2_add(&th3.z, &th->x, &ib);
    ec_compute_montgomery_coefficient(&mont3, &th3);

    fp2_copy(&Montgomery_coefs[0], &mont1.z);
    fp2_copy(&Montgomery_coefs[1], &mont2.z);
    fp2_copy(&Montgomery_coefs[2], &mont3.z);
    fp2_batched_inv(Montgomery_coefs, 3);

    uint32_t coef_is_zero = fp2_is_zero(&Montgomery_coefs[0]);
    if (coef_is_zero)
        return 0;

    // Compute all A coefficients on E
    fp2_mul(&Montgomery_coefs[0], &mont1.x, &Montgomery_coefs[0]); // A1
    fp2_mul(&Montgomery_coefs[1], &mont2.x, &Montgomery_coefs[1]); // A2
    fp2_mul(&Montgomery_coefs[2], &mont3.x, &Montgomery_coefs[2]); // A3
    fp2_neg(&Montgomery_coefs[3], &Montgomery_coefs[0]);           // A1'
    fp2_neg(&Montgomery_coefs[4], &Montgomery_coefs[1]);           // A2'
    fp2_neg(&Montgomery_coefs[5], &Montgomery_coefs[2]);           // A3'

    // Compute all A coefficients on E^p
    if (weil_rest == 2) {
        for (uint8_t i = 0; i < 6; i++) {
            ec_fp2_frob(&Montgomery_coefs[6 + i], &Montgomery_coefs[i]); // A_i'' = A_i^p
        }

    }

    ec_find_max_coefficient_in_list(curveA, Montgomery_coefs, 6 * weil_rest);

    return 1;
}

bool
ec_normalize_montgomery(fp2_t *curveA, const ec_xz_point_t *four_torsion, int weil_rest)
{

    if (!ec_theta_to_montgomery(curveA, four_torsion, weil_rest))
        return 0;

    return 1;
}