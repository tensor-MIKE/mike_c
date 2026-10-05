#include <fp.h>
#include <ec.h>
#include <assert.h>
#include <constants.h>

// Adapted from SQIsign 3.0

int
ec_n_bits(const digit_t *k, const int nwords)
{
    int n = 0;
    digit_t bit;
    for (int i = 0; i < nwords; i++) {
        for (int j = 0; j < RADIX; j++) {
            bit = -((k[i] >> j) & 1);
            n = n ^ ((n ^ (i * RADIX + j)) & bit);
        }
    }
    return n;
}

uint32_t
ec_x_is_on_curve(const fp2_t *x, const ec_curve_t *E)
{
    fp2_t t0, t1;
    fp2_sqr(&t0, x);        // t0 = x^2
    fp2_mul(&t1, &E->A, x); // t1 = A*x
    fp2_add(&t0, &t0, &t1); // t0 = x^2+A*x
    fp2_set_one(&t1);       // t1 =1
    fp2_add(&t0, &t0, &t1); // t0 = x^2+A*x+1
    fp2_mul(&t1, &t0, x);   // y2 = x*(x^2+A*x+1)
    return fp2_is_square(&t1);
}

// Helper function which given a point of order k*2^n with n maximal
// and k odd, computes a point of order 2^f
static inline void
ec_clear_cofactor_for_maximal_even_order(ec_xz_point_t *P, const ec_curve_t *curve, const int f)
{
    // clear out the odd cofactor to get a point of order 2^n
    ec_xz_point_t R0;

    // I don't understand why we need the -1, but that is life.
    ec_xMUL(P, &R0, P, p_cofactor_for_2f, P_COFACTOR_FOR_2F_BITLENGTH - 1, curve);

    // clear the power of two to get a point of order 2^f
    for (int i = 0; i < TORSION_EVEN_POWER - f; i++) {
        ec_xDBL(P, P, curve);
    }
}

// Helper function which finds a point x(P) = n * A
static void
ec_find_nA_x_coord(fp2_t *x, const ec_curve_t *curve, const uint8_t start)
{
    assert(!fp2_is_square(&curve->A)); // Only to be called when A is a NQR

    // when A is NQR we allow x(P) to be a multiple n*A of A
    uint8_t n = start;
    if (n == 1) {
        fp2_copy(x, &curve->A);
    } else {
        fp2_mul_small(x, &curve->A, n);
    }

    while (!ec_x_is_on_curve(x, curve)) {
        fp2_add(x, x, &curve->A);
        n++;
    }
}

// The entangled basis generation does not allow A = 0
// so we simply return the one we have already precomputed
static void
ec_basis_E0_2f(ec_xz_point_t *P, ec_xz_point_t *Q, const ec_curve_t *curve, const int f)
{
    assert(fp2_is_zero(&curve->A));

    // Set P, Q to precomputed (X : 1) values
    fp2_copy(&P->x, &BASIS_E0_PX);
    fp2_copy(&Q->x, &BASIS_E0_QX);
    fp2_set_one(&P->z);
    fp2_set_one(&Q->z);

    // clear the power of two to get a point of order 2^f
    for (int i = 0; i < TORSION_EVEN_POWER - f; i++) {
        ec_xDBL(P, P, curve);
        ec_xDBL(Q, Q, curve);
    }
}

// Helper function which finds an NQR -1 / (1 + i*b) for entangled basis generation
static void
ec_find_nqr_factor(fp2_t *x, const ec_curve_t *curve, const uint8_t start)
{
    // factor = -1/(1 + i*b) for b in Fp will be NQR whenever 1 + b^2 is NQR
    // in Fp, so we find one of these and then invert (1 + i*b).
    uint32_t found = 0;
    uint16_t n = start;

    bool qr_b = 1;
    fp_t b, tmp;
    fp2_t z, t0, t1;

    do {
        while (qr_b) {
            // find b with 1 + b^2 a non-quadratic residue
            fp_set_small(&tmp, (uint32_t)n * n + 1);
            qr_b = fp_is_square(&tmp);
            n++; // keeps track of b = n - 1
        }

        // for Px := -A/(1 + i*b) to be on the curve
        // is equivalent to A^2*(z-1) - z^2 NQR for z = 1 + i*b
        // thus prevents unnecessary inversion pre-check

        // t0 = z - 1 = i*b
        // t1 = z = 1 + i*b
        fp_set_small(&b, (uint32_t)n - 1);
        fp2_set_zero(&t0);
        fp2_set_one(&z);
        fp_copy(&z.im, &b);
        fp_copy(&t0.im, &b);

        // A^2*(z-1) - z^2
        fp2_sqr(&t1, &curve->A);
        fp2_mul(&t0, &t0, &t1); // A^2 * (z - 1)
        fp2_sqr(&t1, &z);
        fp2_sub(&t0, &t0, &t1); // A^2 * (z - 1) - z^2
        found = !fp2_is_square(&t0);

        qr_b = 1;
    } while (!found);

    // set Px to -A/(1 + i*b)
    fp2_copy(x, &z);
    fp2_inv(x);
    fp2_mul(x, x, &curve->A);
    fp2_neg(x, x);
}

// Computes a basis E[2^f] = <P, Q> where the point Q is above (0 : 0 : 1)
void
ec_basis_2f(ec_xz_point_t *P, ec_xz_point_t *Q, const ec_curve_t *curve, const int f)
{
    // Should always be normalised
    // ec_compute_A24(curve);

    if (fp2_is_zero(&curve->A)) {
        ec_basis_E0_2f(P, Q, curve, f);
        return;
    }

    bool hint_A = fp2_is_square(&curve->A);

    // Compute the points P, Q
    if (!hint_A) {
        // when A is NQR we allow x(P) to be a multiple n*A of A
        ec_find_nA_x_coord(&P->x, curve, 1);
    } else {
        // when A is QR we instead have to find (1 + b^2) a NQR
        // such that x(P) = -A / (1 + i*b)
        ec_find_nqr_factor(&P->x, curve, 1);
    }

    fp2_set_one(&P->z);
    fp2_add(&Q->x, &curve->A, &P->x);
    fp2_neg(&Q->x, &Q->x);
    fp2_set_one(&Q->z);

    // clear out the odd cofactor to get a point of order 2^f
    ec_clear_cofactor_for_maximal_even_order(P, curve, f);
    ec_clear_cofactor_for_maximal_even_order(Q, curve, f);
}

void
ec_jac_basis_2f(ec_jac_point_t *P, ec_jac_point_t *Q, const ec_curve_t *E, const int f)
{

    ec_xz_point_t Px, Qx;
    ec_jac_point_t PQ, tmp;

    ec_basis_2f(&Px, &Qx, E, f);

    ec_xz_to_jac(P, &Px, E);
    ec_xz_to_jac(&PQ, &Qx, E);

    ec_ADD(&tmp, P, &PQ, E);

    // Ensure Q above (0,0)
    ec_select_jac_point(Q, &tmp, &PQ, fp2_is_zero(&E->A));
}
