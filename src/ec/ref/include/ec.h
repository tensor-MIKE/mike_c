/** @file
 *
 * @authors Pierrick Dartois, based on
 * work by Luca De Feo and Francisco RH
 * for the SQIsign v2.0 NIST submission
 *
 * @brief Elliptic curve stuff
 */

#ifndef EC_H
#define EC_H

#include <fp2.h>
#include <stdio.h>
#include <fp_constants.h>
#include <ec_params.h>

#include "e0_basis.h"

/** @defgroup ec Elliptic curves
 * @{
 */

/** @defgroup ec_t Data structures
 * @{
 */

/** @brief Projective point in (X:Z) Montgomery
 * coordinates.
 *
 * @typedef ec_xz_point_t
 *
 * @struct ec_xz_point_t
 *
 * A projective point (X:Z) representing a point
 * on the Kummer line.
 */
typedef struct ec_xz_point_t
{
    fp2_t x;
    fp2_t z;
} ec_xz_point_t;

/** @brief Projective point in Jacobian coordinates
 * in the Montgomery model y^2 = x^3 + Ax^2 + x
 *
 * @typedef ec_jac_point_t
 *
 * @struct ec_jac_point_t
 *
 * A projective point in (X:Y:Z) coordinates such that
 * x = X/Z^2 and y = X/Z^3.
 */
typedef struct ec_jac_point_t
{
    fp2_t x;
    fp2_t y;
    fp2_t z;
} ec_jac_point_t;

/** @brief Projective point in Jacobian coordinates
 * in the Weierstrass model y^2 = x^3 + ax + b
 * with a = 1-A^2/3 and b = A/3(2A^2/9-1)
 *
 * @typedef ec_jac_ws_point_t
 *
 * @struct ec_jac_ws_point_t
 *
 * A projective point in (X:Y:Z:T) coordinates such that
 * x = X/Z^2, y = X/Z^2 and T=a*Z^4
 */
typedef struct ec_jac_ws_point_t
{
    fp2_t x;
    fp2_t y;
    fp2_t z;
    fp2_t t;
} ec_jac_ws_point_t;

/** @brief Addition components
 *
 * @typedef ec_bary_coords_t
 *
 * @struct ec_bary_coords_t
 *
 * 3 components u,v,w that define the (X:Z) coordinates of both
 * addition and substraction of two distinct points with
 * P+Q =(u-v:w) and P-Q = (u+v:w)
 */
typedef struct ec_bary_coords_t
{
    fp2_t u;
    fp2_t v;
    fp2_t w;
} ec_bary_coords_t;

/** @brief An elliptic curve
 *
 * @typedef ec_curve_t
 *
 * @struct ec_curve_t
 *
 * An elliptic curve in projective Montgomery form
 */
typedef struct ec_curve_t
{
    fp2_t A;
    fp2_t A24; // A24 = (A+2)/4
} ec_curve_t;

typedef struct ec_ws_curve_t
{
    // Coefficients for the Weierstrass Jacobian model
    fp2_t ao3; // ao3 = A/3
    fp2_t a;   // a = 1-A^2/3
} ec_ws_curve_t;

// end ec_t
/** @}
 */

/** @defgroup ec_curve_t Curves
 * @{
 */

// Initalisation for curves
void ec_curve_init_precomputed(ec_curve_t *E);

/**
 * @brief Initialize an elliptic curve from a coefficient
 *
 * @param A an fp2_t
 * @param E the elliptic curve to initialize
 */
void ec_curve_init_from_A(ec_curve_t *E, const fp2_t *A);

static inline void
ec_copy_xz_point(ec_xz_point_t *P, const ec_xz_point_t *Q)
{
    fp2_copy(&P->x, &Q->x);
    fp2_copy(&P->z, &Q->z);
}

static inline void
ec_copy_jac_point(ec_jac_point_t *P, const ec_jac_point_t *Q)
{
    fp2_copy(&P->x, &Q->x);
    fp2_copy(&P->y, &Q->y);
    fp2_copy(&P->z, &Q->z);
}

static inline void
ec_copy_curve(ec_curve_t *E1, const ec_curve_t *E2)
{
    fp2_copy(&(E1->A), &(E2->A));
    fp2_copy(&(E1->A24), &(E2->A24));
}

// Selecting points and curves
static inline void
ec_select_jac_point(ec_jac_point_t *Q, const ec_jac_point_t *P1, const ec_jac_point_t *P2, const uint32_t option)
{ // Select points in constant time
  // If option = 0 then Q <- P1, else if option = 0xFF...FF then Q <- P2
    fp2_select(&(Q->x), &(P1->x), &(P2->x), option);
    fp2_select(&(Q->y), &(P1->y), &(P2->y), option);
    fp2_select(&(Q->z), &(P1->z), &(P2->z), option);
}

static inline void
ec_cswap_jac_points(ec_jac_point_t *P, ec_jac_point_t *Q, const uint32_t option)
{ // Swap points in constant time
  // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then P <- Q and Q <- P
    fp2_cswap(&(P->x), &(Q->x), option);
    fp2_cswap(&(P->y), &(Q->y), option);
    fp2_cswap(&(P->z), &(Q->z), option);
}

static inline void
ec_cswap_xz_points(ec_xz_point_t *P, ec_xz_point_t *Q, const uint32_t option)
{ // Swap points in constant time
  // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then P <- Q and Q <- P
    fp2_cswap(&(P->x), &(Q->x), option);
    fp2_cswap(&(P->z), &(Q->z), option);
}

/**
 * @brief Compute Jacobian Weierstrass coefficients ao3 and a
 *
 * @param Ews a curve with Jacobian Weierstrass coefficients ao3 and a
 * @param E a curve
 */
void ec_compute_ws(ec_ws_curve_t *Ews, const ec_curve_t *E);

/**
 * @brief Compute Montgomery coefficients A24 = (A+2)/4
 *
 * @param E a curve
 */
void ec_compute_A24(ec_curve_t *E);

// Conversion functions

/**
 * @brief Converts a point from xz-only projective coordinates
 * to Jacobian coordinates in the Montgomery model, choosing sign arbitrarily
 *
 * @param Q output converted point
 * @param P input point to convert
 * @param E curve
 */
void ec_xz_to_jac(ec_jac_point_t *Q, const ec_xz_point_t *P, const ec_curve_t *E);

/**
 * @brief Converts points from Montgomery to Weierstrass Jacobian
 * coordinates
 *
 * @param Q output converted point
 * @param P input point to convert
 * @param E parent curve
 *
 * Computes the additional coordinate T = a*Z^4
 */
void ec_jac_to_ws(ec_jac_ws_point_t *Q, const ec_jac_point_t *P, const ec_ws_curve_t *E);

/**
 * @brief Converts points from Weierstrass to Montgomery Jacobian
 * coordinates
 *
 * @param Q output converted point
 * @param P input point to convert
 * @param E parent curve
 *
 * Forgets the T coordinate.
 */
void ec_ws_to_jac(ec_jac_point_t *Q, const ec_jac_ws_point_t *P, const ec_ws_curve_t *E);

/**
 * @brief Converts points from  Jacobian (X:Y:Z) to (X:Z) Montgomery
 * coordinates
 *
 * @param Q output converted point
 * @param P input point to convert
 *
 * Forgets the Y coordinate.
 */
void ec_jac_to_xz(ec_xz_point_t *Q, const ec_jac_point_t *P);

/**
 * @brief Normalize xz to x/z,1 or 1,0
 *
 * @param P: (input/output) xz point
 *
 */
void ec_xz_normalize(ec_xz_point_t *P);

/** @}
 */

/** @defgroup ec_point_t Point operations
 * @{
 */

/**
 * @brief Point equality
 *
 * @param P a point in Jacobian coordinates
 * @param Q a point in Jacobian coordinates
 * @return 0xFFFFFFFF if equal, zero otherwise
 */
uint32_t ec_jac_is_equal(const ec_jac_point_t *P, const ec_jac_point_t *Q);

/**
 * @brief Point at infinity equality
 *
 * @param P a point in Jacobian coordinates
 * @return 0xFFFFFFFF if point at infinity, zero otherwise
 */
uint32_t ec_jac_is_zero(const ec_jac_point_t *P);

/**
 * @brief Point equality
 *
 * @param P a point in (X:Z) Montgomery coordinates
 * @param Q a point in (X:Z) Montgomery coordinates
 * @return 0xFFFFFFFF if equal, zero otherwise
 */
uint32_t ec_xz_is_equal(const ec_xz_point_t *P, const ec_xz_point_t *Q);

/**
 * @brief Point at infinity equality
 *
 * @param P a point in (X:Z) Montgomery coordinates
 * @return 0xFFFFFFFF if point at infinity, zero otherwise
 */
uint32_t ec_xz_is_zero(const ec_xz_point_t *P);

// Jacobian arithmetic

/**
 * @brief Negation of a point in Jacobian coordinates in the
 * Montgomery model.
 *
 * @param Q output negated point
 * @param P input point to negate
 */
void ec_jac_neg(ec_jac_point_t *Q, const ec_jac_point_t *P);

/**
 * @brief Sums two points in Jacobian coordinates in the
 * Montgomery model.
 *
 * @param R output sum of points P+Q
 * @param P input point to sum
 * @param Q input point to sum
 * @param E parent elliptic curve
 *
 * Complete algorithm that can handle all edge cases (P==0, Q==0, P=+/- Q...)
 */
void ec_ADD(ec_jac_point_t *R, const ec_jac_point_t *P, const ec_jac_point_t *Q, const ec_curve_t *E);

/**
 * @brief Doubles a point in Jacobian coordinates in the Montgomery model.
 *
 * @param Q output point duplication 2P
 * @param P input point to double
 * @param E parent elliptic curve
 *
 * Complete algorithm that can handle the case P == 0.
 */
void ec_DBL(ec_jac_point_t *Q, const ec_jac_point_t *P, const ec_curve_t *E);

/**
 * @brief Doubles a point in Jacobian coordinates in the Weierestrass model.
 *
 * @param Q output point duplication 2P
 * @param P input point to double
 *
 * Complete algorithm that can handle the case P == 0.
 */
void ec_DBL_ws(ec_jac_ws_point_t *Q, const ec_jac_ws_point_t *P);

/**
 * @brief Doubles a point multiple time in Jacobian coordinates in
 * the Jacobian model.
 *
 * @param Q output point(s) duplication [2^n]P
 * @param P input point(s) to double n times
 * @param n number of times to double
 * @param E parent elliptic curve
 * @param number_of_points length of the point arrays P and Q (must be equal)
 *
 * Complete algorithm that can handle the case P == 0. Optimised
 * using the Weierstrass model when multiple iterations are needed.
 */
void ec_jac_dbl_iter(ec_jac_point_t *Q,
                     const ec_jac_point_t *P,
                     const unsigned int n,
                     const ec_curve_t *E,
                     const int number_of_points);

/**
 * @brief The Montgomery ladder with Jacobian projective
 * coordinates in the Montgomery model
 *
 * @param Q output point [k]P
 * @param P input point
 * @param k scalar multiple
 * @param kbits number of bits of k
 * @param curve parent curve
 *
 * Uses xMUL and y reconstruction (revover_y) and converts back
 * to Jacobian coordinates.
 */
void ec_jac_MUL(ec_jac_point_t *Q, const ec_jac_point_t *P, const digit_t *k, const int kbits, const ec_curve_t *curve);

/**
 * @brief The Montgomery ladder with Jacobian projective coordinates in the Montgomery model
 * for scalar k>>k_right_shift
 *
 * @param Q output point [k]P
 * @param P input point
 * @param k scalar multiple
 * @param kbits number of bits of k
 * @param curve parent curve
 * @param k_right_shift Compute ladder with scalar k>>k_right_shift
 *
 * Uses xMUL and y reconstruction (revover_y) and converts back
 * to Jacobian coordinates.
 */
void ec_jac_MUL_shifted(ec_jac_point_t *Q,
                        const ec_jac_point_t *P,
                        const digit_t *k,
                        const int kbits,
                        const ec_curve_t *curve,
                        const int k_right_shift);

/**
 * @brief Take P and Q in E distinct, two jac_point_t, and returns
 * the barycentric coordinates i.e. three components u,v and w in
 * the base fp such that the (X:Z) coordinates of P+Q are (u-d_P*d_Q*v:w)
 * and of P-Q are (u+d_P*d_Q*v:w) on E, where d_P = i if P is on the twist of E
 * and 1 otherwise and d_Q = i if Q is on the twist of E and 1 otherwise
 * (see Lemma 1 and Proposition 2 of https://eprint.iacr.org/2025/736.pdf).
 *
 * @param uvw output barycentric coordinates
 * @param P input point
 * @param Q input point
 * @param E parent elliptic curve (or its twist)
 *
 * Now also works when P == 0 or Q == 0 (but not both).
 */
void ec_jac_to_bary_coords(ec_bary_coords_t *uvw,
                           const ec_jac_point_t *P,
                           const ec_jac_point_t *Q,
                           const ec_curve_t *E);

// xz-aritmetic

/**
 * @brief Doubles a point in (X:Z) coordinates in the Montgomery model.
 *
 * @param Q output point duplication 2P
 * @param P input point to double
 * @param E parent elliptic curve
 *
 * Complete algorithm that can handle the case P == 0.
 */
void ec_xDBL(ec_xz_point_t *Q, const ec_xz_point_t *P, const ec_curve_t *E);

/**
 * @brief Differential addition of points P, Q, P-Q in (X:Z) coordinates
 * in the Montgomery model.
 *
 * @param R output sum of points P+Q
 * @param P input point to sum
 * @param Q input point to sum
 * @param PQ difference P-Q
 *
 * Correct in the case P == 0 or Q == 0 but not P == Q (PQ == 0).
 */
void ec_xADD(ec_xz_point_t *R, const ec_xz_point_t *P, const ec_xz_point_t *Q, const ec_xz_point_t *PQ);

/**
 * @brief Simultaneous doubling and differential addition of
 * points in (X:Z) coordinates in the Montgomery model.
 *
 * @param R output point duplication 2P
 * @param S sum of points P+Q
 * @param P input point to sum
 * @param Q input point to sum
 * @param PQ difference P-Q
 * @param E parent curve
 *
 * Correct in the case P == 0 or Q == 0 but not P == Q (PQ == 0).
 */
void ec_xDBLADD(ec_xz_point_t *R,
                ec_xz_point_t *S,
                const ec_xz_point_t *P,
                const ec_xz_point_t *Q,
                const ec_xz_point_t *PQ,
                const ec_curve_t *E);

/**
 * @brief The Montgomery ladder with (X:Z) Montgomery coordinates
 *
 * @param Q output point [k]P
 * @param R output point [k+1]P
 * @param P input point
 * @param k scalar multiple
 * @param kbits number of bits of k
 * @param curve parent curve
 *
 * Uses https://eprint.iacr.org/2017/212, Algorithm 6.
 */
void ec_xMUL(ec_xz_point_t *Q,
             ec_xz_point_t *R,
             const ec_xz_point_t *P,
             const digit_t *k,
             const int kbits,
             const ec_curve_t *curve);

/**
 * @brief The Montgomery ladder with (X:Z) Montgomery coordinates
 *
 * @param Q output point [k]P
 * @param R output point [k+1]P
 * @param P input point
 * @param k scalar multiple is k>>k_right_shift
 * @param kbits number of bits of k
 * @param curve parent curve
 * @param k_right_shift use ladder with k>>k_right_shift as scalar
 *
 * Uses https://eprint.iacr.org/2017/212, Algorithm 6.
 */
void ec_xMUL_shifted(ec_xz_point_t *Q,
                     ec_xz_point_t *R,
                     const ec_xz_point_t *P,
                     const digit_t *k,
                     const int kbits,
                     const ec_curve_t *curve,
                     const int k_right_shift);

/**
 * @brief Computation of P + [m]*Q
 *
 * @param R computed P + m * Q
 * @param k an unsigned multi-precision integer
 * @param kbits bit length of k
 * @param P a point
 * @param Q a point
 * @param PQ the difference P-Q
 * @param E the curve
 * @return 0 if there was an error, 1 otherwise
 */
int ec_ladder3pt(ec_xz_point_t *R,
                 const ec_xz_point_t *P,
                 const ec_xz_point_t *Q,
                 const ec_xz_point_t *PQ,
                 const digit_t *k,
                 const int kbits,
                 const ec_curve_t *E);
/**
 * @brief Computes a full projective jacobian point Q=(X_Q:Y_Q:Z_Q)
 * from the Montgomery (X:Z) coordinates of Q and P+Q and the full jacobian
 * coordinates (X_P:Y_P:Z_P) of P.
 *
 * @param R output point to recover
 * @param P input point in Montgomery (X:Y:Z) coordinates
 * @param Q input point in Montgomery (X:Z) coordinates
 * @param PQ sum P+Q in Montgomery (X:Z) coorfinates
 * @param E parent curve
 *
 * Uses Okeya–Sakurai's algorithm (DOI 10.1007/3-540-44709-1_12,
 * see also https://eprint.iacr.org/2017/212, Algorithm 5). Complete
 * algorithm that handles edge cases (P, Q or PQ ==0).
 */
void ec_recover_y(ec_jac_point_t *R,
                  const ec_jac_point_t *P,
                  const ec_xz_point_t *Q,
                  const ec_xz_point_t *PQ,
                  const ec_curve_t *E);

/** @}
 */

/** @defgroup Basis Basis generation
 * @{
 */

/**
 * @brief Determines if the x-coordinate defines an Fp-rational point.
 *
 * @param x input affine x-coordinate
 * @param E Montgomery curve y^2 = x^3 + A*x^2 + x
 *
 * @return 0xF...F if x defines an Fp-rational point and 0 otherwise.
 *
 * Also computes y2 to compute square roots afterwards.
 */
uint32_t ec_x_is_on_curve(const fp2_t *x, const ec_curve_t *E);

/**
 * @brief Computes basis <P,Q> of E[2^f]
 *
 * @param P output point
 * @param Q output point
 * @param E parent curve
 * @param f 2-valaution
 */
void ec_jac_basis_2f(ec_jac_point_t *P, ec_jac_point_t *Q, const ec_curve_t *E, const int f);

/**
 * @brief Computes basis <P,Q> of E[2^f]
 *
 * @param P output point
 * @param Q output point
 * @param curve parent curve
 * @param f 2-valaution
 */
void ec_basis_2f(ec_xz_point_t *P, ec_xz_point_t *Q, const ec_curve_t *curve, const int f);

// internal function, but not only used in tests
int ec_n_bits(const digit_t *k, const int nwords);
/** @}
 */

/** @}
 */
#endif
