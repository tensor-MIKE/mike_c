#ifndef _ISOG_H_
#define _ISOG_H_
#include <mike_namespace.h>
#include <ec.h>

/**
 * @ingroup ec
 * @defgroup iso Isogenies
 * @{
 */

/** @brief data struct consisting of constant required to evaluate points through isogenies of degree 4
 *
 * @typedef ec_kps4_t
 *
 * @struct ec_kps4_t
 */
typedef struct ec_kps4_t
{
    fp2_t K[3];
} ec_kps4_t;

/** @brief An isogeny of degree a power of 2
 *
 * @typedef ec_isog_even_t
 *
 * @struct ec_isog_even_t
 */
typedef struct ec_isog_even_t
{
    ec_curve_t curve;     ///< The domain curve
    ec_xz_point_t kernel; ///< A kernel generator
    unsigned length;      ///< The length as a 2-isogeny walk
} ec_isog_even_t;

/** @brief XZ basis
 *
 * @typedef ec_xz_basis_t
 *
 * @struct  ec_xz_basis_t
 *
 * 3 points P, Q, P-Q
 */
typedef struct ec_xz_basis_t
{
    ec_xz_point_t P;
    ec_xz_point_t Q;
    ec_xz_point_t PmQ;
} ec_xz_basis_t;

/**
 * @brief Four torsion test
 *
 * @param P a point
 * @param E the elliptic curve
 * @return 0xFFFFFFFF if P is 4-torsion but not zero, zero otherwise
 */
uint32_t ec_xz_is_four_torsion(const ec_xz_point_t *P, const ec_curve_t *E);

/**
 * @brief Two torsion test
 *
 * @param P a point
 * @param E the elliptic curve
 * @return 0xFFFFFFFF if P is 2-torsion but not zero, zero otherwise
 */
uint32_t ec_xz_is_two_torsion(const ec_xz_point_t *P, const ec_curve_t *E);

/**
 * @brief Given the point (A+2 : 4C) for a curve, compute the curve coefficients A
 *
 * @param curveA A coefficient of a curve to compute
 * @param A24 the value (A+2 : 4C)
 */
static inline void
ec_A24_to_curve_A(fp2_t *curveA, const ec_xz_point_t *A24)
{
    // (A:C) = ((A+2C)*2-4C : 4C)
    fp2_add(curveA, &A24->x, &A24->x);
    fp2_sub(curveA, curveA, &A24->z);
    fp2_add(curveA, curveA, curveA);
    // curve normalization
    fp2_t inv;
    fp2_copy(&inv, &A24->z);
    fp2_inv(&inv);
    fp2_mul(curveA, curveA, &inv);
}

/**
 * @brief Given a curve the point (A+2 : 4C) compute the curve coefficients (A : C)
 *
 * @param A24 Output: A point to be set to the value (A+2 : 4)
 * @param curve a curve
 */
static inline void
ec_curve_to_A24(ec_xz_point_t *A24, const ec_curve_t *curve)
{
    fp2_t cst;
    fp2_set_small(&cst, 2);
    fp2_add(&A24->x, &curve->A, &cst);
    fp2_set_small(&A24->z, 4);
    fp2_inv(&A24->z);
    fp2_mul(&A24->x, &A24->z, &A24->x);
    fp2_set_one(&A24->z);
}

/**
 * @brief Doubling on curves represented by as point
 *
 * @param Q Output: double of P
 * @param A24 Output: Point representing a curve as in ec_curve_to_A24
 * @param P input point
 */
void ec_xDBL_A24(ec_xz_point_t *Q, const ec_xz_point_t *P, const ec_xz_point_t *A24);

// Normalize the Montgomery coefficient of a curve

/**
 * @brief Transform a curves A coefficient and a basis B into the curve's normalised A: the biggest A coefficient
 *
 * @param curveA: (output) maximal A coefficent
 * @param four_torsion: (input) a 4 torsion point of E, not above (0,0).
 * @param normalize: (input) 1 = normalize on E, 2 = normalize on E and E^p
 *
 * @return 0 if not valid.
 */
bool ec_normalize_montgomery(fp2_t *curveA, const ec_xz_point_t *four_torsion, int normalize);

/**
 * @brief Evaluate isogeny of even degree and return normalised codomain
 * Returns 0 if successful and -1 if kernel has the wrong order or includes (0:1).
 *
 * @param imageA Output: A coefficient of the codomain
 * @param phi isogeny, conissting in points
 * @param normalize 0, no normalisation, 1 normalised, 2 normalisation on E and its codomain E^p
 *
 * @return 0 if there was no error, 0xFFFFFFFF otherwise. Some input errors are caught only by asserts
 */
uint32_t ec_iso_isogeny_2chain(fp2_t *imageA, ec_isog_even_t *phi, int normalize);

/**
 * @brief Evaluate isogeny of even degree on list of points.
 * Returns 0 if successful and -1 if kernel has the wrong order or includes (0:1).
 *
 * @param imageA Output: A coefficients of the output
 * @param curve Domain curve
 * Compute_A24 must be called before any computation is on the output curve
 * @param kernel Kernel of the isogeny to be computed
 * @param isog_len Length of the isogeny  (2^isog_len must divide order of kernel)
 * @param to_eval Input/Output Points to be mapped through the isogeny
 * @param to_eval_len Number of points to be mapped through
 * @param normalize If 1 normalize the output A coefficient, else not. Must be 0 if to_eval_len>0
 *
 * @return 0 if there was no error, 0xFFFFFFFF otherwise. Some input errors are caught only by asserts
 */
uint32_t ec_iso_isogeny_2chain_with_strategy(fp2_t *imageA,
                                             const ec_curve_t *curve,
                                             const ec_xz_point_t *kernel,
                                             const int isog_len,
                                             ec_xz_point_t *to_eval,
                                             int to_eval_len,
                                             int normalize);

/**
 * @brief Evaluate a point through a 4 isogeny
 *
 * @param R image point
 * @param Q input point
 * @param lenQ number of evalaution
 * @param kps precomputed data to evaluate the isogeny
 *
 */
void ec_iso_xeval_4(ec_xz_point_t *R, const ec_xz_point_t *Q, const int lenQ, const ec_kps4_t *kps);

/**
 * @brief compute a 4 isogeny chain of kernel P
 *
 * @param kps   precomputed data
 * @param B     codomain (A:C) coefficient
 * @param P     kernel of isogeny
 *
 */
void ec_iso_xisog_4(ec_kps4_t *kps, ec_xz_point_t *B, const ec_xz_point_t P);

/**
 * @}
 */

#endif