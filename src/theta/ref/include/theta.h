#ifndef THETA_H
#define THETA_H

#include <fp.h>
#include <stdint.h>
#include <stdio.h>

extern const uint8_t theta_to_weil[16]; // = {0, 1, 2, 3, 2, 4, 5, 6, 1, 7, 4, 8, 3, 8, 6, 9};
extern const uint8_t weil_to_theta[10]; // = {0, 1, 2, 3, 5, 6, 7, 9, 11 ,15};

/** @defgroup Theta Theta structures
 * @{
 */

/** @defgroup Struct Structures
 * @{
 */

/** @brief Projective theta point in dimension 4
 *
 * @typedef theta_point_dim4_t
 *
 * @struct theta_point_dim4_t
 *
 * Given by 16 coordinates in level 2
 */
typedef fp_t theta_point_dim4_t[16];

/** @brief theta point in dimension 4 of Weil type
 *
 * @typedef theta_point_weil_t
 *
 * @struct theta_point_weil_t
 *
 * Given by 10 distinct coordinates in level 2
 */
typedef fp_t theta_point_weil_t[10];

/** @brief Projective theta point in dimension 2
 *
 * @typedef theta_point_dim2_t
 *
 * @struct theta_point_dim2_t
 *
 * Given by 4 coordinates in level 2
 */
typedef fp_t theta_point_dim2_t[4];

/** @brief Theta structure in dimension 4 for Weil type
 *
 * @typedef theta_struct_weil_t
 *
 * @struct theta_struct_weil_t
 */
typedef struct theta_struct_weil_t
{
    theta_point_weil_t inv_dual_null_point; // Stores the inverse dual theta null point
                                            // by default (for isogeny evaluations).
} theta_struct_weil_t;

/** @brief Theta structure in dimension 4 for Weil type
 *
 * @typedef theta_struct_weil_precomp_t
 *
 * @struct theta_struct_weil_precomp_t
 */
typedef struct theta_struct_weil_precomp_t
{
    theta_point_weil_t inv_null_point_DBL;
    theta_point_weil_t inv_HS_null_point_DBL;

} theta_struct_weil_precomp_t;

/** @brief Theta structure in dimension 2
 *
 * @typedef theta_struct_dim2_t
 *
 * @struct theta_struct_dim2_t
 */
typedef struct theta_struct_dim2_t
{
    theta_point_dim2_t inv_dual_null_point; // Stores the inverse dual theta null point
                                            // by default (for isogeny evaluations).
    theta_point_dim2_t dual_null_point;
} theta_struct_dim2_t;

/** @brief precomputation used for Theta structure in dimension 4
 *
 * @typedef theta_struct_dim2_precomp_t
 *
 * @struct theta_struct_dim2_precomp_t
 */
typedef struct theta_struct_dim2_precomp_t
{
    theta_point_dim2_t inv_null_point_DBL;
    theta_point_dim2_t inv_HS_null_point_DBL;
} theta_struct_dim2_precomp_t;

/** // end Struct
 * @}
 */

/** @defgroup  dim4  Theta points dim 4
 * @{
 */

/** @defgroup  Theta-points 
 * @{
 */
/**
 * @brief transform a Weil-point into a 16-element theta point.
 *
 * @param out (output) the theta_point_dim4_t point.
 * @param in (input) the theta_point_weil_t point.
 *
 */
static inline void
theta_weil_to_dim4_point(theta_point_dim4_t *out, const theta_point_weil_t *in)
{
    for (uint8_t i = 0; i < 10; i++) {
        fp_copy(&(*out)[weil_to_theta[i]], &(*in)[i]);
    }

    // Complete missing coordinates
    fp_copy(&(*out)[4], &(*out)[2]);
    fp_copy(&(*out)[8], &(*out)[1]);
    fp_copy(&(*out)[10], &(*out)[5]);
    fp_copy(&(*out)[12], &(*out)[3]);
    fp_copy(&(*out)[13], &(*out)[11]);
    fp_copy(&(*out)[14], &(*out)[7]);
}

/**
 * @brief transform a dim 4 theta-point into a Weil point
 *
 * @param out (output) the theta_point_weil_t point.
 * @param in (input) the theta_point_dim4_t point.
 *
 */
static inline void
theta_dim4_to_weil(theta_point_weil_t *out, const theta_point_dim4_t *in)
{
    for (uint8_t i = 0; i < 10; i++) {
        fp_copy(&(*out)[i], &(*in)[weil_to_theta[i]]);
    }
}

/**
 * @brief Hadamard transform on thetapoint of dim 4.
 *
 * @param out (output) The hadamard transformed thetapoint.
 * @param in (input) The input thetapoint of dim 4.
 *
 */
static inline void
theta_hadamard_dim4(theta_point_dim4_t *out, const theta_point_dim4_t *in)
{
    fp_t arr[16], arr2[16];

    // len = 1
    fp_add(&arr[0], &(*in)[0], &(*in)[1]);
    fp_sub(&arr[1], &(*in)[0], &(*in)[1]);
    fp_add(&arr[2], &(*in)[2], &(*in)[3]);
    fp_sub(&arr[3], &(*in)[2], &(*in)[3]);
    fp_add(&arr[4], &(*in)[4], &(*in)[5]);
    fp_sub(&arr[5], &(*in)[4], &(*in)[5]);
    fp_add(&arr[6], &(*in)[6], &(*in)[7]);
    fp_sub(&arr[7], &(*in)[6], &(*in)[7]);
    fp_add(&arr[8], &(*in)[8], &(*in)[9]);
    fp_sub(&arr[9], &(*in)[8], &(*in)[9]);
    fp_add(&arr[10], &(*in)[10], &(*in)[11]);
    fp_sub(&arr[11], &(*in)[10], &(*in)[11]);
    fp_add(&arr[12], &(*in)[12], &(*in)[13]);
    fp_sub(&arr[13], &(*in)[12], &(*in)[13]);
    fp_add(&arr[14], &(*in)[14], &(*in)[15]);
    fp_sub(&arr[15], &(*in)[14], &(*in)[15]);

    // len = 2
    fp_add(&arr2[0], &arr[0], &arr[2]);
    fp_sub(&arr2[2], &arr[0], &arr[2]);
    fp_add(&arr2[1], &arr[1], &arr[3]);
    fp_sub(&arr2[3], &arr[1], &arr[3]);
    fp_add(&arr2[4], &arr[4], &arr[6]);
    fp_sub(&arr2[6], &arr[4], &arr[6]);
    fp_add(&arr2[5], &arr[5], &arr[7]);
    fp_sub(&arr2[7], &arr[5], &arr[7]);
    fp_add(&arr2[8], &arr[8], &arr[10]);
    fp_sub(&arr2[10], &arr[8], &arr[10]);
    fp_add(&arr2[9], &arr[9], &arr[11]);
    fp_sub(&arr2[11], &arr[9], &arr[11]);
    fp_add(&arr2[12], &arr[12], &arr[14]);
    fp_sub(&arr2[14], &arr[12], &arr[14]);
    fp_add(&arr2[13], &arr[13], &arr[15]);
    fp_sub(&arr2[15], &arr[13], &arr[15]);

    // len = 4
    fp_add(&arr[0], &arr2[0], &arr2[4]);
    fp_sub(&arr[4], &arr2[0], &arr2[4]);
    fp_add(&arr[1], &arr2[1], &arr2[5]);
    fp_sub(&arr[5], &arr2[1], &arr2[5]);
    fp_add(&arr[2], &arr2[2], &arr2[6]);
    fp_sub(&arr[6], &arr2[2], &arr2[6]);
    fp_add(&arr[3], &arr2[3], &arr2[7]);
    fp_sub(&arr[7], &arr2[3], &arr2[7]);
    fp_add(&arr[8], &arr2[8], &arr2[12]);
    fp_sub(&arr[12], &arr2[8], &arr2[12]);
    fp_add(&arr[9], &arr2[9], &arr2[13]);
    fp_sub(&arr[13], &arr2[9], &arr2[13]);
    fp_add(&arr[10], &arr2[10], &arr2[14]);
    fp_sub(&arr[14], &arr2[10], &arr2[14]);
    fp_add(&arr[11], &arr2[11], &arr2[15]);
    fp_sub(&arr[15], &arr2[11], &arr2[15]);

    // len = 8
    fp_add(&(*out)[0], &arr[0], &arr[8]);
    fp_sub(&(*out)[8], &arr[0], &arr[8]);
    fp_add(&(*out)[1], &arr[1], &arr[9]);
    fp_sub(&(*out)[9], &arr[1], &arr[9]);
    fp_add(&(*out)[2], &arr[2], &arr[10]);
    fp_sub(&(*out)[10], &arr[2], &arr[10]);
    fp_add(&(*out)[3], &arr[3], &arr[11]);
    fp_sub(&(*out)[11], &arr[3], &arr[11]);
    fp_add(&(*out)[4], &arr[4], &arr[12]);
    fp_sub(&(*out)[12], &arr[4], &arr[12]);
    fp_add(&(*out)[5], &arr[5], &arr[13]);
    fp_sub(&(*out)[13], &arr[5], &arr[13]);
    fp_add(&(*out)[6], &arr[6], &arr[14]);
    fp_sub(&(*out)[14], &arr[6], &arr[14]);
    fp_add(&(*out)[7], &arr[7], &arr[15]);
    fp_sub(&(*out)[15], &arr[7], &arr[15]);
}

/**
 * @brief Hadamard transform on compress weil points.
 *
 * @param[out] out The hadamard transformed vector.
 * @param[in] in The vector to hadamard transform.
 */
static inline void
theta_hadamard_weil(theta_point_weil_t *out, const theta_point_weil_t *in)
{

    fp_t t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
    fp_t s0, s1, s2, s3, s4, s5;

    // 10a
    fp_add(&t0, &(*in)[0], &(*in)[9]);
    fp_sub(&t1, &(*in)[0], &(*in)[9]);
    fp_add(&t2, &(*in)[5], &(*in)[7]);
    fp_sub(&t3, &(*in)[5], &(*in)[7]);
    fp_add(&t4, &(*in)[2], &(*in)[8]);
    fp_sub(&t5, &(*in)[2], &(*in)[8]);
    fp_add(&t6, &(*in)[1], &(*in)[6]);
    fp_sub(&t7, &(*in)[1], &(*in)[6]);
    fp_add(&t8, &(*in)[3], &(*in)[4]);
    fp_sub(&t9, &(*in)[3], &(*in)[4]);

    // 6a
    fp_add(&s0, &t0, &t2);
    fp_sub(&s1, &t0, &t2);
    fp_add(&s2, &t1, &t3);
    fp_sub(&s3, &t1, &t3);
    fp_add(&s4, &t6, &t4);
    fp_sub(&s5, &t6, &t4);

    // 6a DBL
    fp_add(&t5, &t5, &t5);
    fp_add(&t7, &t7, &t7);
    fp_add(&t8, &t8, &t8);
    fp_add(&t9, &t9, &t9);
    fp_add(&s4, &s4, &s4);
    fp_add(&s5, &s5, &s5);

    fp_add(&t0, &s0, &t8);
    fp_sub(&t1, &s0, &t8);

    // write in output
    fp_add(&(*out)[0], &t0, &s4);
    fp_sub(&(*out)[9], &t0, &s4);
    fp_add(&(*out)[1], &s2, &t5);
    fp_sub(&(*out)[6], &s2, &t5);
    fp_add(&(*out)[2], &s3, &t7);
    fp_sub(&(*out)[8], &s3, &t7);
    fp_add(&(*out)[3], &s1, &t9);
    fp_sub(&(*out)[4], &s1, &t9);
    fp_add(&(*out)[5], &t1, &s5);
    fp_sub(&(*out)[7], &t1, &s5);
}

/**
 * @brief square all coefficient of thetapoint of dim 4.
 *
 * @param out (output) The squared thetapoint
 * @param in (input) The input thetapoint of dim 4.
 *
 */
static inline void
theta_squared_dim4(theta_point_dim4_t *out, const theta_point_dim4_t *in)
{
    for (uint8_t i = 0; i < 16; i++) {
        fp_sqr(&(*out)[i], &(*in)[i]);
    }
}
/**
 * @brief square all coefficient of weil thetapoint of dim 4.
 *
 * @param out (output) The squared thetapoint
 * @param in (input) The input thetapoint of dim 4.
 *
 */
static inline void
theta_squared_weil(theta_point_weil_t *out, const theta_point_weil_t *in)
{
    for (uint8_t i = 0; i < 10; i++) {
        fp_sqr(&(*out)[i], &(*in)[i]);
    }
}

/**
 * @brief perform dot product of 2 thetapoints of dim 4.
 *
 * @param out (output) The product thetapoint.
 * @param in1 (input) The first thetapoint of dim 4.
 * @param in2 (input) The second thetapoint of dim 4.
 *
 */
static inline void
theta_dot_prod_dim4(theta_point_dim4_t *out, const theta_point_dim4_t *in1, const theta_point_dim4_t *in2)
{
    for (uint8_t i = 0; i < 16; i++) {
        fp_mul(&(*out)[i], &(*in1)[i], &(*in2)[i]);
    }
}

/**
 * @brief perform dot product of 1 thetapoints of dim 4 and a weil point.
 *
 * @param out (output) The product thetapoint.
 * @param in1 (input) The thetapoint of dim 4.
 * @param in2 (input) The weil thetapoint of dim 4.
 *
 */
static inline void
theta_dot_prod_weil(theta_point_dim4_t *out, const theta_point_dim4_t *in1, const theta_point_weil_t *in2)
{
    for (uint8_t i = 0; i < 16; i++) {
        fp_mul(&(*out)[i], &(*in1)[i], &(*in2)[theta_to_weil[i]]);
    }
}

/**
 * @brief Copy a thetapoint to another.
 *
 * @param out (output) where to copy the thetapoint.
 * @param in (input) The thetapoint to copy.
 *
 */
static inline void
theta_copy_theta_dim4(theta_point_dim4_t *out, const theta_point_dim4_t *in)
{
    for (uint8_t i = 0; i < 16; i++) {
        fp_copy(&(*out)[i], &(*in)[i]);
    }
}

/**
 * @brief Copy a theta weil point to another.
 *
 * @param out (output) where to copy the thetapoint.
 * @param in (input) The thetapoint to copy.
 *
 */
static inline void
theta_copy_weil(theta_point_weil_t *out, const theta_point_weil_t *in)
{
    for (uint8_t i = 0; i < 10; i++) {
        fp_copy(&(*out)[i], &(*in)[i]);
    }
}

/**
 * @brief invert a theta point.
 *
 * @param out (output) where to write the inverted thetapoint.
 * @param in (input) The thetapoint to invert.
 *
 */
static inline void
theta_invert_theta_dim4(theta_point_dim4_t *out, const theta_point_dim4_t *in)
{
    fp_t inv[16];
    for (uint8_t i = 0; i < 16; i++) {
        fp_copy(&inv[i], &(*in)[i]);
    }
    fp_proj_batched_inv(inv, 16);
    for (uint8_t i = 0; i < 16; i++) {
        fp_copy(&(*out)[i], &inv[i]);
    }
}

/**
 * @brief invert a theta weil point.
 *
 * @param out (output) where to write the inverted thetapoint.
 * @param in (input) The thetapoint to invert.
 *
 */
static inline void
theta_invert_weil(theta_point_weil_t *out, const theta_point_weil_t *in)
{
    fp_t inv[10];
    for (uint8_t i = 0; i < 10; i++) {
        fp_copy(&inv[i], &(*in)[i]);
    }
    fp_proj_batched_inv(inv, 10);
    for (uint8_t i = 0; i < 10; i++) {
        fp_copy(&(*out)[i], &inv[i]);
    }
}

/**
 * @brief Copy a weil-theta-structure.
 *
 * @param out (output) where to copy the theta-structure.
 * @param in (input)  the theta-structure to copy
 *
 */
static inline void
theta_copy_weil_struct(theta_struct_weil_t *out, const theta_struct_weil_t *in)
{
    theta_copy_weil(&out->inv_dual_null_point, &in->inv_dual_null_point);
}

/**
 * @brief Copy a weil theta-structure precomp
 *
 * @param out (output) where to copy the weil theta-structure precomp.
 * @param in (input)  the weil theta-structure precomp to copy
 *
 */
static inline void
theta_copy_weil_struct_precomp(theta_struct_weil_precomp_t *out, const theta_struct_weil_precomp_t *in)
{
    theta_copy_weil(&out->inv_null_point_DBL, &in->inv_null_point_DBL);
    theta_copy_weil(&out->inv_HS_null_point_DBL, &in->inv_HS_null_point_DBL);
}

/**
 * @brief Perform the right precomputation to DBL points.
 *
 * @param[out] theta_precomp the theta-structure precomputation
 * @param[in] theta_struct the theta-structure where we need to do the precomputation
 *
 */
void theta_struct_weil_arith_precomp(theta_struct_weil_precomp_t *theta_precomp,
                                     const theta_struct_weil_t *theta_struct);

/** // end utility functions for dim 4 theta structures
 * @}
 */

/** @defgroup Algebra Algebra
 * @{
 */

/** @brief double a point on an an Weil type abelian variety.
 * @note Unless precomputed, weil theta structure are in dual mode.
 * @param[out] out: the doubled point
 * @param[out] in: the point to double
 * @param[out] theta_struct: the precomputation of the underlying theta structure.
 */
void theta_DBL_weil(theta_point_dim4_t *out,
                    const theta_point_dim4_t *in,
                    const theta_struct_weil_precomp_t *theta_struct);

/** @brief Compute [2^n]P for P a point on an Weil type abelian variety.
 * @note Unless precomputed, weil theta structure are in dual mode.
 * @param[out] out: [2^n]P
 * @param[out] in: P
 * @param[out] theta_struct: the precomputation of the theta structure.
 * @param[out] n: the number of itarations
 */
void theta_DBL_iter_weil(theta_point_dim4_t *out,
                         const theta_point_dim4_t *in,
                         const theta_struct_weil_precomp_t *theta_struct,
                         const int n);

/** // end Algebra
 * @}
 */

/** @defgroup isogenies 
 * @{
 */

/**
 * @brief Given B = < T_3,T_4 > half of a maximal isotropic subgroup of A[8] expressed in a consistent theta structure,
 * compute the codomain of the isogeny f: A --> B with [4]T_2,[4]T_3 in ker(f) and A,B Weil type isogenies.
 *
 * @param[out] codomain (output) The codomain of the isogeny f: A --> B of kernel [4]B.
 * @param[in] ker (input) {T_0,...,T_3} a maximal isotropic subgroup of A[8] expressed in a consistent theta structure.
 * @param[in] dual_domain (input) indicate wether the input point is in dual mode or not.
 *
 */
void isogeny_weil_compute_codomain(theta_struct_weil_t *codomain,
                                   const theta_point_dim4_t *ker,
                                   const bool dual_domain);

/**
 * @brief Given P in consistent theta coordinates and the codomain of f, compute f(P).
 *
 * @param imP (output) The image of f(P).
 * @param P (input) a point P in theta coordinates consistent with f.
 * @param codomain (input) the codomain of f.
 * @param dual_domain (input) indicate wether the input point is in dual mode or not.
 *
 */
void isogeny_weil_eval(theta_point_dim4_t *imP,
                       const theta_point_dim4_t *P,
                       const theta_struct_weil_t *codomain,
                       const bool dual_domain);

/**
 * @brief Given T3, T4 points defining the kernel of a chain between Weil type isogeny, compute the codomain
 *
 * @param[inout] codomain input the domain of an isogeny chain between Weil type isogenies and return the codomain
 * @param[in] codomain_precomp The input domain precomputation of the isogeny chain between Weil type isogenies
 * @param[in] T3 the 3rd point of the kernel detailed in Theorem 4.2. Point of order len + 2
 * @param[in] T4 the 4th point of the kernel detailed in Theorem 4.2. Point of order len + 2
 * @param[in] len lenth of the chain
 *
 */
void isogeny_chain_weil(theta_struct_weil_t *codomain,
                        theta_struct_weil_precomp_t *codomain_precomp,
                        const theta_point_dim4_t *T3,
                        const theta_point_dim4_t *T4,
                        const unsigned int len);

/** // end isogeny functions
 * @}
 */

/** // end dim 4
 * @}
 */

/** @defgroup  dim2  Theta points dim 2
 * @{
 */

/** @defgroup  Theta-points 
 * @{
 */

/**
 * @brief Hadamard transform on a vector of 4 elements.
 *
 * @param out (output) The hadamard transformed vector.
 * @param in (input) The input vector of 4 elements.
 *
 */
static inline void
theta_hadamard_dim2(theta_point_dim2_t *out, const theta_point_dim2_t *in)
{
    fp_t tmp[4];

    fp_add(&tmp[0], &(*in)[0], &(*in)[1]); 
    fp_sub(&tmp[1], &(*in)[0], &(*in)[1]);
    fp_add(&tmp[2], &(*in)[2], &(*in)[3]);
    fp_sub(&tmp[3], &(*in)[2], &(*in)[3]);

    fp_add(&(*out)[0], &tmp[0], &tmp[2]);
    fp_add(&(*out)[1], &tmp[1], &tmp[3]);
    fp_sub(&(*out)[2], &tmp[0], &tmp[2]);
    fp_sub(&(*out)[3], &tmp[1], &tmp[3]);
}
/**
 * @brief square all coefficient of thetapoint of dim 2.
 *
 * @param out (output) The squared thetapoint
 * @param in (input) The input thetapoint of dim 2.
 *
 */
static inline void
theta_squared_dim2(theta_point_dim2_t *out, const theta_point_dim2_t *in)
{
    for (uint8_t i = 0; i < 4; i++) {
        fp_sqr(&(*out)[i], &(*in)[i]);
    }
}

/**
 * @brief perform dot product of 2 thetapoints of dim 2.
 *
 * @param out (output) The product thetapoint.
 * @param in1 (input) The first thetapoint of dim 2.
 * @param in2 (input) The second thetapoint of dim 2.
 *
 */
static inline void
theta_dot_prod_dim2(theta_point_dim2_t *out, const theta_point_dim2_t *in1, const theta_point_dim2_t *in2)
{
    for (uint8_t i = 0; i < 4; i++) {
        fp_mul(&(*out)[i], &(*in1)[i], &(*in2)[i]);
    }
}

/**
 * @brief Copy a thetapoint to another.
 *
 * @param out (output) where to copy the thetapoint.
 * @param in (input) The thetapoint to copy.
 *
 */
static inline void
theta_copy_dim2(theta_point_dim2_t *out, const theta_point_dim2_t *in)
{
    for (uint8_t i = 0; i < 4; i++) {
        fp_copy(&(*out)[i], &(*in)[i]);
    }
}

/**
 * @brief invert a theta weil point.
 *
 * @param out (output) where to write the inverted thetapoint.
 * @param in (input) The thetapoint to invert.
 *
 */
static inline void
theta_invert_dim2(theta_point_dim2_t *out, const theta_point_dim2_t *in)
{
    fp_t inv[4];
    for (uint8_t i = 0; i < 4; i++) {
        fp_copy(&inv[i], &(*in)[i]);
    }
    fp_proj_batched_inv(inv, 4);
    for (uint8_t i = 0; i < 4; i++) {
        fp_copy(&(*out)[i], &inv[i]);
    }
}

/**
 * @brief Copy a dim2 theta-structure.
 *
 * @param out (output) where to copy the theta-structure.
 * @param in (input)  the theta-structure to copy
 *
 */
static inline void
copy_theta_struct_dim2(theta_struct_dim2_t *out, const theta_struct_dim2_t *in)
{
    theta_copy_dim2(&out->inv_dual_null_point, &(in->inv_dual_null_point));
    theta_copy_dim2(&out->dual_null_point, &in->dual_null_point);
}

/**
 * @brief Copy a dim2 theta-structure precomp.
 *
 * @param out (output) where to copy the theta-structure precomp.
 * @param in (input)  the theta-structure precomp to copy
 *
 */
static inline void
copy_theta_struct_precomp_dim2(theta_struct_dim2_precomp_t *out, const theta_struct_dim2_precomp_t *in)
{
    theta_copy_dim2(&out->inv_null_point_DBL, &(in->inv_null_point_DBL));
    theta_copy_dim2(&out->inv_HS_null_point_DBL, &in->inv_HS_null_point_DBL);
}

/**
 * @brief Perform the right precomputation to double points.
 *
 * @param[out] theta_precomp the theta-structure precomputation
 * @param[in]  theta_struct the theta-structure where we need to compute the precomputation
 *
 */
void theta_struct_dim2_arith_precomp(theta_struct_dim2_precomp_t *theta_precomp,
                                     const theta_struct_dim2_t *theta_struct);

/** // end 
 * @}
 */


/** @defgroup Algebra Algebra dim 2
 * @{
 */

/** @brief double a point on an abelian variety of dim 2.
 * @note dim2 theta structure are in dual mode.
 * @param[out] out: the doubled point
 * @param[out] in: the point to double
 * @param[out] theta_struct: the underlying theta structure. Note that points are in dual mode
 */
void theta_dim2_DBL(theta_point_dim2_t *out,
                    const theta_point_dim2_t *in,
                    const theta_struct_dim2_precomp_t *theta_struct);

/** @brief Compute [2^n]P for P a point on an abelian variety of dim 2.
 * @note dim2 theta structure are in dual mode.
 * @param[out] out: [2^n]P
 * @param[out] in: P
 * @param[out] theta_struct: the underlying theta structure. Note that points are in dual mode
 * @param[out] n: the number of itarations
 */
void theta_dim2_DBL_iter(theta_point_dim2_t *out,
                         const theta_point_dim2_t *in,
                         const theta_struct_dim2_precomp_t *theta_struct,
                         const int n);

/** @brief Compute P + Q given P, Q and P-Q on an abelian variety of dim 2.
 * @note dim2 theta structure are in dual mode.
 * @param[out] out: P + Q
 * @param[out] P:  P
 * @param[out] Q:  Q
 * @param[out] inv_PmQ:  the projective inverse of P - Q
 * @param[out] theta_struct: the underlying theta structure. Note that points are in dual mode
 */
void theta_dim2_diff_ADD(theta_point_dim2_t *out,
                         const theta_point_dim2_t *P,
                         const theta_point_dim2_t *Q,
                         const theta_point_dim2_t *inv_PmQ,
                         const theta_struct_dim2_precomp_t *theta_struct);

/** // end Algebra
 * @}
 */

/** @defgroup isogenies 
 * @{
 */

/**
 * @brief Given B = < T_1,T_2 > basis of A[8] expressed in a consistent theta structure, compute the codomain of the
 * isogeny f: A --> B
 *
 * @param codomain (output) The codomain of the isogeny f: A --> B of kernel [4]B.
 * @param ker (input) {T_0,...,T_3} a maximal isotropic subgroup of A[8] expressed in a consistent theta structure.
 * NOTE: input in normal mode.
 */
void isogeny_dim2_compute_codomain(theta_struct_dim2_t *codomain, const theta_point_dim2_t *ker);

/**
 * @brief Given P in consistent theta coordinates and the codomain of f, compute f(P).
 *
 * @param imP (output) The image of f(P).
 * @param P (input) a point P in theta coordinates consistent with f.
 * @param codomain (input) the codomain of f.
 * NOTE: input in normal mode, and outputs in dual form, because this is only used in gluing and we need maximum control
 * because of dual stuff.
 */
void isogeny_dim2_eval(theta_point_dim2_t *imP, const theta_point_dim2_t *P, const theta_struct_dim2_t *codomain);

/** // end isogeny functions
 * @}
 */


#endif
