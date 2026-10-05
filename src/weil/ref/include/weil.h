#ifndef WEIL_H
#define WEIL_H

#include <theta.h>
#include <fp.h>
#include <fp2.h>
#include <ec.h>
#include <ec_isog.h>
#include <ec_params.h>

extern const uint8_t gluing_change_struct_dim4[16]; // = { 0, 5, 10, 15, 1, 4, 11, 14, 8, 13, 2, 7, 9, 12, 3, 6 };

/** @defgroup struct MIKE structures
 * @{
 */

/** @brief MIKE gluing basis
 *
 * @typedef mike_gluing_basis_t
 *
 * @struct mike_gluing_basis_t
 */
typedef struct mike_gluing_basis_t
{
    ec_jac_point_t P_8;
    ec_jac_point_t Q_8;

    ec_jac_point_t P_16;
    ec_jac_point_t Q_16;

    ec_jac_point_t T3_1_32;
    ec_jac_point_t T3_2_32;
    ec_jac_point_t T4_1_32;
    ec_jac_point_t T4_2_32;

} mike_gluing_basis_t;


/** @brief MIKE isogeny chain input data
 *
 * @typedef mike_chain_basis_t
 *
 * @struct mike_chain_basis_t
 */
typedef struct mike_chain_basis_t
{

    // Point T3 (coordinates on E only)
    ec_jac_point_t T3_1;
    ec_jac_point_t T3_2;

    // Point T4 (coordinates on E only)
    ec_jac_point_t T4_1;
    ec_jac_point_t T4_2;
} mike_chain_basis_t;

// Structure for Scholten gluing

/** @brief Coordinates used to construct the theta structure on the Weil restriction compatible with Scholten gluing
 *
 * @typedef weil_restriction_theta_basis_t
 *
 * @struct weil_restriction_theta_basis_t
 */
typedef fp_t weil_restriction_theta_basis_t[4];

/** @brief Precomputation of the Scholten gluing
 *
 * @typedef scholten_precomp_t
 *
 * @struct scholten_precomp_t
 */
typedef struct scholten_precomp_t
{
    fp_t A;
    fp_t B;
    fp_t C;
    fp_t D;
    fp_t E;

} scholten_precomp_t;

/** @brief  Scholten gluing
 *
 * @typedef scholten_gluing_t
 *
 * @struct scholten_gluing_t
 */
typedef struct scholten_gluing_t
{
    scholten_precomp_t precomp;
    theta_point_dim2_t inv_image_aux_point;
    ec_jac_point_t aux_point; //= P_8
    ec_curve_t E;

} scholten_gluing_t;

/** @brief  Weil gluing
 *
 * @typedef weil_gluing_t
 *
 * @struct weil_gluing_t
 */
typedef struct diag_gluing_t
{
    theta_struct_dim2_t S_U;
    theta_struct_dim2_precomp_t S_U_precomp;
    theta_struct_dim2_t S_V;
    theta_struct_dim2_precomp_t S_V_precomp;

} diag_gluing_t;

/** @brief  Weil gluing
 *
 * @typedef weil_gluing_t
 *
 * @struct weil_gluing_t
 */
typedef struct weil_gluing_t
{
    theta_struct_weil_t weil_codom;
    theta_struct_weil_precomp_t weil_precomp;
    theta_point_dim4_t inv_aux_point;

} weil_gluing_t;

/** @brief  MIKE gluing
 *
 * @typedef mike_gluing_t
 *
 * @struct mike_gluing_t
 */
typedef struct mike_gluing_t
{
    scholten_gluing_t scholten_gluing;
    diag_gluing_t phi_UV;
    weil_gluing_t weil_gluing;
    ec_jac_point_t aux_point_1; //= T3_1_32
    ec_jac_point_t aux_point_2; //= T3_2_32

} mike_gluing_t;

/** @brief MIKE relative invariants (I_1/I_0,...,I_9/I_0) from Definition 4.7
 *
 * @typedef mike_rel_invariants_t
 *
 * @struct mike_rel_invariants_t
 */
typedef fp_t mike_rel_invariants_t[9];

/** @brief MIKE relative invariants (J_2,...,J_5) from Definition 4.10
 *
 * @typedef mike_abs_invariants_t
 *
 * @struct mike_abs_invariants_t
 */
typedef fp_t mike_abs_invariants_t[4];

/** // end struct
 * @}
 */

/** @defgroup Gluing
 * @{
 */

/** @defgroup Scholten Scholten gluing
 * @{
 */

/**
 * @brief Given P a non canonical 4 torsion point over E, compute m the coefficients defining the theta coordinates over
 * the Weil restriction E x E^p, as described in Lemma F.10 of the paper.
 *
 * @param[out] coeff_basis pointer to [m0,m1,m2,m3] defining the change of theta structure on the Weil restriction;
 * @param[in] P A 4 torsion point different from (1: ± 1).
 */
void scholten_variable_compute(weil_restriction_theta_basis_t *coeff_basis, const ec_xz_point_t *P);

/**
 * @brief lift a point P, representing (P, pi(P)) in the theta model of the Weil restriction E x E^p
 *
 * @param[out] out representation of the point in the theta model;
 * @param[in] coeff_basis pointer to [m0,m1,m2,m3] defining the change of theta structure on the Weil restriction;
 * @param[in] P a point representing (P, pi(P))
 */
void scholten_jacobian_to_theta(theta_point_dim2_t *out,
                                const weil_restriction_theta_basis_t *coeff_basis,
                                const ec_xz_point_t *P);

/**
 * @brief Computes the 5 precomputation coefficients needed to compute the Scholten gluing using Theorem F.11
 *
 * @param[out] scholten the 5 coordinates used to for fast Scholten gluing
 * @param[in] coeff_basis pointer to [m0,m1,m2,m3] defining the change of theta structure on the Weil restriction;
 */
void scholten_precompute(scholten_precomp_t *scholten, const weil_restriction_theta_basis_t *coeff_basis);

/**
 * @brief Computes the formula of Theorem F.11. Given (u:v:w) = (P + Q)(P - Q), compute F1(P)xF1(Q) on the Scholten.
 *
 * @param[out] out F1(P)xF1(Q) on Scholten in the dual theta model.
 * @param[in] uvw barycentric coordinates of (P + Q)(P - Q)
 * @param[in] scholten the 5 coordinates used to for fast Scholten gluing
 * @param[in] twist Indicating wether we are evaluating points on the twist or not [0] --> (P, P^p) and [1] --> (P,
 * -P^p).
 */
void scholten_eval(theta_point_dim2_t *out,
                   const ec_bary_coords_t *uvw,
                   const scholten_precomp_t *scholten,
                   bool twist);

/**
 * @brief Perform all computation to perform the Scholten Gluing.
 *
 * @param[out] phi the Scholten gluing
 * @param[in] T1_8 first point of 8 torsion above ker(phi)
 * @param[in] T2_8 second point of 8 torsion above ker(phi) such that [4]T2_8 = (0:0:1)
 * @param[in] E the elliptic curve over which T1_8 and T2_8 are defined.
 *
 * @returns 0xFFFFFFFF if the inputs are valid, 0 otherwise.
 */
uint32_t schoten_gluing_compute_and_verify(scholten_gluing_t *phi,
                                           const ec_jac_point_t *T1_8,
                                           const ec_jac_point_t *T2_8,
                                           const ec_curve_t *E);

/**
 * @brief Given P, evaluate F1(P).
 *
 * @note formulae do not work if P = P_8.
 *
 * @param[out] out F1(P)
 * @param[in] P Point we want to compute
 * @param[in] phi Scholten gluing
 * @param[in] twist Indicating wether we are evaluating points on the twist or not [0] --> (P, P^p) and [1] --> (P,
 * -P^p).
 *
 */
void scholten_gluing_eval(theta_point_dim2_t *out, const ec_jac_point_t *P, const scholten_gluing_t *phi, bool twist);

// end scholten
/** @}
 */

/** @defgroup diag Diagonal isogeny
 * @{
 */

/**
 * @brief Given a point P on the Scholten in the dual theta model compatible with the Scholten gluing, apply the right
 * change of theta coordinates so that the. point end up in tha theta structure compatible with the evalaution of phi_V
 *
 * @param[out] out output
 * @param[in] in input
 *
 */
void gluing_change_theta_structure_dim2_to_V_compatible_theta_struct(theta_point_dim2_t *out,
                                                                     const theta_point_dim2_t *in);

/**
 * @brief Given a point P on the Scholten in the dual theta model compatible with the Scholten gluing, apply the right
 * change of theta coordinates so that the. point end up in tha theta structure compatible with the evalaution of phi_U
 *
 * @param[out] out output
 * @param[in] in input
 *
 */
void gluing_change_theta_structure_dim2_to_U_compatible_theta_struct(theta_point_dim2_t *out,
                                                                     const theta_point_dim2_t *in);

/**
 * @brief Given V1, V2 on the Scholten, 8 torsion points, compute the diagonal isogeny phi_UV = Diag(phi_U, phi_V).
 *
 * @param[out] phiUV the diagonal isogeny
 * @param[in] V1_8 8 torsion points above ker(phi_V) in the dual theta model compatible with the Scholten gluing.
 * @param[in] V2_8 8 torsion points above ker(phi_V) in the dual theta model compatible with the Scholten gluing.
 */
void gluing_diag_isogeny_compute(diag_gluing_t *phiUV, const theta_point_dim2_t *V1_8, const theta_point_dim2_t *V2_8);

/**
 * @brief in1, in2 a couple of theta points in the dual theta model compatible with Scholten gluing, compute phi_U(in1)
 * x phi_V(in2). Outputs are in their respective dual form
 *
 * @param[out] out1 phi_U(in1) in the dual model
 * @param[out] out2 phi_V(in2) in the dual model
 * @param[in] in1 Point in the dual theta model compatible with the Scholten gluing.
 * @param[in] in2 Point in the dual theta model compatible with the Scholten gluing.
 * @param[in] phiUV diagonal isogeny
 */
void gluing_diag_isogeny_eval(theta_point_dim2_t *out1,
                              theta_point_dim2_t *out2,
                              const theta_point_dim2_t *in1,
                              const theta_point_dim2_t *in2,
                              const diag_gluing_t *phiUV);

// end diagonal gluing
/** @}
 */

/** @defgroup Weil-gluing Weil gluing
 * @{
 */

//

/**
 * @brief Take 2 dim 2 theta points in dual mode and apply the right change of variables to get points compatible with
// the chain in dim4.
 *
 * @param[out] out tensor product of in1 x in2 compatible with the rest of the chain.
 * @param[in] in1 Point in the dual theta model compatible with phi_U.
 * @param[in] in2 Point in the dual theta model compatible with phi_V.
 */
void gluing_from_couple_dim2_to_dim4_compatible_with_isogeny(theta_point_dim4_t *out,
                                                             theta_point_dim2_t *in1,
                                                             theta_point_dim2_t *in2);

/**
 * @brief compute the projective inverse of the gluing codomain (works around the zeros)
 *
 * @param[out] inv_vec projective inverse
 * @param[in] vec input codomain
 */
void gluing_special_inv(theta_point_weil_t *inv_vec, const theta_point_weil_t *vec);

/**
 * @brief Compute the Weil gluing, as described in Algorithm 7.
 *
 * @param[out] weil_gluing gluing isogeny
 * @param[in] T3 Point above the weil gluing kernel, whose theta action is (0, 4).
 * @param[in] T4 Point above the weil gluing kernel, whose theta action is (0, 8).
 * @param[in] T34 sum of T3 + T4
 * @param[in] T3_p_2T4  T3 + [2]T4
 *
 */
void gluing_compute_weil_gluing(weil_gluing_t *weil_gluing,
                                const theta_point_dim4_t *T3,
                                const theta_point_dim4_t *T4,
                                const theta_point_dim4_t *T34,
                                const theta_point_dim4_t *T3_p_2T4);

/**
 * @brief Evalaute through the  Weil gluing.
 *
 * @param[out] out gluing isogeny
 * @param[in] P_p_T3 P + T3 in a compatible theta structure.
 * @param[in] P_m_T3 P - T3 in a compatible theta structure.
 * @param[in] gluing sum of T3 + T4
 *
 */
void gluing_evaluate_weil_gluing(theta_point_dim4_t *out,
                                 const theta_point_dim4_t *P_p_T3,
                                 const theta_point_dim4_t *P_m_T3,
                                 const weil_gluing_t *gluing);

// end Weil-gluing
/** @}
 */

/** @defgroup Global Global gluing
 * @{
 */

/**
 * @brief computes the gluing basis from PQ32;
 *
 * @param[out] basis pointer to a structure encapsulating the gluing kernel basis
 * @param[in] E pointer to the public key curve (supersingular and defined over Fp2 strictly)
 * @param[in] PQ_32 basis of E[2^{e+2}] in Jacobian coordinates (table of length 2)
 * @param[in] x_bshift_mod32 integer congruent to the secret key bitshifted by 1 mod 32 (i.e. (x_sec - 1) >> 1 )
 *
 */
void gluing_basis_compute(mike_gluing_basis_t *basis,
                          const ec_curve_t *E,
                          const ec_jac_point_t *PQ_32,
                          const uint8_t x_bshift_mod32);

/**
 * @brief Computes the gluing F_glue=F3*F2*F1 (Weil*Diagonal*Scholten) as expressed
 * in Theorem F.9.
 *
 * @param[out] gluing pointer to a structure encapsulating the gluing data needed for evaluation
 * (see gluing_mike_eval)
 * @param[in] E pointer to the public key curve (supersingular and defined over Fp2 strictly)
 * @param[in] basis pointer to a structure encapsulating the gluing kernel basis
 */
uint32_t gluing_mike_compute(mike_gluing_t *gluing, const ec_curve_t *E, const mike_gluing_basis_t *basis);

/**
 * @brief Evaluates the gluing F_glue=F3*F2*F1 as outputed by gluing_mike_compute at
 * points of the form (T0,-sigma(T0),T1,-sigma(T1)) in W(E)^2.
 *
 * @param[out] imgT pointer to the gluing image point F_glue(T)
 * @param[in] T0 pointer to an elliptic curve point in Jacobian coordinates
 * @param[in] T1 pointer to an elliptic curve point in Jacobian coordinates
 * @param[in] gluing pointer to a structure encapsulating the gluing data
 */
void gluing_mike_eval(theta_point_dim4_t *imgT,
                      const ec_jac_point_t *T0,
                      const ec_jac_point_t *T1,
                      const mike_gluing_t *gluing);

// end Global
/** @}
 */

// end Gluing
/** @}
 */

/** @defgroup invariants
 * @{
 */

/**
 * @brief Given a MIKE adapted theta null point on an abelian fourfold of Weil type,
 * computes its absolute MIKE invariants (given by Definition 4.10). This is Algoritm 3.
 *
 * @param[out] abs_inv pointer to the absolute MIKE invariants (J2,...,J5)
 * @param[in] theta pointer to a theta null point of Weil type adapted to MIKE isogenies
 *
 * Total cost: 1I + 28M + 22S + 64a
 */
void mike_absolute_invariants(mike_abs_invariants_t *abs_inv, const theta_point_weil_t *theta);

// end Invariants
/** @}
 */

/** @defgroup Chain MIKE 4-dimensional chain
 * @{
 */
/**
 * @brief Computes the T3, T4 points of the basis of isogeny corresponding to Theorem 4.2.
 *
 * @param[out] basis: the T3, T4 points. They are points of
 * @param[in] P: A point of 2^pow torsion not above (0:0:1).
 * @param[in] Q: A point of 2^pow torsion above (0:0:1).
 * @param[in] x_div2: the secret coordinates. It is used shifted right by 1 bit
 * @param[in] pow: the 2-power of the torsion points.
 * @param[in] E: the curve above which P,Q are defined over.
 *
 */
void mike_compute_basis(mike_chain_basis_t *basis,
                        const ec_jac_point_t *P,
                        const ec_jac_point_t *Q,
                        const digit_t *x_div2,
                        const int pow,
                        const ec_curve_t *E);

/**
 * @brief Computes MIKE shared secret key from input 4-dimensional isogeny chain data
 * @param[out] mike_abs_inv: the absolute invariant (J2,..., J5) of the codomain of the isogeny.
 * @param[in] mike_basis: the T3 and T4 points of the basis of the 4D chain. They are points of 2^(length + 2) torsion
 * @param[in] gluing_basis: the basis of the gluing
 * @param[in] E: the codomain curve
 * @param[in] length: the length of the chain
 *
 * @returns 0xFFFFFFFF if the inputs are valid, 0 otherwise.
 */
uint32_t mike_isogeny_chain_dim4(mike_abs_invariants_t *mike_abs_inv,
                                 const mike_chain_basis_t *mike_basis,
                                 const mike_gluing_basis_t *gluing_basis,
                                 const ec_curve_t *E,
                                 const unsigned int length);

// end Chain
/** @}
 */

#endif