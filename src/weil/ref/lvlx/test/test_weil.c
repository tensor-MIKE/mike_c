#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <rng.h>
#include <inttypes.h>
#include <bench_test_arguments.h>

#include <weil.h>

// compute scalar product mod 2
static inline uint8_t
scalprod(uint8_t i, uint8_t j)
{
    uint8_t ret = 0;
    ret ^= (i & 1) & (j & 1);
    ret ^= ((i >> 1) & 1) & ((j >> 1) & 1);
    ret ^= ((i >> 2) & 1) & ((j >> 2) & 1);
    ret ^= ((i >> 3) & 1) & ((j >> 3) & 1);
    return ret;
}

// helper function used to sample a supersingular curve over Fp2
static void
test_sample_random_supersingular_curve(ec_curve_t *E)
{

    ec_curve_t E0;
    ec_jac_point_t jacP, jacQ, jacPQ;
    ec_xz_point_t P, Q, PQ;
    ec_curve_init_precomputed(&E0);

    ec_jac_basis_2f(&jacP, &jacQ, &E0, TORSION_EVEN_POWER);

    ec_jac_neg(&jacQ, &jacQ);
    ec_ADD(&jacPQ, &jacP, &jacQ, &E0);
    ec_jac_neg(&jacQ, &jacQ);

    ec_jac_to_xz(&P, &jacP);
    ec_jac_to_xz(&Q, &jacQ);
    ec_jac_to_xz(&PQ, &jacPQ);

    digit_t num[NWORDS_ORDER] = { 0 };
    for (uint8_t i = 0; i < NWORDS_ORDER; i++) {
        num[i] = 0;
    }

    ec_xz_point_t ker;

    uint8_t rand[4];
    randombytes(rand, 4);
    num[0] = rand[0] ^ (rand[1] << 8) ^ (rand[2] << 16) ^ (rand[3] << 24);

    if (!ec_ladder3pt(&ker, &P, &Q, &PQ, num, ec_n_bits(num, NWORDS_ORDER), &E0)) {
        fprintf(stderr, "ec_ladder3pt failed in test_sample_random_supersingular_curve\n");
        exit(1);
    }

    ec_isog_even_t phi;
    ec_copy_curve(&phi.curve, &E0);
    ec_copy_xz_point(&phi.kernel, &ker);
    phi.length = TORSION_EVEN_POWER;

    fp2_t A;
    ec_iso_isogeny_2chain(&A, &phi, 2);
    ec_curve_init_from_A(E, &A);
}

int
test_theta_structure_scholten(int iterations)
{

    for (int g = 0; g < iterations; g++) {
        ec_curve_t E;
        ec_jac_point_t jacP, jacQ;

        test_sample_random_supersingular_curve(&E);
        ec_jac_basis_2f(&jacP, &jacQ, &E, TORSION_EVEN_POWER);

        // Start true test
        // ec_jac_point_t jacP8, jacQ8;
        ec_jac_point_t jacP4, jacQ4, jacPQ4;
        ec_xz_point_t xz_P4, xz_Q4, xz_PQ4;

        ec_jac_dbl_iter(&jacP4, &jacP, TORSION_EVEN_POWER - 2, &E, 1);
        ec_jac_dbl_iter(&jacQ4, &jacQ, TORSION_EVEN_POWER - 2, &E, 1);

        ec_ADD(&jacPQ4, &jacP4, &jacQ4, &E);

        ec_jac_to_xz(&xz_P4, &jacP4);
        ec_jac_to_xz(&xz_Q4, &jacQ4);
        ec_jac_to_xz(&xz_PQ4, &jacPQ4);

        // testing change of theta coordinates
        weil_restriction_theta_basis_t m;
        theta_point_dim2_t theta_P4, theta_Q4, theta_PQ4;

        // compute value necessary for theta struct on Weil restriction
        scholten_variable_compute(&m, &xz_P4);

        // lift 4-torsion points
        scholten_jacobian_to_theta(&theta_P4, &m, &xz_P4);
        scholten_jacobian_to_theta(&theta_Q4, &m, &xz_Q4);
        scholten_jacobian_to_theta(&theta_PQ4, &m, &xz_PQ4);

        if (!fp_is_zero(&theta_P4[1]) || !fp_is_zero(&theta_P4[3])) {
            printf("theta structure is not compatible with P4 points\n");
            return 1;
        }

        if (!fp_is_zero(&theta_Q4[2]) || !fp_is_zero(&theta_Q4[3])) {
            printf("theta structure is not compatible with Q4 points\n");
            return 1;
        }

        if (!fp_is_zero(&theta_PQ4[1]) || !fp_is_zero(&theta_PQ4[2])) {
            printf("theta structure is not compatible with P4 + Q4 points\n");
            return 1;
        }
    }

    printf("computation theta structure Weil restriction...........................PASSED\n");

    return 0;
}

int
test_scholten_formulae(int iterations)
{

    for (int g = 0; g < iterations; g++) {

        ec_curve_t E;
        ec_jac_point_t jacP, jacQ;

        test_sample_random_supersingular_curve(&E);
        ec_jac_basis_2f(&jacP, &jacQ, &E, TORSION_EVEN_POWER);

        // Start true test
        ec_jac_point_t jacP8, jacQ8;
        ec_jac_point_t jacP4;
        ec_xz_point_t xz_P8, xz_Q8, xz_P4;

        ec_jac_dbl_iter(&jacP8, &jacP, TORSION_EVEN_POWER - 3, &E, 1);
        ec_jac_dbl_iter(&jacQ8, &jacQ, TORSION_EVEN_POWER - 3, &E, 1);

        ec_DBL(&jacP4, &jacP8, &E);

        ec_jac_to_xz(&xz_P4, &jacP4);
        ec_jac_to_xz(&xz_P8, &jacP8);
        ec_jac_to_xz(&xz_Q8, &jacQ8);

        // testing change of theta coordinates
        theta_point_dim2_t theta_P8_dual, theta_Q8_dual;
        theta_point_dim2_t theta_P8_super, theta_Q8_super;

        // Do precomputation
        weil_restriction_theta_basis_t Weil_basis;
        // compute value necessary for theta struct on Weil restriction
        scholten_variable_compute(&Weil_basis, &xz_P4);
        // Superglue formula
        scholten_precomp_t precomp;
        scholten_precompute(&precomp, &Weil_basis);

        // lift 8-torsion points
        scholten_jacobian_to_theta(&theta_P8_dual, &Weil_basis, &xz_P8);
        scholten_jacobian_to_theta(&theta_Q8_dual, &Weil_basis, &xz_Q8);

        theta_squared_dim2(&theta_P8_dual, &theta_P8_dual);
        fp_neg(&theta_P8_dual[3],
               &theta_P8_dual[3]); // last coefficient was complex, so once squared, need to sign-flip
        theta_hadamard_dim2(&theta_P8_dual, &theta_P8_dual);
        theta_squared_dim2(&theta_Q8_dual, &theta_Q8_dual);
        fp_neg(&theta_Q8_dual[3],
               &theta_Q8_dual[3]); // last coefficient was complex, so once squared, need to sign-flip
        theta_hadamard_dim2(&theta_Q8_dual, &theta_Q8_dual);
        ec_bary_coords_t uvw;

        // fast eval P
        fp2_copy(&uvw.u, &xz_P8.x);
        fp2_set_zero(&uvw.v);
        fp2_copy(&uvw.w, &xz_P8.z);

        scholten_eval(&theta_P8_super, &uvw, &precomp, false);
        // fast eval Q
        fp2_copy(&uvw.u, &xz_Q8.x);
        fp2_set_zero(&uvw.v);
        fp2_copy(&uvw.w, &xz_Q8.z);

        scholten_eval(&theta_Q8_super, &uvw, &precomp, false);

        fp_t tmp, tmp2, tmp3, tmp4;
        // Check P8
        fp_copy(&tmp, &theta_P8_dual[0]);
        fp_inv(&tmp);
        fp_copy(&tmp2, &theta_P8_super[0]);
        fp_inv(&tmp2);
        for (int i = 0; i < 4; i++) {
            fp_mul(&tmp3, &theta_P8_dual[i], &tmp);
            fp_mul(&tmp4, &theta_P8_super[i], &tmp2);
            if (!fp_is_equal(&tmp3, &tmp4)) {
                printf("P8 is not the same between canonical lift and superglue\n");
                return 1;
            }
        }
        // Check Q8
        fp_copy(&tmp, &theta_Q8_dual[0]);
        fp_inv(&tmp);
        fp_copy(&tmp2, &theta_Q8_super[0]);
        fp_inv(&tmp2);
        for (int i = 0; i < 4; i++) {
            fp_mul(&tmp3, &theta_Q8_dual[i], &tmp);
            fp_mul(&tmp4, &theta_Q8_super[i], &tmp2);
            if (!fp_is_equal(&tmp3, &tmp4)) {
                printf("Q8 is not the same between canonical lift and superglue\n");
                return 1;
            }
        }

        // test P8 + Q8
        ec_jac_point_t P_p_Q8, P_m_Q8;
        ec_xz_point_t xz_P_p_Q8, xz_P_m_Q8;
        ec_ADD(&P_p_Q8, &jacP8, &jacQ8, &E);
        ec_jac_neg(&jacQ8, &jacQ8);
        ec_ADD(&P_m_Q8, &jacP8, &jacQ8, &E);
        ec_jac_neg(&jacQ8, &jacQ8);

        ec_jac_to_xz(&xz_P_p_Q8, &P_p_Q8);
        ec_jac_to_xz(&xz_P_m_Q8, &P_m_Q8);

        scholten_jacobian_to_theta(&theta_P8_dual, &Weil_basis, &xz_P_p_Q8);
        scholten_jacobian_to_theta(&theta_Q8_dual, &Weil_basis, &xz_P_m_Q8);

        theta_point_dim2_t theta_prod_PQ8_dual, theta_prod_PQ8_super;
        theta_dot_prod_dim2(&theta_prod_PQ8_dual, &theta_P8_dual, &theta_Q8_dual);
        fp_neg(&theta_prod_PQ8_dual[3],
               &theta_prod_PQ8_dual[3]); // last coefficient was complex, so once squared, need to sign-flip
        theta_hadamard_dim2(&theta_prod_PQ8_dual, &theta_prod_PQ8_dual);

        ec_jac_to_bary_coords(&uvw, &jacP8, &jacQ8, &E);
        scholten_eval(&theta_prod_PQ8_super, &uvw, &precomp, false);

        // Check Q8
        fp_copy(&tmp, &theta_prod_PQ8_dual[0]);
        fp_inv(&tmp);
        fp_copy(&tmp2, &theta_prod_PQ8_super[0]);
        fp_inv(&tmp2);
        for (int i = 0; i < 4; i++) {
            fp_mul(&tmp3, &theta_prod_PQ8_dual[i], &tmp);
            fp_mul(&tmp4, &theta_prod_PQ8_super[i], &tmp2);

            if (!fp_is_equal(&tmp3, &tmp4)) {
                printf("P8 * Q8 is not the same between canonical lift and superglue\n");
                return 1;
            }
        }
    }

    printf("compatibility generic and superglue....................................PASSED\n");

    return 0;
}

int
test_scholten_gluing(int iterations)
{

    for (int g = 0; g < iterations; g++) {

        ec_curve_t E;
        ec_jac_point_t jacP, jacQ;

        test_sample_random_supersingular_curve(&E);
        ec_jac_basis_2f(&jacP, &jacQ, &E, TORSION_EVEN_POWER);

        // start true test //
        // Start true test
        ec_jac_point_t jacP8, jacQ8;
        ec_jac_point_t jacP4, jacQ4, jacPQ4;

        ec_jac_dbl_iter(&jacP8, &jacP, TORSION_EVEN_POWER - 3, &E, 1);
        ec_jac_dbl_iter(&jacQ8, &jacQ, TORSION_EVEN_POWER - 3, &E, 1);
        ec_DBL(&jacP4, &jacP8, &E);
        ec_DBL(&jacQ4, &jacQ8, &E);
        ec_ADD(&jacPQ4, &jacP4, &jacQ4, &E);

        scholten_gluing_t psi;
        schoten_gluing_compute_and_verify(&psi, &jacP8, &jacQ8, &E);

        theta_point_dim2_t theta_scholten_T4_1, theta_scholten_T4_2, theta_scholten_T4_1_2;
        scholten_gluing_eval(&theta_scholten_T4_1, &jacP4, &psi, false);
        scholten_gluing_eval(&theta_scholten_T4_2, &jacQ4, &psi, false);
        scholten_gluing_eval(&theta_scholten_T4_1_2, &jacPQ4, &psi, false);

        fp_t tmp, tmp2, tmp3, tmp4;
        fp_copy(&tmp, &theta_scholten_T4_1[1]);
        fp_copy(&tmp3, &theta_scholten_T4_2[2]);
        fp_inv(&tmp);
        fp_inv(&tmp3);

        for (int i = 0; i < 4; i++) {
            fp_mul(&tmp2, &theta_scholten_T4_1[i ^ 1], &tmp);
            fp_mul(&tmp4, &theta_scholten_T4_2[i ^ 2], &tmp3);
            if (!fp_is_equal(&tmp2, &tmp4)) {
                printf("Scholten codomain does not follow the theta group action.\n");
                return 1;
            }
        }
    }
    printf("Schoten gluing.........................................................PASSED\n");

    return 0;
}

int
test_diagonal_isogenies(int iterations)
{

    for (int g = 0; g < iterations; g++) {

        ec_curve_t E;
        ec_jac_point_t jacP, jacQ;

        test_sample_random_supersingular_curve(&E);
        ec_jac_basis_2f(&jacP, &jacQ, &E, TORSION_EVEN_POWER);

        // start true test //
        ec_jac_point_t jacP32, jacQ32;
        ec_jac_point_t jacP16, jacQ16;
        ec_jac_point_t jacP8, jacQ8;

        ec_jac_dbl_iter(&jacP32, &jacP, TORSION_EVEN_POWER - 5, &E, 1);
        ec_jac_dbl_iter(&jacQ32, &jacQ, TORSION_EVEN_POWER - 5, &E, 1);
        ec_DBL(&jacP16, &jacP32, &E);
        ec_DBL(&jacQ16, &jacQ32, &E);
        ec_DBL(&jacP8, &jacP16, &E);
        ec_DBL(&jacQ8, &jacQ16, &E);

        scholten_gluing_t scholten_gluing;
        schoten_gluing_compute_and_verify(&scholten_gluing, &jacP8, &jacQ8, &E);

        theta_point_dim2_t V_8_dual[2], U_8_dual[2];
        theta_point_dim2_t V_8[2];
        scholten_gluing_eval(&V_8_dual[1], &jacP16, &scholten_gluing, true); // V2
        scholten_gluing_eval(&V_8_dual[0], &jacQ16, &scholten_gluing, true); // V1

        // From V, compute U
        theta_copy_dim2(&U_8_dual[0], &V_8_dual[1]);
        fp_neg(&U_8_dual[0][3], &U_8_dual[0][3]);
        theta_copy_dim2(&U_8_dual[1], &V_8_dual[0]);
        fp_neg(&U_8_dual[1][3], &U_8_dual[1][3]);

        gluing_change_theta_structure_dim2_to_U_compatible_theta_struct(&U_8_dual[0], &U_8_dual[0]);
        gluing_change_theta_structure_dim2_to_U_compatible_theta_struct(&U_8_dual[1], &U_8_dual[1]);
        gluing_change_theta_structure_dim2_to_V_compatible_theta_struct(&V_8[0], &V_8_dual[0]);
        gluing_change_theta_structure_dim2_to_V_compatible_theta_struct(&V_8[1], &V_8_dual[1]);

        theta_struct_dim2_t S_U, S_V;

        isogeny_dim2_compute_codomain(&S_U, U_8_dual);
        isogeny_dim2_compute_codomain(&S_V, V_8);

        theta_point_dim2_t U_8_1, U_8_2;
        theta_point_dim2_t V_8_1, V_8_2;
        isogeny_dim2_eval(&U_8_1, &U_8_dual[0], &S_U);
        isogeny_dim2_eval(&U_8_2, &U_8_dual[1], &S_U);
        isogeny_dim2_eval(&V_8_1, &V_8[0], &S_V);
        isogeny_dim2_eval(&V_8_2, &V_8[1], &S_V);

        for (int i = 0; i < 4; i++) {
            if (!fp_is_equal(&U_8_1[i], &U_8_1[i ^ 1])) {
                printf("S_U is not well defined\n");
                return 1;
            }
            if (!fp_is_equal(&U_8_2[i], &U_8_2[i ^ 2])) {
                printf("S_U is not well defined\n");
                return 1;
            }
        }
        for (int i = 0; i < 4; i++) {
            if (!fp_is_equal(&V_8_1[i], &V_8_1[i ^ 1])) {
                printf("S_V is not well defined\n");
                return 1;
            }
            if (!fp_is_equal(&V_8_2[i], &V_8_2[i ^ 2])) {
                printf("S_V is not well defined\n");
                return 1;
            }
        }
    }

    printf("Diagonal isogenies.....................................................PASSED\n");

    return 0;
}

int
test_into_dim4(int iterations)
{

    for (int g = 0; g < iterations; g++) {

        ec_curve_t E;
        ec_jac_point_t jacP, jacQ;

        test_sample_random_supersingular_curve(&E);
        ec_jac_basis_2f(&jacP, &jacQ, &E, TORSION_EVEN_POWER);

        // -------------------------------------------------------------- //
        // start true test //
        ec_jac_point_t jacP32, jacQ32;
        ec_jac_point_t jacP16, jacQ16;
        ec_jac_point_t jacP8, jacQ8;

        ec_jac_dbl_iter(&jacP32, &jacP, TORSION_EVEN_POWER - 5, &E, 1);
        ec_jac_dbl_iter(&jacQ32, &jacQ, TORSION_EVEN_POWER - 5, &E, 1);
        ec_DBL(&jacP16, &jacP32, &E);
        ec_DBL(&jacQ16, &jacQ32, &E);
        ec_DBL(&jacP8, &jacP16, &E);
        ec_DBL(&jacQ8, &jacQ16, &E);

        ec_jac_point_t T3_1_32, T3_2_32, T4_1_32, T4_2_32;
        ec_jac_point_t T3_1_16, T3_2_16, T4_1_16, T4_2_16;
        ec_jac_point_t T34_1_32, T34_2_32;
        ec_jac_point_t T34_1_16, T34_2_16;
        ec_jac_point_t tmp;

        // uint8_t sec_mod_32 = 11;
        // uint8_t sec_p_1_over_2 = (sec_mod_32 + 1)>>1; //6
        // uint8_t sec_m_1_over_2 = (sec_mod_32 - 1)>>1; //5

        ec_ADD(&tmp, &jacQ32, &jacQ16, &E);      // 3*Q_32
        ec_ADD(&T3_2_32, &tmp, &jacQ16, &E);     // 5*Q_32
        ec_ADD(&T3_1_32, &T3_2_32, &jacQ32, &E); // 6*Q_32

        ec_ADD(&tmp, &jacP32, &jacP16, &E);      // 3*P_32
        ec_ADD(&T4_2_32, &tmp, &jacP16, &E);     // 5*P_32
        ec_ADD(&T4_1_32, &T4_2_32, &jacP32, &E); // 6*P_32
        ec_jac_neg(&T4_1_32, &T4_1_32);          // -5*P_32
        ec_jac_neg(&T4_2_32, &T4_2_32);          // -6*P_32

        ec_ADD(&T34_1_32, &T3_1_32, &T4_1_32, &E); // 6(Q_32 - P_32)
        ec_ADD(&T34_2_32, &T3_2_32, &T4_2_32, &E); // 5(Q_32 - P_32)

        ec_DBL(&T3_1_16, &T3_1_32, &E);   // 5*Q_16
        ec_DBL(&T3_2_16, &T3_2_32, &E);   // 6*Q_16 = 3*Q_8
        ec_DBL(&T4_1_16, &T4_1_32, &E);   // 6*P_16 = 3*P_8
        ec_DBL(&T4_2_16, &T4_2_32, &E);   // 6*P_16 = 3*P_8
        ec_DBL(&T34_1_16, &T34_1_32, &E); // 6(Q_32 - P_32)
        ec_DBL(&T34_2_16, &T34_2_32, &E); // 5(Q_32 - P_32)

        // T1 = (6Q,  5Q)
        // T2 = (-6P, -5P)

        scholten_gluing_t scholten_gluing;
        schoten_gluing_compute_and_verify(&scholten_gluing, &jacP8, &jacQ8, &E);

        theta_point_dim2_t V_8_dual[2], U_8_dual[2];

        theta_point_dim2_t theta_T3_16[2], theta_T4_16[2], theta_T34_16[2];
        theta_point_dim2_t theta_T3_8[2], theta_T4_8[2], theta_T34_8[2];

        scholten_gluing_eval(&V_8_dual[1], &jacP16, &scholten_gluing, true); // V2
        scholten_gluing_eval(&V_8_dual[0], &jacQ16, &scholten_gluing, true); // V1

        scholten_gluing_eval(&theta_T3_16[0], &T3_1_32, &scholten_gluing, true);
        scholten_gluing_eval(&theta_T3_16[1], &T3_2_32, &scholten_gluing, true);
        scholten_gluing_eval(&theta_T4_16[0], &T4_1_32, &scholten_gluing, true);
        scholten_gluing_eval(&theta_T4_16[1], &T4_2_32, &scholten_gluing, true);

        scholten_gluing_eval(&theta_T3_8[0], &T3_1_16, &scholten_gluing, true);
        scholten_gluing_eval(&theta_T3_8[1], &T3_2_16, &scholten_gluing, true);
        scholten_gluing_eval(
            &theta_T4_8[0], &T4_1_16, &scholten_gluing, true); // rubish because barycentric coordinates
        scholten_gluing_eval(&theta_T4_8[1], &T4_2_16, &scholten_gluing, true);

        scholten_gluing_eval(&theta_T34_16[0], &T34_1_32, &scholten_gluing, true);
        scholten_gluing_eval(&theta_T34_16[1], &T34_2_32, &scholten_gluing, true);
        scholten_gluing_eval(&theta_T34_8[0], &T34_1_16, &scholten_gluing, true);
        scholten_gluing_eval(&theta_T34_8[1], &T34_2_16, &scholten_gluing, true);

        diag_gluing_t phiUV;
        gluing_diag_isogeny_compute(&phiUV, &V_8_dual[0], &V_8_dual[1]);

        //------------------------------------------------------------
        // test diagonal isogenies
        //------------------------------------------------------------
        theta_point_dim2_t phi_U_1_4, phi_U_2_4, phi_V_1_4, phi_V_2_4;

        // From V, compute U
        theta_copy_dim2(&U_8_dual[0], &V_8_dual[1]);
        fp_neg(&U_8_dual[0][3], &U_8_dual[0][3]);
        theta_copy_dim2(&U_8_dual[1], &V_8_dual[0]);
        fp_neg(&U_8_dual[1][3], &U_8_dual[1][3]);

        gluing_diag_isogeny_eval(&phi_U_1_4, &phi_V_1_4, &U_8_dual[0], &V_8_dual[0], &phiUV);
        gluing_diag_isogeny_eval(&phi_U_2_4, &phi_V_2_4, &U_8_dual[1], &V_8_dual[1], &phiUV);

        for (int i = 0; i < 4; i++) {
            if (!fp_is_equal(&phi_U_1_4[i], &phi_U_1_4[i ^ 1])) {
                printf("S_U is not well defined\n");
                return 1;
            }
            if (!fp_is_equal(&phi_U_2_4[i], &phi_U_2_4[i ^ 2])) {
                printf("S_U is not well defined\n");
                return 1;
            }
        }
        for (int i = 0; i < 4; i++) {
            if (!fp_is_equal(&phi_V_1_4[i], &phi_V_1_4[i ^ 1])) {
                printf("S_V is not well defined\n");
                return 1;
            }
            if (!fp_is_equal(&phi_V_2_4[i], &phi_V_2_4[i ^ 2])) {
                printf("S_V is not well defined\n");
                return 1;
            }
        }
        //------------------------------------------------------------
        // test dim4 theta structure
        //------------------------------------------------------------
        // evaluate points through diagonal isogenies.

        gluing_diag_isogeny_eval(&theta_T3_16[0], &theta_T3_16[1], &theta_T3_16[0], &theta_T3_16[1], &phiUV);
        gluing_diag_isogeny_eval(&theta_T4_16[0], &theta_T4_16[1], &theta_T4_16[0], &theta_T4_16[1], &phiUV);
        gluing_diag_isogeny_eval(&theta_T34_16[0], &theta_T34_16[1], &theta_T34_16[0], &theta_T34_16[1], &phiUV);

        gluing_diag_isogeny_eval(&theta_T3_8[0], &theta_T3_8[1], &theta_T3_8[0], &theta_T3_8[1], &phiUV);
        gluing_diag_isogeny_eval(&theta_T4_8[0], &theta_T4_8[1], &theta_T4_8[0], &theta_T4_8[1], &phiUV);
        gluing_diag_isogeny_eval(&theta_T34_8[0], &theta_T34_8[1], &theta_T34_8[0], &theta_T34_8[1], &phiUV);

        // theta_T4_8[0] is rubbish because of the barycentric coordinates, but we can compute it by doubling
        // theta_T4_16[0] instead.
        theta_dim2_DBL(&theta_T4_8[0], &theta_T4_16[0], &phiUV.S_U_precomp);

        theta_point_dim4_t T3_4, T4_4, T34_4;

        gluing_from_couple_dim2_to_dim4_compatible_with_isogeny(&T3_4, &theta_T3_8[0], &theta_T3_8[1]);
        gluing_from_couple_dim2_to_dim4_compatible_with_isogeny(&T4_4, &theta_T4_8[0], &theta_T4_8[1]);
        gluing_from_couple_dim2_to_dim4_compatible_with_isogeny(&T34_4, &theta_T34_8[0], &theta_T34_8[1]);

        // test that T_i are indeed 4 torsion points above the right
        for (uint8_t i = 0; i < 16; i++) {
            if (scalprod(i, 4)) {
                if (!fp_is_zero(&T3_4[i])) {
                    printf("T3_4 is not a 4 torsion point above the right base point\n");
                    return 1;
                }
            }
            if (scalprod(i, 8)) {
                if (!fp_is_zero(&T4_4[i])) {
                    printf("T4_4 is not a 4 torsion point above the right base point\n");
                    return 1;
                }
            }
            if (scalprod(i, 12)) {
                if (!fp_is_zero(&T34_4[i])) {
                    printf("T34_4 is not a 4 torsion point above the right base point\n");
                    return 1;
                }
            }
        }
    }

    printf("Into dim 4 points......................................................PASSED\n");

    return 0;
}

int
test_weil_gluing(int iterations)
{
    for (int g = 0; g < iterations; g++) {

        ec_curve_t E;
        ec_jac_point_t jacP, jacQ;

        test_sample_random_supersingular_curve(&E);
        ec_jac_basis_2f(&jacP, &jacQ, &E, TORSION_EVEN_POWER);

        // -------------------------------------------------------------- //
        // start true test //
        ec_jac_point_t jacP32, jacQ32;
        ec_jac_point_t jacP16, jacQ16;
        ec_jac_point_t jacP8, jacQ8;

        ec_jac_dbl_iter(&jacP32, &jacP, TORSION_EVEN_POWER - 5, &E, 1);
        ec_jac_dbl_iter(&jacQ32, &jacQ, TORSION_EVEN_POWER - 5, &E, 1);
        ec_DBL(&jacP16, &jacP32, &E);
        ec_DBL(&jacQ16, &jacQ32, &E);
        ec_DBL(&jacP8, &jacP16, &E);
        ec_DBL(&jacQ8, &jacQ16, &E);

        ec_jac_point_t T3_1_32, T3_2_32, T4_1_32, T4_2_32;
        ec_jac_point_t T3_1_16, T3_2_16, T4_1_16, T4_2_16;
        ec_jac_point_t T34_1_32, T34_2_32;
        ec_jac_point_t T34_1_16, T34_2_16;
        ec_jac_point_t tmp;

        // let sec_mod_32 = 11;
        //  (sec_mod_32 + 1)>>1; //6
        //  (sec_mod_32 - 1)>>1; //5

        ec_ADD(&tmp, &jacQ32, &jacQ16, &E);      // 3*Q_32
        ec_ADD(&T3_2_32, &tmp, &jacQ16, &E);     // 5*Q_32
        ec_ADD(&T3_1_32, &T3_2_32, &jacQ32, &E); // 6*Q_32

        ec_ADD(&tmp, &jacP32, &jacP16, &E);      // 3*P_32
        ec_ADD(&T4_2_32, &tmp, &jacP16, &E);     // 5*P_32
        ec_ADD(&T4_1_32, &T4_2_32, &jacP32, &E); // 6*P_32
        ec_jac_neg(&T4_1_32, &T4_1_32);          // -5*P_32
        ec_jac_neg(&T4_2_32, &T4_2_32);          // -6*P_32

        ec_ADD(&T34_1_32, &T3_1_32, &T4_1_32, &E); // 6(Q_32 - P_32)
        ec_ADD(&T34_2_32, &T3_2_32, &T4_2_32, &E); // 5(Q_32 - P_32)

        ec_DBL(&T3_1_16, &T3_1_32, &E);   // 5*Q_16
        ec_DBL(&T3_2_16, &T3_2_32, &E);   // 6*Q_16 = 3*Q_8
        ec_DBL(&T4_1_16, &T4_1_32, &E);   // 6*P_16 = 3*P_8
        ec_DBL(&T4_2_16, &T4_2_32, &E);   // 6*P_16 = 3*P_8
        ec_DBL(&T34_1_16, &T34_1_32, &E); // 6(Q_32 - P_32)
        ec_DBL(&T34_2_16, &T34_2_32, &E); // 5(Q_32 - P_32)

        // T1 = (6Q,  5Q)
        // T2 = (-6P, -5P)

        scholten_gluing_t scholten_gluing;
        schoten_gluing_compute_and_verify(&scholten_gluing, &jacP8, &jacQ8, &E);

        theta_point_dim2_t V_8_dual[2];

        theta_point_dim2_t theta_T3_16[2], theta_T4_16[2], theta_T34_16[2], theta_T3_2T4_16[2];
        // theta_point_dim2_t theta_T3_8[2], theta_T4_8[2], theta_T34_8[2];

        scholten_gluing_eval(&V_8_dual[1], &jacP16, &scholten_gluing, true); // V2
        scholten_gluing_eval(&V_8_dual[0], &jacQ16, &scholten_gluing, true); // V1

        scholten_gluing_eval(&theta_T3_16[0], &T3_1_32, &scholten_gluing, true);
        scholten_gluing_eval(&theta_T3_16[1], &T3_2_32, &scholten_gluing, true);
        scholten_gluing_eval(&theta_T4_16[0], &T4_1_32, &scholten_gluing, true);
        scholten_gluing_eval(&theta_T4_16[1], &T4_2_32, &scholten_gluing, true);
        scholten_gluing_eval(&theta_T34_16[0], &T34_1_32, &scholten_gluing, true);
        scholten_gluing_eval(&theta_T34_16[1], &T34_2_32, &scholten_gluing, true);

        diag_gluing_t phiUV;
        gluing_diag_isogeny_compute(&phiUV, &V_8_dual[0], &V_8_dual[1]);

        //------------------------------------------------------------
        // test dim4 theta structure
        //------------------------------------------------------------
        // evaluate points through diagonal isogenies.

        gluing_diag_isogeny_eval(&theta_T3_16[0], &theta_T3_16[1], &theta_T3_16[0], &theta_T3_16[1], &phiUV);
        gluing_diag_isogeny_eval(&theta_T4_16[0], &theta_T4_16[1], &theta_T4_16[0], &theta_T4_16[1], &phiUV);
        gluing_diag_isogeny_eval(&theta_T34_16[0], &theta_T34_16[1], &theta_T34_16[0], &theta_T34_16[1], &phiUV);

        // compute 4th point of dim 4 using diff_add.
        theta_point_dim2_t inv_T3_16[2];
        theta_invert_dim2(&inv_T3_16[0], &theta_T3_16[0]);
        theta_invert_dim2(&inv_T3_16[1], &theta_T3_16[1]);
        theta_dim2_diff_ADD(&theta_T3_2T4_16[0], &theta_T34_16[0], &theta_T4_16[0], &inv_T3_16[0], &phiUV.S_U_precomp);
        theta_dim2_diff_ADD(&theta_T3_2T4_16[1], &theta_T34_16[1], &theta_T4_16[1], &inv_T3_16[1], &phiUV.S_V_precomp);

        // value computed to test the gluing
        theta_point_dim2_t theta_T3_8[2];
        theta_dim2_DBL(&theta_T3_8[0], &theta_T3_16[0], &phiUV.S_U_precomp);
        theta_dim2_DBL(&theta_T3_8[1], &theta_T3_16[1], &phiUV.S_V_precomp);

        theta_point_dim2_t inv_T34_16[2], theta_T3_m_4_16[2];
        theta_invert_dim2(&inv_T34_16[0], &theta_T34_16[0]);
        theta_invert_dim2(&inv_T34_16[1], &theta_T34_16[1]);

        theta_dim2_diff_ADD(
            &theta_T3_m_4_16[0], &theta_T3_16[0], &theta_T4_16[0], &inv_T34_16[0], &phiUV.S_U_precomp); // T3 - T4
        theta_dim2_diff_ADD(
            &theta_T3_m_4_16[1], &theta_T3_16[1], &theta_T4_16[1], &inv_T34_16[1], &phiUV.S_V_precomp); // T3 - T4

        theta_point_dim4_t T3_8, T4_8, T34_8, T3_2T4_8;
        theta_point_dim4_t domain_SUxSV;
        gluing_from_couple_dim2_to_dim4_compatible_with_isogeny(
            &domain_SUxSV, &phiUV.S_U.dual_null_point, &phiUV.S_V.dual_null_point);

        gluing_from_couple_dim2_to_dim4_compatible_with_isogeny(&T3_8, &theta_T3_16[0], &theta_T3_16[1]);
        gluing_from_couple_dim2_to_dim4_compatible_with_isogeny(&T4_8, &theta_T4_16[0], &theta_T4_16[1]);
        gluing_from_couple_dim2_to_dim4_compatible_with_isogeny(&T34_8, &theta_T34_16[0], &theta_T34_16[1]);
        gluing_from_couple_dim2_to_dim4_compatible_with_isogeny(&T3_2T4_8, &theta_T3_2T4_16[0], &theta_T3_2T4_16[1]);

        theta_point_dim4_t T3_4, T3_m_T4;
        gluing_from_couple_dim2_to_dim4_compatible_with_isogeny(&T3_4, &theta_T3_8[0], &theta_T3_8[1]);
        gluing_from_couple_dim2_to_dim4_compatible_with_isogeny(&T3_m_T4, &theta_T3_m_4_16[0], &theta_T3_m_4_16[1]);

        weil_gluing_t glue_3;

        gluing_compute_weil_gluing(&glue_3, &T3_8, &T4_8, &T34_8, &T3_2T4_8);

        theta_point_dim4_t inv_dual_codomain;
        theta_weil_to_dim4_point(&inv_dual_codomain, &glue_3.weil_codom.inv_dual_null_point);

        theta_point_dim4_t codomain_sq, codomain_inv;
        theta_squared_dim4(&codomain_sq, &domain_SUxSV);
        theta_hadamard_dim4(&codomain_sq, &codomain_sq);
        theta_dot_prod_dim4(&codomain_sq, &codomain_sq, &inv_dual_codomain);

        theta_point_weil_t codom;
        gluing_special_inv(&codom, &glue_3.weil_codom.inv_dual_null_point);
        theta_weil_to_dim4_point(&codomain_inv, &codom);

        fp_t temp, temp2, temp3, temp4;
        fp_copy(&temp, &codomain_sq[0]);
        fp_copy(&temp3, &codomain_inv[0]);
        fp_inv(&temp);
        fp_inv(&temp3);

        for (uint8_t i = 0; i < 16; i++) {
            fp_mul(&temp2, &temp, &codomain_sq[i]);
            fp_mul(&temp4, &temp3, &codomain_inv[i]);
            if (!fp_is_equal(&temp2, &temp4)) {
                printf("Weil gluing does not compute the right dual theta null point \n");
                return 1;
            }
        }

        // test the the codomain is well computed.
        // compute phi(T_3_8):

        theta_point_dim4_t final_T3, final_T4;
        gluing_evaluate_weil_gluing(&final_T3, &domain_SUxSV, &T3_4, &glue_3);
        gluing_evaluate_weil_gluing(&final_T4, &T3_m_T4, &T34_8, &glue_3);

        for (uint8_t i = 0; i < 16; i++) {
            if (!fp_is_equal(&final_T3[i], &final_T3[i ^ 4])) {
                printf("final T3 point does not have the right shape \n");
                return 1;
            }
        }

        for (uint8_t i = 0; i < 16; i++) {
            if (!fp_is_equal(&final_T4[i], &final_T4[i ^ 8])) {
                printf("final T4 point does not have the right shape \n");
                return 1;
            }
        }
    }

    printf("Weil gluing............................................................PASSED\n");

    return 0;
}

int
test_mike_gluing(int iterations)
{

    for (int g = 0; g < iterations; g++) {

        ec_curve_t E;
        ec_jac_point_t jacP, jacQ;

        test_sample_random_supersingular_curve(&E);
        ec_jac_basis_2f(&jacP, &jacQ, &E, TORSION_EVEN_POWER);

        // -------------------------------------------------------------- //
        // start true test //
        ec_jac_point_t jacP64, jacQ64;
        ec_jac_point_t jacP32, jacQ32;
        ec_jac_point_t jacP16, jacQ16;
        ec_jac_point_t jacP8, jacQ8;

        ec_jac_dbl_iter(&jacP64, &jacP, TORSION_EVEN_POWER - 6, &E, 1);
        ec_jac_dbl_iter(&jacQ64, &jacQ, TORSION_EVEN_POWER - 6, &E, 1);

        ec_DBL(&jacP32, &jacP64, &E);
        ec_DBL(&jacQ32, &jacQ64, &E);
        ec_DBL(&jacP16, &jacP32, &E);
        ec_DBL(&jacQ16, &jacQ32, &E);
        ec_DBL(&jacP8, &jacP16, &E);
        ec_DBL(&jacQ8, &jacQ16, &E);

        mike_gluing_basis_t mike_basis;

        ec_copy_jac_point(&mike_basis.P_8, &jacP8);
        ec_copy_jac_point(&mike_basis.Q_8, &jacQ8);
        ec_copy_jac_point(&mike_basis.P_16, &jacP16);
        ec_copy_jac_point(&mike_basis.Q_16, &jacQ16);

        // Compute gluing basis using x_sec = 11 md 32
        ec_jac_point_t tmp;
        ec_ADD(&tmp, &jacQ32, &jacQ16, &E);                            // 3*Q_32
        ec_ADD(&mike_basis.T3_2_32, &tmp, &jacQ16, &E);                // 5*Q_32
        ec_ADD(&mike_basis.T3_1_32, &mike_basis.T3_2_32, &jacQ32, &E); // 6*Q_32

        ec_ADD(&tmp, &jacP32, &jacP16, &E);                            // 3*P_32
        ec_ADD(&mike_basis.T4_2_32, &tmp, &jacP16, &E);                // 5*P_32
        ec_ADD(&mike_basis.T4_1_32, &mike_basis.T4_2_32, &jacP32, &E); // 6*P_32

        ec_jac_neg(&mike_basis.T4_1_32, &mike_basis.T4_1_32); // -5*P_32
        ec_jac_neg(&mike_basis.T4_2_32, &mike_basis.T4_2_32); // -6*P_32

        mike_gluing_t mike_glue;

        gluing_mike_compute(&mike_glue, &E, &mike_basis);

        theta_point_dim4_t T3_4, T4_4;

        gluing_mike_eval(&T3_4, &mike_basis.T3_1_32, &mike_basis.T3_2_32, &mike_glue);
        gluing_mike_eval(&T4_4, &mike_basis.T4_1_32, &mike_basis.T4_2_32, &mike_glue);

        theta_hadamard_dim4(&T3_4, &T3_4);
        theta_hadamard_dim4(&T4_4, &T4_4);

        for (uint8_t i = 0; i < 16; i++) {

            if (!fp_is_equal(&T3_4[i], &T3_4[i ^ 4])) {
                printf("T3_4 by isog eval does not have the right shape \n");
                return 1;
            }
        }
        for (uint8_t i = 0; i < 16; i++) {

            if (!fp_is_equal(&T4_4[i], &T4_4[i ^ 8])) {
                printf("T4_4 by isog eval does not have the right shape \n");
                return 1;
            }
        }

        theta_hadamard_dim4(&T3_4, &T3_4);
        theta_hadamard_dim4(&T4_4, &T4_4);

        //------------------------------------------------------------
        // test doubling on dim4
        //------------------------------------------------------------

        theta_struct_weil_precomp_t Weil_fold;

        theta_copy_weil_struct_precomp(&Weil_fold, &mike_glue.weil_gluing.weil_precomp);

        fp_t tmp_1, tmp_2, tmp_3, tmp_4;

        // test values above using
        ec_jac_point_t T3_1_64, T3_2_64;
        ec_jac_point_t T4_1_64, T4_2_64;
        ec_ADD(&tmp, &jacQ64, &jacQ32, &E);      // 3*Q_64
        ec_ADD(&T3_2_64, &tmp, &jacQ32, &E);     // 5*Q_64
        ec_ADD(&T3_1_64, &T3_2_64, &jacQ64, &E); // 6*Q_64

        ec_ADD(&tmp, &jacP64, &jacP32, &E);      // 3*P_64
        ec_ADD(&T4_2_64, &tmp, &jacP32, &E);     // 5*P_64
        ec_ADD(&T4_1_64, &T4_2_64, &jacP64, &E); // 6*P_64

        theta_point_dim4_t T3_8, T4_8;

        gluing_mike_eval(&T3_8, &T3_1_64, &T3_2_64, &mike_glue);
        gluing_mike_eval(&T4_8, &T4_1_64, &T4_2_64, &mike_glue);

        theta_point_dim4_t T3_2, T4_2, T3_0, T4_0;
        theta_DBL_weil(&T3_2, &T3_8, &Weil_fold); // 4 torsion
        theta_DBL_weil(&T4_2, &T4_8, &Weil_fold); // 4 torsion

        theta_hadamard_dim4(&T3_2, &T3_2);
        theta_hadamard_dim4(&T4_2, &T4_2);

        for (uint8_t i = 0; i < 16; i++) {

            if (!fp_is_equal(&T3_2[i], &T3_2[i ^ 4])) {
                printf("T3_4 by doubling does not have the right shape \n");
                return 1;
            }
        }
        for (uint8_t i = 0; i < 16; i++) {

            if (!fp_is_equal(&T4_2[i], &T4_2[i ^ 8])) {
                printf("T4_4 by doubling does not have the right shape \n");
                return 1;
            }
        }

        theta_hadamard_dim4(&T3_2, &T3_2);
        theta_hadamard_dim4(&T4_2, &T4_2);

        theta_DBL_weil(&T3_2, &T3_2, &Weil_fold); // 2 torsion
        theta_DBL_weil(&T3_0, &T3_2, &Weil_fold); // 0 torsion
        theta_DBL_weil(&T4_2, &T4_2, &Weil_fold); // 2 torsion
        theta_DBL_weil(&T4_0, &T4_2, &Weil_fold); // 0 torsion

        fp_copy(&tmp_1, &T4_2[0]);
        fp_inv(&tmp_1);
        fp_copy(&tmp_3, &T4_0[0]);
        fp_inv(&tmp_3);
        fp_t tmp_5;

        for (uint8_t i = 0; i < 16; i++) {

            fp_mul(&tmp_2, &tmp_1, &T4_2[i]);
            fp_mul(&tmp_4, &tmp_3, &T4_0[i]);

            if (scalprod(i, 8) == 0) {
                fp_sub(&tmp_5, &tmp_2, &tmp_4);
            } else {
                fp_add(&tmp_5, &tmp_2, &tmp_4);
            }

            if (!fp_is_zero(&tmp_5)) {
                printf("T4 points are not above a 2 torsion point\n");
                return 1;
            }
        }

        fp_copy(&tmp_1, &T3_2[0]);
        fp_inv(&tmp_1);
        fp_copy(&tmp_3, &T3_0[0]);
        fp_inv(&tmp_3);
        for (uint8_t i = 0; i < 16; i++) {

            fp_mul(&tmp_2, &tmp_1, &T3_2[i]);
            fp_mul(&tmp_4, &tmp_3, &T3_0[i]);

            if (scalprod(i, 4) == 0) {
                fp_sub(&tmp_5, &tmp_2, &tmp_4);
            } else {
                fp_add(&tmp_5, &tmp_2, &tmp_4);
            }

            if (!fp_is_zero(&tmp_5)) {
                printf("T3 points are not above a 2 torsion point\n");
                return 1;
            }
        }
    }

    printf("Mike gluing............................................................PASSED\n");

    return 0;
}

int
test_gluing_basis_gen(int iterations)
{

    for (int g = 0; g < iterations; g++) {

        ec_curve_t E;
        ec_jac_point_t jacP, jacQ;

        test_sample_random_supersingular_curve(&E);
        ec_jac_basis_2f(&jacP, &jacQ, &E, TORSION_EVEN_POWER);

        // -------------------------------------------------------------- //
        // start true test //
        ec_jac_point_t jacP64, jacQ64;
        ec_jac_point_t jacPQ32[2];

        ec_jac_dbl_iter(&jacP64, &jacP, TORSION_EVEN_POWER - 6, &E, 1);
        ec_jac_dbl_iter(&jacQ64, &jacQ, TORSION_EVEN_POWER - 6, &E, 1);

        ec_DBL(&jacPQ32[0], &jacP64, &E);
        ec_DBL(&jacPQ32[1], &jacQ64, &E);

        uint8_t x_32 = ((3 + 8 * g) & 63) >> 1;

        mike_gluing_basis_t mike_basis;

        gluing_basis_compute(&mike_basis, &E, jacPQ32, x_32);

        mike_gluing_t mike_glue;

        gluing_mike_compute(&mike_glue, &E, &mike_basis);

        theta_point_dim4_t T3_4, T4_4;

        gluing_mike_eval(&T3_4, &mike_basis.T3_1_32, &mike_basis.T3_2_32, &mike_glue);
        gluing_mike_eval(&T4_4, &mike_basis.T4_1_32, &mike_basis.T4_2_32, &mike_glue);

        theta_hadamard_dim4(&T3_4, &T3_4);
        theta_hadamard_dim4(&T4_4, &T4_4);

        for (uint8_t i = 0; i < 16; i++) {

            if (!fp_is_equal(&T3_4[i], &T3_4[i ^ 4])) {
                printf("T3_4 by isog eval does not have the right shape \n");
                return 1;
            }
        }
        for (uint8_t i = 0; i < 16; i++) {

            if (!fp_is_equal(&T4_4[i], &T4_4[i ^ 8])) {
                printf("T4_4 by isog eval does not have the right shape \n");
                return 1;
            }
        }

        theta_hadamard_dim4(&T3_4, &T3_4);
        theta_hadamard_dim4(&T4_4, &T4_4);
    }

    printf("Gluing basis...........................................................PASSED\n");

    return 0;
}

int
test_dim4_chain(int iterations)
{

    for (int g = 0; g < iterations; g++) {

        ec_curve_t E;
        ec_jac_point_t jacP, jacQ;
        ec_jac_point_t jacP_len, jacQ_len;

        test_sample_random_supersingular_curve(&E);
        ec_jac_basis_2f(&jacP, &jacQ, &E, TORSION_EVEN_POWER);

        int LENGTH = 192;

        digit_t sk[3];
        randombytes((unsigned char *)sk, 3 * sizeof(digit_t));
        sk[0] = sk[0] - (sk[0] & 7) + 3;

        ec_jac_dbl_iter(&jacP_len, &jacP, TORSION_EVEN_POWER - LENGTH, &E, 1);
        ec_jac_dbl_iter(&jacQ_len, &jacQ, TORSION_EVEN_POWER - LENGTH, &E, 1);

        mike_chain_basis_t mike_basis;
        mike_gluing_basis_t gluing_basis;

        // compute x >> 1
        uint8_t sk_bshift1_mod32 = ((sk[0] >> 1) & 31);

        // ec_jac_point_t *P, const ec_jac_point_t *Q, const digit_t *x,const int nwords , const int pow , ec_curve_t *E
        mike_compute_basis(&mike_basis, &jacP_len, &jacQ_len, sk, LENGTH, &E);

        ec_jac_point_t jacPQ32[2];
        ec_jac_dbl_iter(&jacPQ32[0], &jacP_len, LENGTH - 5, &E, 1);
        ec_jac_dbl_iter(&jacPQ32[1], &jacQ_len, LENGTH - 5, &E, 1);

        mike_compute_basis(&mike_basis, &jacP_len, &jacQ_len, sk, LENGTH, &E);

        gluing_basis_compute(&gluing_basis, &E, jacPQ32, sk_bshift1_mod32 ); 


        ec_jac_point_t basis_DBL[4]; 

        ec_jac_dbl_iter(&basis_DBL[0], &mike_basis.T3_1, LENGTH - 5 , &E, 1); 
        ec_jac_dbl_iter(&basis_DBL[1], &mike_basis.T3_2, LENGTH - 5 , &E, 1); 
        ec_jac_dbl_iter(&basis_DBL[2], &mike_basis.T4_1, LENGTH - 5 , &E, 1); 
        ec_jac_dbl_iter(&basis_DBL[3], &mike_basis.T4_2, LENGTH - 5 , &E, 1); 

        if(!ec_jac_is_equal(&basis_DBL[0], &gluing_basis.T3_1_32))
        {
            printf("mike basis and gluing basis are not consistent: T3_1");
            return 1;
        }
        if(!ec_jac_is_equal(&basis_DBL[1], &gluing_basis.T3_2_32))
        {
            printf("mike basis and gluing basis are not consistent: T3_2");
            return 1;
        }
        if(!ec_jac_is_equal(&basis_DBL[2], &gluing_basis.T4_1_32))
        {
            printf("mike basis and gluing basis are not consistent: T4_1");
            return 1;
        }
        if(!ec_jac_is_equal(&basis_DBL[3], &gluing_basis.T4_2_32))
        {
            printf("mike basis and gluing basis are not consistent: T4_2");
            return 1;
        }

        // Test dim 4
        mike_abs_invariants_t invariant;
        mike_isogeny_chain_dim4(&invariant, &mike_basis, &gluing_basis, &E, LENGTH - 2);

        // if it passes, then the chain should be correct
    }

    printf("Dim 4 chain............................................................PASSED\n");
    return 0;
}

int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int iterations = 10 * MIKE_TEST_REPS;
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

    print_seed(seed);

#if defined(TARGET_BIG_ENDIAN)
    for (int i = 0; i < 12; i++) {
        seed[i] = BSWAP32(seed[i]);
    }
#endif

    randombytes_init((unsigned char *)seed, NULL, 256);

    res = test_theta_structure_scholten(iterations);
    res |= test_scholten_formulae(iterations);
    res |= test_scholten_gluing(iterations);
    res |= test_diagonal_isogenies(iterations);
    res |= test_into_dim4(iterations);
    res |= test_weil_gluing(iterations);
    res |= test_mike_gluing(iterations);
    res |= test_gluing_basis_gen(iterations);
    res |= test_dim4_chain(10);
    // res |= 1; // Replace by unit tests (output 0 on success) here

    if (res) {
        printf("Tests failed!\n");
    } else {
        printf("All Weil module tests passed.\n");
    }

    return res;
}
