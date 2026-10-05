#include <weil.h>
#include "assert.h"
#include "stdint.h"

const uint8_t gluing_change_struct_dim4[16] = { 0, 5, 10, 15, 1, 4, 11, 14, 8, 13, 2, 7, 9, 12, 3, 6 };

void
scholten_variable_compute(weil_restriction_theta_basis_t *coeff_basis, const ec_xz_point_t *P)
{
    fp_t x_n, z_n;
    fp_t a, b;

    fp_sqr(&a, &P->x.re);
    fp_sqr(&b, &P->x.im);
    fp_add(&x_n, &a, &b); // n(x)
    fp_sqr(&a, &P->z.re);
    fp_sqr(&b, &P->z.im);
    fp_add(&z_n, &a, &b); // n(z)

    fp_add(&(*coeff_basis)[0], &x_n, &z_n); // coeff_basis_0 = n(x) + n(z)
    fp_sub(&(*coeff_basis)[2], &z_n, &x_n); // coeff_basis_2 = n(z) - n(x)

    fp_mul(&a, &P->x.re, &P->z.re);
    fp_mul(&b, &P->x.im, &P->z.im);
    fp_add(&(*coeff_basis)[1], &a, &b);
    fp_add(&(*coeff_basis)[1], &(*coeff_basis)[1], &(*coeff_basis)[1]); // coeff_basis_1 = 2tr(x_P*conj(z_P))

    fp_mul(&a, &P->x.im, &P->z.re);
    fp_mul(&b, &P->x.re, &P->z.im);
    fp_sub(&(*coeff_basis)[3], &a, &b);
    fp_add(&(*coeff_basis)[3], &(*coeff_basis)[3], &(*coeff_basis)[3]); // coeff_basis_3 = 2x_P*comj(z_P)/i
}

void
scholten_jacobian_to_theta(theta_point_dim2_t *out, const weil_restriction_theta_basis_t *m, const ec_xz_point_t *P)
{
    fp_t n_x, n_z;
    fp_t n_add, n_sub;
    fp_t t_xz, ti_xz;
    fp_t a, b;

    fp_sqr(&a, &P->x.re);
    fp_sqr(&b, &P->x.im);
    fp_add(&n_x, &a, &b); // n(x)

    fp_sqr(&a, &P->z.re);
    fp_sqr(&b, &P->z.im);
    fp_add(&n_z, &a, &b); // n(z)

    fp_add(&n_add, &n_x, &n_z); // n(x) + n(z)
    fp_sub(&n_sub, &n_x, &n_z); // n(x) - n(z)

    fp_mul(&a, &P->x.re, &P->z.re);
    fp_mul(&b, &P->x.im, &P->z.im);
    fp_add(&t_xz, &a, &b);
    fp_add(&t_xz, &t_xz, &t_xz); // 2(x*conj(z) + conj(x}*z)

    fp_mul(&a, &P->x.im, &P->z.re);
    fp_mul(&b, &P->x.re, &P->z.im);
    fp_sub(&ti_xz, &a, &b);
    fp_add(&ti_xz, &ti_xz, &ti_xz); // 2(x*conj(z) - conj(x}*z)/i

    fp_mul(&(*out)[0], &(*m)[0], &n_add);
    fp_mul(&a, &(*m)[1], &t_xz);
    fp_sub(&(*out)[0], &(*out)[0], &a);

    fp_mul(&(*out)[1], &(*m)[0], &t_xz);
    fp_mul(&a, &(*m)[1], &n_add);
    fp_sub(&(*out)[1], &(*out)[1], &a);

    fp_mul(&(*out)[2], &(*m)[2], &n_sub);
    fp_mul(&a, &(*m)[3], &ti_xz);
    fp_sub(&(*out)[2], &(*out)[2], &a); // sub from hiden factor by i

    fp_mul(&(*out)[3], &(*m)[2], &ti_xz);
    fp_mul(&a, &(*m)[3], &n_sub);
    fp_add(&(*out)[3], &(*out)[3], &a); // note, it is off by a factor i
}

void
scholten_precompute(scholten_precomp_t *scholten, const weil_restriction_theta_basis_t *coeff_basis)
{
    fp_t m0, m1, m2;

    fp_mul(&scholten->D, &(*coeff_basis)[0], &(*coeff_basis)[1]);
    fp_add(&scholten->D, &scholten->D, &scholten->D);
    fp_add(&scholten->D, &scholten->D, &scholten->D);
    fp_mul(&scholten->E, &(*coeff_basis)[2], &(*coeff_basis)[3]);
    fp_add(&scholten->E, &scholten->E, &scholten->E);
    fp_add(&scholten->E, &scholten->E, &scholten->E);

    fp_sqr(&m0, &(*coeff_basis)[0]);
    fp_sqr(&m1, &(*coeff_basis)[1]);
    fp_sqr(&m2, &(*coeff_basis)[2]);

    fp_add(&scholten->A, &m1, &m2);
    fp_sub(&scholten->B, &m0, &m1);
    fp_sub(&scholten->C, &m0, &m2);
}

void
scholten_eval(theta_point_dim2_t *out, const ec_bary_coords_t *uvw, const scholten_precomp_t *scholten, bool twist)
{

    fp2_t u2, v2, w2;
    fp2_t uw, u2_v2_p_w2, u2_v2_m_w2;
    fp_t u_n, v_n, w_n;
    fp_t uv_p_w_n, uv_m_w_n, u_nw_n, v_nw_n;

    fp_t U, V;
    fp_t a, b, c;
    fp2_t tmp;

    // compute u^2 and n(u)
    fp_sqr(&a, &uvw->u.re);
    fp_sqr(&b, &uvw->u.im);
    fp_add(&u_n, &a, &b);     // u_r^2 + u_i^2 = n(u)
    fp_add(&u_n, &u_n, &u_n); // 2n(u)

    fp_sub(&u2.re, &a, &b); // u_r^2 - u_i^2
    fp_mul(&a, &uvw->u.re, &uvw->u.im);
    fp_add(&u2.im, &a, &a); // 2u_ru_i

    // compute v^2 and n(v)
    fp_sqr(&a, &uvw->v.re);
    fp_sqr(&b, &uvw->v.im);
    fp_add(&v_n, &a, &b);     // v_r^2 + v_i^2 = n(v)
    fp_add(&v_n, &v_n, &v_n); // 2n(v)

    fp_sub(&v2.re, &a, &b); // v_r^2 - v_i^2
    fp_mul(&a, &uvw->v.re, &uvw->v.im);
    fp_add(&v2.im, &a, &a); // 2w_rw_i

    // compute w^2 and n(w)
    fp_sqr(&a, &uvw->w.re);
    fp_sqr(&b, &uvw->w.im);
    fp_add(&w_n, &a, &b);     // w_r^2 + w_i^2 = n(w)
    fp_add(&w_n, &w_n, &w_n); // 2n(w)

    fp_sub(&w2.re, &a, &b); // w_r^2 - w_i^2
    fp_mul(&a, &uvw->w.re, &uvw->w.im);
    fp_add(&w2.im, &a, &a); // 2w_rw_i

    // compute u*w
    fp2_mul(&uw, &uvw->u, &uvw->w);

    // compute product of norms
    fp_mul(&u_nw_n, &u_n, &w_n); // 4n(u)n(w)
    fp_mul(&v_nw_n, &v_n, &w_n); // 4n(v)n(w)

    // compute u^2 - v^2 ± w^2
    fp2_sub(&tmp, &u2, &v2);
    fp2_add(&u2_v2_p_w2, &tmp, &w2); // u^2 - v^2 + w^2
    fp2_sub(&u2_v2_m_w2, &tmp, &w2); // u^2 - v^2 - w^2

    // compute n(u^2 - v^2 ± w^2)
    fp_sqr(&a, &u2_v2_p_w2.re);
    fp_sqr(&b, &u2_v2_p_w2.im);
    fp_add(&uv_p_w_n, &a, &b); // n(u^2 - v^2 + w^2)

    fp_sqr(&a, &u2_v2_m_w2.re);
    fp_sqr(&b, &u2_v2_m_w2.im);
    fp_add(&uv_m_w_n, &a, &b); // n(u^2 - v^2 - w^2)

    // compute U + iV = (u^2 - v^2 + w^2)* bar(u * w)
    fp_mul(&a, &u2_v2_p_w2.re, &uw.re);
    fp_mul(&b, &u2_v2_p_w2.im, &uw.im);
    fp_add(&U, &a, &b);

    fp_mul(&a, &u2_v2_p_w2.im, &uw.re);
    fp_mul(&b, &u2_v2_p_w2.re, &uw.im);
    fp_sub(&V, &a, &b);

    fp_mul(&a, &U, &scholten->D);
    fp_mul(&b, &V, &scholten->E);

    // First component
    fp_mul(&(*out)[0], &scholten->A, &uv_p_w_n);
    fp_mul(&c, &scholten->C, &u_nw_n);
    fp_add(&(*out)[0], &(*out)[0], &c);
    fp_sub(&(*out)[0], &(*out)[0], &a);
    fp_sub(&(*out)[0], &(*out)[0], &b);

    // Third component
    fp_mul(&(*out)[2], &scholten->C, &uv_p_w_n);
    fp_mul(&c, &scholten->A, &u_nw_n);
    fp_add(&(*out)[2], &(*out)[2], &c);
    fp_sub(&(*out)[2], &(*out)[2], &a);
    fp_add(&(*out)[2], &(*out)[2], &b);

    // Second component
    fp_mul(&(*out)[1], &scholten->B, &uv_m_w_n);

    // Last component
    fp_mul(&(*out)[3], &scholten->B, &v_nw_n);

    // In mike, twist is always true, but used for some debug
    if (twist) {
        fp_neg(&(*out)[3], &(*out)[3]);
    }
}

uint32_t
schoten_gluing_compute_and_verify(scholten_gluing_t *phi,
                                  const ec_jac_point_t *T1_8,
                                  const ec_jac_point_t *T2_8,
                                  const ec_curve_t *E)
{

    ec_copy_curve(&phi->E, E);
    ec_copy_jac_point(&phi->aux_point, T1_8);

    ec_xz_point_t xz_T1, xz_T2;
    ec_jac_point_t T1, T2;
    uint32_t res = -1;
    fp2_t Q4z2, mQ4z2, testA;

    // Test the kernel (this ensures, by Hass-bound, that E is supersingular)

    // Check A is defined over Fp2
    res &= ~fp_is_zero(&E->A.im);

    ec_DBL(&T1, T1_8, E); // 4 torsion point.
    ec_jac_to_xz(&xz_T1, &T1);
    ec_DBL(&T1, &T1, E); // 2 torsion point

    // T1 = (*:0:*)
    res &= fp2_is_zero(&T1.y);

    ec_DBL(&T2, T2_8, E);

    // T2 == (± 1:*:1) <==> T2.x == T2.z^2 or T2.x == -T2.z^2
    fp2_sqr(&Q4z2, &T2.z);
    fp2_neg(&mQ4z2, &Q4z2);
    res &= (fp2_is_equal(&T2.x, &Q4z2) | fp2_is_equal(&T2.x, &mQ4z2));

    // Check that the Montgomery coeff is the maximal one.
    fp2_copy(&testA, &E->A);
    ec_normalize_montgomery(&testA, &xz_T1, 2);
    res &= fp2_is_equal(&testA, &E->A);

    // Compute values change of theta structure.
    weil_restriction_theta_basis_t weil_basis;
    scholten_variable_compute(&weil_basis, &xz_T1);
    scholten_precompute(&phi->precomp, &weil_basis);

    // compute codomain
    theta_point_dim2_t HSK1, HSK2;

    ec_jac_to_xz(&xz_T1, T1_8);
    ec_jac_to_xz(&xz_T2, T2_8);

    scholten_jacobian_to_theta(&HSK1, &weil_basis, &xz_T1);
    scholten_jacobian_to_theta(&HSK2, &weil_basis, &xz_T2);
    theta_squared_dim2(&HSK1, &HSK1);
    theta_squared_dim2(&HSK2, &HSK2);
    fp_neg(&HSK1[3], &HSK1[3]);
    fp_neg(&HSK2[3], &HSK2[3]);
    theta_hadamard_dim2(&HSK1, &HSK1);
    theta_hadamard_dim2(&HSK2, &HSK2);

    // HSK1 = (Ax: Bx: Cy: 0)
    // HSK2 = (Az: Bw: Cz: 0)
    fp_mul(&phi->inv_image_aux_point[2], &HSK1[0], &HSK2[2]); // AxCz
    fp_copy(&phi->inv_image_aux_point[3], &phi->inv_image_aux_point[2]);

    fp_mul(&phi->inv_image_aux_point[0], &HSK1[2], &HSK2[0]); // CyAz
    fp_copy(&phi->inv_image_aux_point[1], &phi->inv_image_aux_point[0]);

    return res;
}

void
scholten_gluing_eval(theta_point_dim2_t *out, const ec_jac_point_t *P, const scholten_gluing_t *phi, bool twist)
{

    ec_bary_coords_t uvw;
    ec_jac_to_bary_coords(&uvw, P, &phi->aux_point, &phi->E);

    // If P = (0:1:0) --> (u:v:w) = (1:0:0)
    if (fp2_is_zero(&uvw.w)) {
        fp2_set_one(&uvw.u);
        fp2_set_zero(&uvw.v);
    }

    theta_point_dim2_t prod;
    scholten_eval(&prod, &uvw, &phi->precomp, twist);

    for (uint8_t i = 0; i < 4; i++) {
        fp_mul(&(*out)[i], &prod[i], &phi->inv_image_aux_point[i]);
    }
}

// Compute H * D * H
void
gluing_change_theta_structure_dim2_to_V_compatible_theta_struct(theta_point_dim2_t *out, const theta_point_dim2_t *in)
{
    fp_t in1[2], in2[2];
    fp_add(&in1[0], &(*in)[0], &(*in)[1]); // a + b
    fp_sub(&in1[1], &(*in)[0], &(*in)[1]); // a - b

    fp_add(&in2[0], &(*in)[2], &(*in)[3]); // c + d
    fp_sub(&in2[1], &(*in)[2], &(*in)[3]); // c - d

    fp_add(&(*out)[0], &in1[0], &in2[1]); //  a + b + c - d
    fp_sub(&(*out)[1], &in1[0], &in2[1]); //  a + b - c + d
    fp_add(&(*out)[2], &in2[0], &in1[1]); //  a - b + c + d
    fp_sub(&(*out)[3], &in2[0], &in1[1]); // -a + b + c + d
}

// Compute D * H
void
gluing_change_theta_structure_dim2_to_U_compatible_theta_struct(theta_point_dim2_t *out, const theta_point_dim2_t *in)
{
    fp_t in1[2], in2[2];
    fp_add(&in1[0], &(*in)[0], &(*in)[1]); // a + b
    fp_sub(&in1[1], &(*in)[0], &(*in)[1]); // a - b

    fp_add(&in2[0], &(*in)[2], &(*in)[3]); // c + d
    fp_sub(&in2[1], &(*in)[2], &(*in)[3]); // c - d

    fp_add(&(*out)[0], &in1[0], &in2[0]); //  a + b + c + d
    fp_add(&(*out)[1], &in1[1], &in2[1]); //  a - b + c - d
    fp_sub(&(*out)[2], &in1[0], &in2[0]); //  a + b - c - d
    fp_sub(&(*out)[3], &in2[1], &in1[1]); // -a + b + c - d
}

void
gluing_diag_isogeny_compute(diag_gluing_t *phiUV, const theta_point_dim2_t *V1_8, const theta_point_dim2_t *V2_8)
{
    theta_point_dim2_t ker_V[2];
    theta_point_dim2_t ker_U[2];

    theta_copy_dim2(&ker_U[0], V2_8);
    fp_neg(&ker_U[0][3], &ker_U[0][3]); // Phi_1(P, pi(P))
    theta_copy_dim2(&ker_U[1], V1_8);
    fp_neg(&ker_U[1][3], &ker_U[1][3]); // Phi_1(Q, pi(Q))

    gluing_change_theta_structure_dim2_to_U_compatible_theta_struct(&ker_U[0], &ker_U[0]);
    gluing_change_theta_structure_dim2_to_U_compatible_theta_struct(&ker_U[1], &ker_U[1]);
    gluing_change_theta_structure_dim2_to_V_compatible_theta_struct(&ker_V[0], V1_8);
    gluing_change_theta_structure_dim2_to_V_compatible_theta_struct(&ker_V[1], V2_8);

    isogeny_dim2_compute_codomain(&phiUV->S_U, ker_U);
    isogeny_dim2_compute_codomain(&phiUV->S_V, ker_V);

    theta_struct_dim2_arith_precomp(&phiUV->S_U_precomp, &phiUV->S_U);
    theta_struct_dim2_arith_precomp(&phiUV->S_V_precomp, &phiUV->S_V);
}

// Apply phi_U to 1st coeff and phi_V to 2nd coeff. Output is in dual theta mode for both
void
gluing_diag_isogeny_eval(theta_point_dim2_t *out1,
                         theta_point_dim2_t *out2,
                         const theta_point_dim2_t *in1,
                         const theta_point_dim2_t *in2,
                         const diag_gluing_t *phiUV)
{

    theta_point_dim2_t inter_1, inter_2;
    // Apply change of coordinates
    // H * H * D * H = D * H
    gluing_change_theta_structure_dim2_to_U_compatible_theta_struct(&inter_1, in1);
    // H * D * H
    gluing_change_theta_structure_dim2_to_V_compatible_theta_struct(&inter_2, in2);

    // Compute the respective dual theta point
    isogeny_dim2_eval(out1, &inter_1, &phiUV->S_U);
    isogeny_dim2_eval(out2, &inter_2, &phiUV->S_V);
}

static inline void
kronecker_product(theta_point_dim4_t *out, theta_point_dim2_t *in1, theta_point_dim2_t *in2)
{
    uint8_t i, j;
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            fp_mul(&(*out)[4 * i + j], &(*in1)[i], &(*in2)[j]);
        }
    }
}

static inline void
gluing_change_theta_struct_dim4_to_weil_gluing_compatible_theta_struct(theta_point_dim4_t *out,
                                                                       const theta_point_dim4_t *in)
{
    fp_t in1[2], in2[2];
    for (uint8_t i = 0; i < 4; i++) {
        fp_add(
            &in1[0], &(*in)[gluing_change_struct_dim4[4 * i]], &(*in)[gluing_change_struct_dim4[4 * i + 1]]); // a + b
        fp_sub(
            &in1[1], &(*in)[gluing_change_struct_dim4[4 * i]], &(*in)[gluing_change_struct_dim4[4 * i + 1]]); // a - b

        fp_add(&in2[0],
               &(*in)[gluing_change_struct_dim4[4 * i + 2]],
               &(*in)[gluing_change_struct_dim4[4 * i + 3]]); // c + d
        fp_sub(&in2[1],
               &(*in)[gluing_change_struct_dim4[4 * i + 2]],
               &(*in)[gluing_change_struct_dim4[4 * i + 3]]); // c - d

        fp_add(&(*out)[4 * i], &in1[0], &in2[0]);     // a + b + c + d
        fp_add(&(*out)[4 * i + 1], &in1[1], &in2[1]); // a - b + c - d
        fp_sub(&(*out)[4 * i + 2], &in1[0], &in2[0]); // a + b - c - d
        fp_sub(&(*out)[4 * i + 3], &in1[1], &in2[1]); // a - b - c + d
    }
}

void
gluing_from_couple_dim2_to_dim4_compatible_with_isogeny(theta_point_dim4_t *out,
                                                        theta_point_dim2_t *in1,
                                                        theta_point_dim2_t *in2)
{
    theta_point_dim2_t dual_in2;
    theta_point_dim4_t tmp;

    theta_hadamard_dim2(&dual_in2, in2);
    kronecker_product(&tmp, in1, &dual_in2);
    gluing_change_theta_struct_dim4_to_weil_gluing_compatible_theta_struct(out, &tmp);
}

// Weil gluing

// special function used to inverse the theta_null_point of the gluing isogeny (because of zeros)
void
gluing_special_inv(theta_point_weil_t *inv_vec, const theta_point_weil_t *vec)
{
    fp_t non_zero_points[7];

    fp_copy(&non_zero_points[0], &(*vec)[0]);
    fp_copy(&non_zero_points[1], &(*vec)[1]);
    fp_copy(&non_zero_points[2], &(*vec)[2]);
    fp_copy(&non_zero_points[3], &(*vec)[3]);
    fp_copy(&non_zero_points[4], &(*vec)[5]);
    fp_copy(&non_zero_points[5], &(*vec)[7]);
    fp_copy(&non_zero_points[6], &(*vec)[9]);

    fp_proj_batched_inv(non_zero_points, 7);

    fp_copy(&(*inv_vec)[0], &non_zero_points[0]);
    fp_copy(&(*inv_vec)[1], &non_zero_points[1]);
    fp_copy(&(*inv_vec)[2], &non_zero_points[2]);
    fp_copy(&(*inv_vec)[3], &non_zero_points[3]);
    fp_copy(&(*inv_vec)[5], &non_zero_points[4]);
    fp_copy(&(*inv_vec)[7], &non_zero_points[5]);
    fp_copy(&(*inv_vec)[9], &non_zero_points[6]);

    fp_set_zero(&(*inv_vec)[4]);
    fp_set_zero(&(*inv_vec)[6]);
    fp_set_zero(&(*inv_vec)[8]);
}

// Solve the HIIP on Weil gluing isogeny
static void
solve_gluing_HIIP(weil_gluing_t *gluing, const theta_point_dim4_t *HSK)
{

    // intermediary values
    fp_t leg_4[3], leg_8[3], leg_12[3];

    // intermediary values
    fp_t aux_point_compress[8];

    fp_t tmp;

    // leg_4
    fp_mul(&leg_4[0], &HSK[0][4], &HSK[0][6]);
    fp_mul(&leg_4[1], &HSK[0][0], &HSK[0][6]);
    fp_mul(&leg_4[2], &HSK[0][0], &HSK[0][2]);

    // leg_8
    fp_mul(&leg_8[0], &HSK[1][8], &HSK[1][9]);
    fp_mul(&leg_8[1], &HSK[1][0], &HSK[1][9]);
    fp_mul(&leg_8[2], &HSK[1][0], &HSK[1][1]);

    // leg_12
    fp_mul(&leg_12[0], &HSK[2][12], &HSK[2][15]);
    fp_mul(&leg_12[1], &HSK[2][0], &HSK[2][15]);
    fp_mul(&leg_12[2], &HSK[2][0], &HSK[2][3]);

    // compute inverse_theta_null_point:
    fp_mul(&tmp, &leg_8[0], &leg_12[0]);
    fp_mul(&gluing->weil_codom.inv_dual_null_point[0], &tmp, &leg_4[0]);
    fp_mul(&gluing->weil_codom.inv_dual_null_point[2], &tmp, &leg_4[1]);
    fp_mul(&gluing->weil_codom.inv_dual_null_point[5], &tmp, &leg_4[2]);

    fp_mul(&tmp, &leg_4[0], &leg_12[0]);
    fp_mul(&gluing->weil_codom.inv_dual_null_point[1], &tmp, &leg_8[1]);
    fp_mul(&gluing->weil_codom.inv_dual_null_point[7], &tmp, &leg_8[2]);

    fp_mul(&tmp, &leg_4[0], &leg_8[0]);
    fp_mul(&gluing->weil_codom.inv_dual_null_point[3], &tmp, &leg_12[1]);
    fp_mul(&gluing->weil_codom.inv_dual_null_point[9], &tmp, &leg_12[2]);

    fp_set_zero(&gluing->weil_codom.inv_dual_null_point[4]);
    fp_set_zero(&gluing->weil_codom.inv_dual_null_point[6]);
    fp_set_zero(&gluing->weil_codom.inv_dual_null_point[8]);

    // Compute the auxiliary point inv_T3

    // reminder of Weil-to-theta
    //{0, 1, 2, 3, 2, 4, 5, 6, 1, 7, 4, 8, 3, 8, 6, 9};
    fp_mul(&aux_point_compress[0], &HSK[0][0], &gluing->weil_codom.inv_dual_null_point[0]);
    fp_mul(&aux_point_compress[1], &HSK[0][1], &gluing->weil_codom.inv_dual_null_point[1]);
    fp_mul(&aux_point_compress[2], &HSK[0][2], &gluing->weil_codom.inv_dual_null_point[2]);
    fp_mul(&aux_point_compress[3], &HSK[0][3], &gluing->weil_codom.inv_dual_null_point[3]);
    fp_mul(&aux_point_compress[4], &HSK[0][8], &gluing->weil_codom.inv_dual_null_point[1]);
    fp_mul(&aux_point_compress[5], &HSK[0][9], &gluing->weil_codom.inv_dual_null_point[7]);
    fp_mul(&aux_point_compress[6],
           &HSK[0][10],
           &gluing->weil_codom.inv_dual_null_point[4]); // will be zero by the geometry
    fp_mul(&aux_point_compress[7], &HSK[0][15], &gluing->weil_codom.inv_dual_null_point[9]);

    // In order to compute aux_point_compress[6], we need to use the fourth point.
    //
    // aux_point_compress[6] =  aux_point_compress[0] * HSK[3][2] * inv_dual_null_point[2] / (HSK[3][8] *
    // inv_dual_null_point[1])
    //
    // But we can instead clear the denominator correction_den = kernel[3][8] * inv_dual_null_point[8]
    // with 7M instead of 1I

    fp_mul(&tmp, &HSK[3][2], &gluing->weil_codom.inv_dual_null_point[2]);
    fp_mul(&aux_point_compress[6], &aux_point_compress[0], &tmp); // will be zero by the geometry

    fp_mul(&tmp, &HSK[3][8], &gluing->weil_codom.inv_dual_null_point[1]); // denominator
    fp_mul(&aux_point_compress[0], &aux_point_compress[0], &tmp);
    fp_mul(&aux_point_compress[1], &aux_point_compress[1], &tmp);
    fp_mul(&aux_point_compress[2], &aux_point_compress[2], &tmp);
    fp_mul(&aux_point_compress[3], &aux_point_compress[3], &tmp);
    fp_mul(&aux_point_compress[4], &aux_point_compress[4], &tmp);
    fp_mul(&aux_point_compress[5], &aux_point_compress[5], &tmp);
    fp_mul(&aux_point_compress[7], &aux_point_compress[7], &tmp);

    // compute the projective inverse of aux_point_compress
    fp_proj_batched_inv(aux_point_compress, 8);

    // uncompress this point into the output.
    for (uint8_t i = 0; i < 2; i++) {
        for (uint8_t j = 0; j < 4; j++) {
            fp_copy(&gluing->inv_aux_point[8 * i + j], &aux_point_compress[4 * i + j]);
            fp_copy(&gluing->inv_aux_point[8 * i + j + 4], &aux_point_compress[4 * i + j]);
        }
    }
}

void
gluing_compute_weil_gluing(weil_gluing_t *gluing,
                           const theta_point_dim4_t *T3,
                           const theta_point_dim4_t *T4,
                           const theta_point_dim4_t *T34,
                           const theta_point_dim4_t *T3_p_2T4)
{

    theta_point_dim4_t HSK[4];
    theta_squared_dim4(&HSK[0], T3);
    theta_squared_dim4(&HSK[1], T4);
    theta_squared_dim4(&HSK[2], T34);
    theta_squared_dim4(&HSK[3], T3_p_2T4);

    for (uint8_t i = 0; i < 4; i++) {
        theta_hadamard_dim4(&HSK[i], &HSK[i]);
    }

    solve_gluing_HIIP(gluing, HSK);

    // Compute the arithmetic values on the gluing codomain. Will NOT be in dual mode as we cannot double in this model
    // have to use special function as we have zeros
    gluing_special_inv(&gluing->weil_precomp.inv_null_point_DBL,&gluing->weil_codom.inv_dual_null_point);

    theta_hadamard_weil(&gluing->weil_precomp.inv_null_point_DBL,
                        &gluing->weil_precomp.inv_null_point_DBL); // null_point is the theta null point

    theta_squared_weil(&gluing->weil_precomp.inv_HS_null_point_DBL, &gluing->weil_precomp.inv_null_point_DBL);
    theta_hadamard_weil(&gluing->weil_precomp.inv_HS_null_point_DBL, &gluing->weil_precomp.inv_HS_null_point_DBL);

    theta_invert_weil(&gluing->weil_precomp.inv_HS_null_point_DBL, &gluing->weil_precomp.inv_HS_null_point_DBL);

    theta_invert_weil(&gluing->weil_precomp.inv_null_point_DBL, &gluing->weil_precomp.inv_null_point_DBL);
    // gluing->weil_codom.arith_precomp = true;
}

void
gluing_evaluate_weil_gluing(theta_point_dim4_t *out,
                            const theta_point_dim4_t *P_p_T3,
                            const theta_point_dim4_t *P_m_T3,
                            const weil_gluing_t *gluing)
{
    theta_point_dim4_t temp;

    theta_dot_prod_dim4(&temp, P_p_T3, P_m_T3);
    theta_hadamard_dim4(&temp, &temp);
    theta_dot_prod_dim4(out, &temp, &gluing->inv_aux_point);
}

uint32_t
gluing_mike_compute(mike_gluing_t *gluing, const ec_curve_t *E, const mike_gluing_basis_t *basis)
{

    ec_jac_point_t T34_1_32, T34_2_32;
    theta_point_dim2_t V1, V2;
    theta_point_dim2_t T3_16[2], T4_16[2], T34_16[2];
    theta_point_dim2_t inv_T3_16[2], T3_16_p_T4_8[2];
    theta_point_dim4_t T3_8, T4_8, T34_8, T3_8_p_T4_4;
    uint32_t res;

    ec_copy_jac_point(&gluing->aux_point_1, &basis->T3_1_32);
    ec_copy_jac_point(&gluing->aux_point_2, &basis->T3_2_32);

    // compute auxiliary points needed for 3rd isogeny
    ec_ADD(&T34_1_32, &basis->T3_1_32, &basis->T4_1_32, E);
    ec_ADD(&T34_2_32, &basis->T3_2_32, &basis->T4_2_32, E);

    // Compute scholten gluing and verify the supersingularity of the public key
    res = schoten_gluing_compute_and_verify(&gluing->scholten_gluing, &basis->P_8, &basis->Q_8, E);

    // evaluate the basis of the second isogeny
    scholten_gluing_eval(&V2, &basis->P_16, &gluing->scholten_gluing, 1); // V2
    scholten_gluing_eval(&V1, &basis->Q_16, &gluing->scholten_gluing, 1); // V1

    // Compute the second isogeny
    gluing_diag_isogeny_compute(&gluing->phi_UV, &V1, &V2);

    // evaluate the basis of the 3rd isogeny throught the first isogeny
    scholten_gluing_eval(&T3_16[0], &basis->T3_1_32, &gluing->scholten_gluing, 1);
    scholten_gluing_eval(&T3_16[1], &basis->T3_2_32, &gluing->scholten_gluing, 1);

    scholten_gluing_eval(&T4_16[0], &basis->T4_1_32, &gluing->scholten_gluing, 1);
    scholten_gluing_eval(&T4_16[1], &basis->T4_2_32, &gluing->scholten_gluing, 1);

    scholten_gluing_eval(&T34_16[0], &T34_1_32, &gluing->scholten_gluing, 1);
    scholten_gluing_eval(&T34_16[1], &T34_2_32, &gluing->scholten_gluing, 1);

    // Apply the second isogeny
    gluing_diag_isogeny_eval(&T3_16[0], &T3_16[1], &T3_16[0], &T3_16[1], &gluing->phi_UV);
    gluing_diag_isogeny_eval(&T4_16[0], &T4_16[1], &T4_16[0], &T4_16[1], &gluing->phi_UV);
    gluing_diag_isogeny_eval(&T34_16[0], &T34_16[1], &T34_16[0], &T34_16[1], &gluing->phi_UV);

    // Compute 4th point of the gluing basis
    theta_invert_dim2(&inv_T3_16[0], &T3_16[0]);
    theta_invert_dim2(&inv_T3_16[1], &T3_16[1]);
    theta_dim2_diff_ADD(&T3_16_p_T4_8[0], &T34_16[0], &T4_16[0], &inv_T3_16[0], &gluing->phi_UV.S_U_precomp);
    theta_dim2_diff_ADD(&T3_16_p_T4_8[1], &T34_16[1], &T4_16[1], &inv_T3_16[1], &gluing->phi_UV.S_V_precomp);

    // Lift couple of dim 2 to dimension 4 points
    gluing_from_couple_dim2_to_dim4_compatible_with_isogeny(&T3_8, &T3_16[0], &T3_16[1]);
    gluing_from_couple_dim2_to_dim4_compatible_with_isogeny(&T4_8, &T4_16[0], &T4_16[1]);
    gluing_from_couple_dim2_to_dim4_compatible_with_isogeny(&T34_8, &T34_16[0], &T34_16[1]);
    gluing_from_couple_dim2_to_dim4_compatible_with_isogeny(&T3_8_p_T4_4, &T3_16_p_T4_8[0], &T3_16_p_T4_8[1]);

    // Compute gluing isogenies
    gluing_compute_weil_gluing(&gluing->weil_gluing, &T3_8, &T4_8, &T34_8, &T3_8_p_T4_4);
    return res;
}

void
gluing_mike_eval(theta_point_dim4_t *imgT,
                 const ec_jac_point_t *T1,
                 const ec_jac_point_t *T2,
                 const mike_gluing_t *gluing)
{

    theta_point_dim2_t T_p_dim2[2], T_m_dim2[2];
    theta_point_dim4_t T_p_dim4, T_m_dim4;

    ec_jac_point_t T1_p_aux_1, T1_m_aux_1;
    ec_jac_point_t T2_p_aux_2, T2_m_aux_2;
    ec_jac_point_t tmp;

    ec_ADD(&T1_p_aux_1, T1, &gluing->aux_point_1, &gluing->scholten_gluing.E); // T1 + aux1
    ec_jac_neg(&tmp, &gluing->aux_point_1);
    ec_ADD(&T1_m_aux_1, T1, &tmp, &gluing->scholten_gluing.E); // T1 - aux1

    ec_ADD(&T2_p_aux_2, T2, &gluing->aux_point_2, &gluing->scholten_gluing.E); // T2 + aux2
    ec_jac_neg(&tmp, &gluing->aux_point_2);
    ec_ADD(&T2_m_aux_2, T2, &tmp, &gluing->scholten_gluing.E); // T2 - aux2

    // Evaluate the first isogeny
    scholten_gluing_eval(&T_p_dim2[0], &T1_p_aux_1, &gluing->scholten_gluing, 1);
    scholten_gluing_eval(&T_p_dim2[1], &T2_p_aux_2, &gluing->scholten_gluing, 1);
    scholten_gluing_eval(&T_m_dim2[0], &T1_m_aux_1, &gluing->scholten_gluing, 1);
    scholten_gluing_eval(&T_m_dim2[1], &T2_m_aux_2, &gluing->scholten_gluing, 1);

    // Evaluate the second isogeny
    gluing_diag_isogeny_eval(&T_p_dim2[0], &T_p_dim2[1], &T_p_dim2[0], &T_p_dim2[1], &gluing->phi_UV);
    gluing_diag_isogeny_eval(&T_m_dim2[0], &T_m_dim2[1], &T_m_dim2[0], &T_m_dim2[1], &gluing->phi_UV);

    // Lift the theta dim 4
    gluing_from_couple_dim2_to_dim4_compatible_with_isogeny(&T_p_dim4, &T_p_dim2[0], &T_p_dim2[1]);
    gluing_from_couple_dim2_to_dim4_compatible_with_isogeny(&T_m_dim4, &T_m_dim2[0], &T_m_dim2[1]);

    gluing_evaluate_weil_gluing(imgT, &T_p_dim4, &T_m_dim4, &gluing->weil_gluing);

    // Apply hadamard here to be in theta model (because we cannot double in dual)
    theta_hadamard_dim4(imgT, imgT);
}
