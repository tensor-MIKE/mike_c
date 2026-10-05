#include <theta.h>

void
theta_struct_dim2_arith_precomp(theta_struct_dim2_precomp_t *theta_precomp, const theta_struct_dim2_t *theta_struct)
{
    // Assume that both inv_dual_null_point and dual_null_point are computed during codomain compute
    // Furthermore, we are in dual mode
    // Cost: 6M + 4S + 8a
    theta_copy_dim2(&theta_precomp->inv_null_point_DBL, &theta_struct->inv_dual_null_point);

    theta_squared_dim2(&theta_precomp->inv_HS_null_point_DBL, &theta_struct->dual_null_point);
    theta_hadamard_dim2(&theta_precomp->inv_HS_null_point_DBL, &theta_precomp->inv_HS_null_point_DBL);

    theta_invert_dim2(&theta_precomp->inv_HS_null_point_DBL, &theta_precomp->inv_HS_null_point_DBL);
}

void
theta_dim2_diff_ADD(theta_point_dim2_t *out,
                    const theta_point_dim2_t *P,
                    const theta_point_dim2_t *Q,
                    const theta_point_dim2_t *inv_PmQ,
                    const theta_struct_dim2_precomp_t *theta_struct)
{

    theta_point_dim2_t tmp1;

    theta_squared_dim2(out, P);
    theta_hadamard_dim2(out, out);
    theta_squared_dim2(&tmp1, Q);
    theta_hadamard_dim2(&tmp1, &tmp1);
    theta_dot_prod_dim2(out, out, &tmp1);
    theta_dot_prod_dim2(out, out, &theta_struct->inv_HS_null_point_DBL);
    theta_hadamard_dim2(out, out);
    theta_dot_prod_dim2(out, out, inv_PmQ);
}

void
theta_dim2_DBL(theta_point_dim2_t *out, const theta_point_dim2_t *in, const theta_struct_dim2_precomp_t *theta_struct)
{
    // Cost (without precomputation): 8 M + 8 S + 16 a
    theta_squared_dim2(out, in);
    theta_hadamard_dim2(out, out);
    theta_squared_dim2(out, out);
    theta_dot_prod_dim2(out, out, &theta_struct->inv_HS_null_point_DBL);
    theta_hadamard_dim2(out, out);
    theta_dot_prod_dim2(out, out, &theta_struct->inv_null_point_DBL);
}

void
theta_dim2_DBL_iter(theta_point_dim2_t *out,
                    const theta_point_dim2_t *in,
                    const theta_struct_dim2_precomp_t *theta_struct,
                    const int n)
{
    theta_copy_dim2(out, in);
    for (int i = 0; i < n; i++) {
        theta_dim2_DBL(out, out, theta_struct);
    }
}
