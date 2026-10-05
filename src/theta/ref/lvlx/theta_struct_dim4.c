#include <theta.h>

// Faster doubling using Weil special curves.
// This does the precomputation in theta dual mode.
void
theta_struct_weil_arith_precomp(theta_struct_weil_precomp_t *theta_precomp, const theta_struct_weil_t *theta_struct)
{
    // Cost: 2*(3*10-4)M + 16 S + 64 a = 52 M + 16 S + 64 a
    theta_copy_weil(&theta_precomp->inv_null_point_DBL, &theta_struct->inv_dual_null_point);

    theta_invert_weil(&theta_precomp->inv_HS_null_point_DBL, &theta_struct->inv_dual_null_point);
    theta_squared_weil(&theta_precomp->inv_HS_null_point_DBL, &theta_precomp->inv_HS_null_point_DBL);
    theta_hadamard_weil(&theta_precomp->inv_HS_null_point_DBL, &theta_precomp->inv_HS_null_point_DBL);

    theta_invert_weil(&theta_precomp->inv_HS_null_point_DBL, &theta_precomp->inv_HS_null_point_DBL);

}

void
theta_DBL_weil(theta_point_dim4_t *out, const theta_point_dim4_t *in, const theta_struct_weil_precomp_t *theta_struct)
{
    // Cost (without precomputation): 32 M + 32 S + 128 a
    theta_squared_dim4(out, in);
    theta_hadamard_dim4(out, out);
    theta_squared_dim4(out, out);
    theta_dot_prod_weil(out, out, &theta_struct->inv_HS_null_point_DBL);
    theta_hadamard_dim4(out, out);
    theta_dot_prod_weil(out, out, &theta_struct->inv_null_point_DBL);
}

void
theta_DBL_iter_weil(theta_point_dim4_t *out, const theta_point_dim4_t *in, const theta_struct_weil_precomp_t *theta_struct, const int n)
{

    theta_copy_theta_dim4(out, in);
    for (int i = 0; i < n; i++) {
        theta_DBL_weil(out, out, theta_struct);
    }
}
