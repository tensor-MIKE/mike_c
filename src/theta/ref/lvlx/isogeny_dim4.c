
#include <theta.h>

const uint8_t theta_to_weil[16] = { 0, 1, 2, 3, 2, 4, 5, 6, 1, 7, 4, 8, 3, 8, 6, 9 };
const uint8_t weil_to_theta[10] = { 0, 1, 2, 3, 5, 6, 7, 9, 11, 15 };

// Draw me a rabbit !
static void
solve_HIIP_weil(theta_point_weil_t *inv_dual, const theta_point_dim4_t *HSK_8)
{

    fp_t red[8], blue[8], tmp;

    // compute red
    fp_copy(&red[0], &HSK_8[0][12]);
    fp_mul(&red[1], &red[0], &HSK_8[0][5]);
    fp_mul(&red[2], &red[1], &HSK_8[1][1]);
    fp_mul(&red[3], &red[2], &HSK_8[0][9]);
    fp_mul(&red[4], &red[3], &HSK_8[0][11]);
    fp_mul(&red[5], &red[4], &HSK_8[1][15]);
    fp_mul(&red[6], &red[5], &HSK_8[1][14]);
    fp_mul(&red[7], &red[6], &HSK_8[0][6]);

    // compute blue
    fp_copy(&blue[7], &HSK_8[0][0]);
    fp_mul(&blue[6], &blue[7], &HSK_8[0][2]);
    fp_mul(&blue[5], &blue[6], &HSK_8[1][6]);
    fp_mul(&blue[4], &blue[5], &HSK_8[1][7]);
    fp_mul(&blue[3], &blue[4], &HSK_8[0][15]);
    fp_mul(&blue[2], &blue[3], &HSK_8[0][13]);
    fp_mul(&blue[1], &blue[2], &HSK_8[1][9]);
    fp_mul(&blue[0], &blue[1], &HSK_8[0][1]);

    // rabbit part
    fp_mul(&(*inv_dual)[0], &red[7], &HSK_8[0][4]);
    fp_mul(&tmp, &HSK_8[0][8], &HSK_8[0][5]);
    fp_mul(&(*inv_dual)[3], &tmp, &blue[1]);

    // merge the rest
    //{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}
    //{0, 1, 2, 3, 2, 4, 5, 6, 1, 7,  4,  8,  3,  8,  6,  9}
    fp_mul(&(*inv_dual)[4], &red[0], &blue[0]);
    fp_mul(&(*inv_dual)[1], &red[1], &blue[1]);
    fp_mul(&(*inv_dual)[7], &red[2], &blue[2]);
    fp_mul(&(*inv_dual)[8], &red[3], &blue[3]);
    fp_mul(&(*inv_dual)[9], &red[4], &blue[4]);
    fp_mul(&(*inv_dual)[6], &red[5], &blue[5]);
    fp_mul(&(*inv_dual)[5], &red[6], &blue[6]);
    fp_mul(&(*inv_dual)[2], &red[7], &blue[7]);
}

void
isogeny_weil_compute_codomain(theta_struct_weil_t *codomain, const theta_point_dim4_t *ker, const bool dual_domain)
{
    theta_point_dim4_t HSK_8[2];
    uint8_t i;

    for (i = 0; i < 2; i++) {
        if (dual_domain)
            theta_hadamard_dim4(&HSK_8[i], &ker[i]);
        else
            theta_copy_theta_dim4(&HSK_8[i], &ker[i]);

        theta_squared_dim4(&HSK_8[i], &HSK_8[i]);
        theta_hadamard_dim4(&HSK_8[i],&HSK_8[i]);
    }

    solve_HIIP_weil(&codomain->inv_dual_null_point, HSK_8);
}

void
isogeny_weil_eval(theta_point_dim4_t *imP,
               const theta_point_dim4_t *P,
               const theta_struct_weil_t *codomain,
               const bool dual_domain)
{
    theta_point_dim4_t point;

    if (dual_domain)
        theta_hadamard_dim4(&point, P);
    else
        theta_copy_theta_dim4(&point, P);

    // Cost 16 M + 16 S + 128 a
    theta_squared_dim4(&point, &point);
    theta_hadamard_dim4(&point, &point);
    theta_dot_prod_weil(imP, &point, &(codomain->inv_dual_null_point));

}
