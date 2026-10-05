
#include <theta.h>
#include <stdint.h>


void isogeny_dim2_compute_codomain(theta_struct_dim2_t *codomain,
                               const theta_point_dim2_t *ker
                                )
{


   theta_point_dim2_t HSK_8[2];
    uint8_t i;
    for (i = 0; i < 2; i++) {
        theta_squared_dim2(&HSK_8[i], &ker[i]);
        theta_hadamard_dim2(&HSK_8[i], &HSK_8[i]);
    }

    fp_t t1, t2; 

    fp_mul(&t1, &HSK_8[0][0], &HSK_8[1][1]);
    fp_mul(&t2, &HSK_8[0][1], &HSK_8[1][0]);
    fp_mul(&(codomain->dual_null_point)[0], &HSK_8[1][0], &t1);
    fp_mul(&(codomain->dual_null_point)[1], &HSK_8[1][1], &t2);
    fp_mul(&(codomain->dual_null_point)[2], &HSK_8[1][2], &t1);
    fp_mul(&(codomain->dual_null_point)[3], &HSK_8[1][3], &t2);
    fp_t t3;
    fp_mul(&t3, &HSK_8[1][2], &HSK_8[1][3]);
    fp_mul(&(codomain->inv_dual_null_point)[0], &t3, &HSK_8[0][1]);
    fp_mul(&(codomain->inv_dual_null_point)[1], &t3, &HSK_8[0][0]);
    fp_copy(&(codomain->inv_dual_null_point)[2], &(codomain->dual_null_point)[3]);
    fp_copy(&(codomain->inv_dual_null_point)[3], &(codomain->dual_null_point)[2]);

}


void isogeny_dim2_eval(theta_point_dim2_t *imP, const theta_point_dim2_t *P, const theta_struct_dim2_t *codomain)
{

    theta_squared_dim2(imP,P);
    theta_hadamard_dim2(imP,imP);
    theta_dot_prod_dim2(imP, imP, &(codomain->inv_dual_null_point));

}
