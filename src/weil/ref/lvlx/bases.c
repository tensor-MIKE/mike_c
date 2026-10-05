#include <weil.h>
#include <assert.h>

void
mike_compute_basis(mike_chain_basis_t *basis,
                   const ec_jac_point_t *P,
                   const ec_jac_point_t *Q,
                   const digit_t *x_div2, // (x_sec - 1)/2
                   const int pow,
                   const ec_curve_t *E)
{

    ec_jac_point_t minus_P;

    ec_jac_MUL_shifted(&basis->T3_2, Q, x_div2, pow, E, 1);
    ec_ADD(&basis->T3_1, &basis->T3_2, Q, E);

    ec_jac_neg(&minus_P, P);
    ec_jac_MUL_shifted(&basis->T4_2, &minus_P, x_div2, pow, E, 1);
    ec_ADD(&basis->T4_1, &basis->T4_2, &minus_P, E);
}

void
gluing_basis_compute(mike_gluing_basis_t *basis,
                     const ec_curve_t *E,
                     const ec_jac_point_t *PQ_32,
                     const uint8_t x_bshift_mod32 // (x_sec - 1)/2
)
{
    ec_jac_point_t mP_32;

    ec_DBL(&basis->P_16, &PQ_32[0], E);
    ec_DBL(&basis->Q_16, &PQ_32[1], E);

    ec_DBL(&basis->P_8, &basis->P_16, E);
    ec_DBL(&basis->Q_8, &basis->Q_16, E);

    digit_t x_bshift_mod32_digit = x_bshift_mod32;

    ec_jac_MUL(&basis->T3_2_32, &PQ_32[1], &x_bshift_mod32_digit, 5, E); //[(sec - 1)/2]Q_32
    ec_ADD(&basis->T3_1_32, &basis->T3_2_32, &PQ_32[1], E);              //[(sec + 1)/2]Q_32

    ec_jac_neg(&mP_32, &PQ_32[0]); // P_2 = -P32

    ec_jac_MUL(&basis->T4_2_32, &mP_32, &x_bshift_mod32_digit, 5, E); //-[(sec - 1)/2]P_32
    ec_ADD(&basis->T4_1_32, &basis->T4_2_32, &mP_32, E);              //-[(sec + 1)/2]P_32
}
