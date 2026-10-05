#include <nike.h>

int
protocols_exchange(shared_t *out, const public_key_t *pk, const secret_key_t *sk)
{

    ec_jac_point_t PQ[2], PQ32[2];
    mike_chain_basis_t mike_basis;
    mike_gluing_basis_t gluing_basis;
    ec_curve_t Epk;
    ec_curve_init_from_A(&Epk, &pk->curveA);
    int ret = -1;
    uint8_t sk_bshift_32 = ((sk->x[0] >> 1) & 31);

    ret &= -((sk_bshift_32 & 3) == 1); 

    // Sample <P,Q> = E[2^(l+2)] a basis of E_pk
    ec_jac_basis_2f(&PQ[0], &PQ[1], &Epk, SECRETKEY_CHAIN_LENGTH + HD_MARGIN);

    // Compute chain basis
    mike_compute_basis(&mike_basis, &PQ[0], &PQ[1], sk->x, SECRETKEY_CHAIN_LENGTH + HD_MARGIN, &Epk);

    // Compute gluing basis
    ec_jac_dbl_iter(&PQ32[0], &PQ[0], SECRETKEY_CHAIN_LENGTH + HD_MARGIN - 5, &Epk, 2);

    gluing_basis_compute(&gluing_basis, &Epk, PQ32, sk_bshift_32);

    // Compute dim 4 chain
    ret &= mike_isogeny_chain_dim4(&out->invariant, &mike_basis, &gluing_basis, &Epk, SECRETKEY_CHAIN_LENGTH);

    return ret;
}
