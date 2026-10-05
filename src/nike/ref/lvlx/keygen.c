#include <nike.h>
#include <e0_basis.h>
#include <ec_isog.h>
#include "rng.h"

static void
nike_use_precomputed(ec_xz_point_t *P, ec_xz_point_t *Q, ec_xz_point_t *PmQ, ec_curve_t *curve)
{

    // set curve to curve of A-invariant A0
    fp2_t Afp2;
    fp_set_zero(&Afp2.im);
    fp_copy(&Afp2.re, &A0);
    ec_curve_init_from_A(curve, &Afp2);
    ec_compute_A24(curve);
    // set sk kernel points using x from preomputation
    fp2_copy(&P->x, &BASIS_sk_PX);
    fp2_set_one(&P->z);
    fp2_copy(&Q->x, &BASIS_sk_QX);
    fp2_set_one(&Q->z);
    fp2_copy(&PmQ->x, &BASIS_sk_PMQX);
    fp2_set_one(&PmQ->z);
    
#ifndef NDEBUG
    int res = 0;
    ec_jac_point_t Pc, Qc, pc, qc;
    ec_xz_to_jac(&Pc, P, curve);
    ec_xz_to_jac(&Qc, Q, curve);
    ec_jac_dbl_iter(&pc, &Pc, SECRETKEY_CHAIN_LENGTH + HD_MARGIN - 1, curve, 1);
    ec_jac_dbl_iter(&qc, &Qc, SECRETKEY_CHAIN_LENGTH + HD_MARGIN - 1, curve, 1);
    res = res | ec_jac_is_zero(&pc);
    res = res | ec_jac_is_zero(&qc);
    ec_DBL(&Pc, &pc, curve);
    ec_DBL(&Qc, &qc, curve);
    res = res | !ec_jac_is_zero(&Pc);
    res = res | !ec_jac_is_zero(&Qc);
    assert(!res);
#endif
}

int
protocols_keygen(public_key_t *pk, secret_key_t *sk)
{
    // sample the secret key 
    unsigned char xbytes[SECRETKEY_BYTES];
    randombytes(xbytes, SECRETKEY_BYTES);
    // shorten to bitlength and set fixed bits
    xbytes[0] = xbytes[0] - (xbytes[0] & 7) + 3;
    xbytes[SECRETKEY_BYTES - 1] = xbytes[SECRETKEY_BYTES - 1] % (1 << (SECRETKEY_BITS & 7));
    // little-endian bytewise to little-endian wordwise
    secret_key_from_bytes(sk, xbytes);

    //  compute kernel isogeny
    ec_isog_even_t iso;
    ec_xz_point_t P, Q, PmQ;
    iso.length = SECRETKEY_CHAIN_LENGTH + HD_MARGIN;
    nike_use_precomputed(&P, &Q, &PmQ, &iso.curve);
    int res = ec_ladder3pt(&iso.kernel, &P, &Q, &PmQ, sk->x, SECRETKEY_BITS, &iso.curve);

#ifndef NDEBUG
    // test kernel order
    int resdbg = 0;
    ec_jac_point_t Pc, pc;
    ec_xz_to_jac(&Pc, &iso.kernel, &iso.curve);
    ec_jac_dbl_iter(&pc, &Pc, SECRETKEY_CHAIN_LENGTH + HD_MARGIN - 1, &iso.curve, 1);
    resdbg = resdbg | ec_jac_is_zero(&pc);
    ec_DBL(&Pc, &pc, &iso.curve);
    resdbg = resdbg | !ec_jac_is_zero(&Pc);
    resdbg = resdbg | !(SECRETKEY_CHAIN_LENGTH + HD_MARGIN == iso.length);
    assert(!resdbg);
#endif
    // evaluate isogeny
    res = res & !ec_iso_isogeny_2chain(&pk->curveA, &iso, 2);
#ifndef NDEBUG
    ec_curve_t curve;
    ec_curve_init_from_A(&curve, &pk->curveA);
    // test -sk gives same pk
    fp2_t cdmA;
    ec_xADD(&PmQ, &P, &Q, &PmQ);
    ec_ladder3pt(&iso.kernel, &P, &Q, &PmQ, sk->x, SECRETKEY_BITS, &iso.curve);
    ec_iso_isogeny_2chain(&cdmA, &iso, 2);
    assert(fp2_is_equal(&cdmA, &pk->curveA));
#endif

    return res;
}
