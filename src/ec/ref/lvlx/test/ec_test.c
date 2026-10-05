#include <assert.h>
#include <stdio.h>
#include <inttypes.h>

#include "test_extras.h"
#include <ec.h>
#include <ec_params.h>
#include <constants.h>
#include <rng.h>
#include <bench_test_arguments.h>

/******************************
Test functions
******************************/

void
ec_jac_normalize(ec_jac_point_t *P)
{
    fp2_t t0, t1;

    fp2_copy(&t0, &P->z);       // t0=z
    fp2_inv(&t0);               // t0=1/z
    fp2_sqr(&t1, &t0);          // t1=1/z^2
    fp2_mul(&P->x, &P->x, &t1); // x=x/z^2
    fp2_mul(&t1, &t1, &t0);     // t1=1/z^3
    fp2_mul(&P->y, &P->y, &t1); // y=y/z^3
    fp2_set_one(&P->z);         // z=1
}

void
ec_ADD_test(ec_jac_point_t *R, const ec_jac_point_t *P, const ec_jac_point_t *Q, const ec_curve_t *E)
{
    fp2_t t0, t1, t2, t3, t4, t5, t6, t7, t8, dx, dy;

    fp2_sqr(&t0, &P->z);      // t0 = z1^2
    fp2_sqr(&t1, &Q->z);      // t1 = z2^2
    fp2_mul(&t2, &P->x, &t1); // t2 = x1*z2^2
    fp2_mul(&t3, &Q->x, &t0); // t3 = x2*z1^2
    fp2_sub(&dx, &t3, &t2);   // dx = x2*z1^2 - x1*z2^2

    fp2_mul(&t0, &t0, &P->z); // t0 = z1^3
    fp2_mul(&t1, &t1, &Q->z); // t1 = z2^3
    fp2_mul(&t4, &P->y, &t1); // t4 = y1*z2^3
    fp2_mul(&t5, &Q->y, &t0); // t5 = y2*z1^3
    fp2_sub(&dy, &t5, &t4);   // dy = y2*z1^3 - y1*z2^3

    fp2_mul(&t6, &P->z, &Q->z); // t6 = z1*z2
    fp2_mul(&R->z, &t6, &dx);   // z3 = z1*z2*dx

    fp2_sqr(&t6, &t6);        // t6 = (z1*z2)^2
    fp2_mul(&t6, &E->A, &t6); // t6 = A*(z1*z2)^2
    fp2_add(&t6, &t6, &t2);   // t6 = A*(z1*z2)^2 + x1*z2^2
    fp2_add(&t6, &t6, &t3);   // t6 = A*(z1*z2)^2 + x1*z2^2 + x2*z1^2
    fp2_sqr(&t7, &dx);        // t7 = (dx)^2
    fp2_sqr(&t8, &dy);        // t8 = (dy)^2
    fp2_mul(&t6, &t7, &t6);   // t6 = (dx)^2*(A*(z1*z2)^2 + x1*z2^2 + x2*z1^2)
    fp2_sub(&R->x, &t8, &t6); // x3 = (dy)^2 - (dx)^2*(A*(z1*z2)^2 + x1*z2^2 + x2*z1^2)

    fp2_mul(&t6, &t2, &t7);   // t6 = x1*z2^2*(dx)^2
    fp2_sub(&t6, &t6, &R->x); // t6 = x1*z2^2*(dx)^2 - x3
    fp2_mul(&t6, &t6, &dy);   // t6 = dy*(x1*z2^2*(dx)^2 - x3)
    fp2_mul(&t7, &t7, &dx);   // t7 = (dx)^3
    fp2_mul(&t8, &t4, &t7);   // t8 = y1*z2^3*(dx)^3
    fp2_sub(&R->y, &t6, &t8); // y3 = dy*(x1*z2^2*(dx)^2 - x3) - y1*z2^3*(dx)^3
}

uint32_t
ec_jac_is_on_curve(const ec_jac_point_t *P, const ec_curve_t *E)
{
    fp2_t t0, t1, t2, x, y2;

    fp2_copy(&t0, &P->z);
    fp2_inv(&t0);             // t0 =1/Z
    fp2_sqr(&t1, &t0);        // t1 = 1/Z^2
    fp2_mul(&x, &P->x, &t1);  // x = X/Z^2
    fp2_mul(&t1, &t1, &t0);   // t1 = 1/Z^3
    fp2_mul(&y2, &P->y, &t1); // y = Y/Z^3
    fp2_sqr(&y2, &y2);        // y2 = y^2

    fp2_sqr(&t0, &x);        // t0 = x^2
    fp2_mul(&t1, &E->A, &x); // t1 = A*x
    fp2_set_one(&t2);        // t2 = 1
    fp2_add(&t2, &t2, &t1);  // t2 = A*x + 1
    fp2_add(&t2, &t2, &t0);  // t2 = x^2 + A*x +1
    fp2_mul(&t2, &t2, &x);   // t2 = x*(x^2 + A*x +1)

    return fp2_is_equal(&t2, &y2);
}

int
ec_test_jacobian(ec_curve_t *curve, unsigned int Ntest)
{
    unsigned int i;

    ec_xz_point_t P, Q, X, P1;
    ec_jac_point_t R, S, T, U, jac_zero;
    ec_jac_ws_point_t Rw;
    ec_ws_curve_t Ews;
    ec_compute_ws(&Ews, curve);

    // set jac_zero to 0
    fp2_set_zero(&jac_zero.x);
    fp2_set_one(&jac_zero.y);
    fp2_set_zero(&jac_zero.z);

    for (i = 0; i < Ntest; i++) {
        ec_xz_random_test(&P, curve);
        ec_xz_random_test(&Q, curve);
        ec_xz_random_test(&X, curve);

        // Convert to Jacobian coordinates.
        ec_xz_to_jac(&R, &P, curve);
        ec_xz_to_jac(&S, &Q, curve);
        ec_xz_to_jac(&T, &X, curve);

        // Convert back to (X:Y:Z)
        ec_jac_to_xz(&P1, &R);
        if (!ec_xz_is_equal(&P1, &P)) {
            printf("Failed Montgomery (X:Y:Z) <--> Jacobian conversion\n");
            return 1;
        }

        ec_ADD(&R, &jac_zero, &jac_zero, curve);
        if (!ec_jac_is_equal(&R, &jac_zero)) {
            printf("Failed 0 + 0 = 0 in jac\n");
            return 1;
        }

        ec_DBL(&R, &jac_zero, curve);
        if (!ec_jac_is_equal(&R, &jac_zero)) {
            printf("Failed 2*0 = 0 in jac\n");
            return 1;
        }

        ec_jac_neg(&R, &S);
        ec_ADD(&R, &S, &R, curve);
        if (!ec_jac_is_equal(&R, &jac_zero)) {
            printf("Failed P - P = 0 in jac\n");
            return 1;
        }

        ec_ADD(&R, &S, &jac_zero, curve);
        if (!ec_jac_is_equal(&R, &S)) {
            printf("Failed P + 0 = P in jac\n");
            return 1;
        }

        ec_ADD(&R, &jac_zero, &S, curve);
        if (!ec_jac_is_equal(&R, &S)) {
            printf("Failed P + 0 = P in jac\n");
            return 1;
        }

        ec_ADD(&R, &S, &jac_zero, curve);
        if (!ec_jac_is_equal(&R, &S)) {
            printf("Failed 0 + P = P in jac\n");
            return 1;
        }

        ec_DBL(&R, &S, curve);
        ec_ADD(&U, &S, &S, curve);
        if (!ec_jac_is_equal(&R, &U)) {
            printf("Failed P + P = 2*P in jac\n");
            return 1;
        }

        ec_ADD(&R, &T, &S, curve);
        ec_ADD(&T, &S, &T, curve);
        if (!ec_jac_is_equal(&R, &T)) {
            printf("Failed P + Q = Q + P in jac\n");
            return 1;
        }

        ec_ADD(&R, &T, &S, curve);
        ec_ADD(&U, &S, &T, curve);
        if (!ec_jac_is_equal(&R, &U)) {
            printf("Failed P + Q = Q + P in jac\n");
            return 1;
        }

        // Double R to make it different than (T + S).
        ec_DBL(&R, &R, curve);
        ec_ADD(&U, &S, &T, curve);
        ec_ADD(&U, &U, &R, curve);
        ec_ADD(&R, &R, &T, curve);
        ec_ADD(&R, &R, &S, curve);
        if (!ec_jac_is_equal(&R, &U)) {
            printf("Failed (P + Q) + R = P + (Q + R) in jac\n");
            return 1;
        }

        ec_ADD(&R, &S, &T, curve);
        ec_DBL(&R, &R, curve);
        ec_DBL(&S, &S, curve);
        ec_DBL(&T, &T, curve);
        ec_ADD(&U, &S, &T, curve);
        if (!ec_jac_is_equal(&R, &U)) {
            printf("Failed 2*(P + Q) = 2*P + 2*Q in jac\n");
            return 1;
        }

        ec_jac_to_ws(&Rw, &jac_zero, &Ews);
        ec_ws_to_jac(&R, &Rw, &Ews);
        if (!ec_jac_is_equal(&R, &jac_zero)) {
            printf("Failed converting to Weierstrass 0\n");
            return 1;
        }

        ec_jac_to_ws(&Rw, &S, &Ews);
        ec_ws_to_jac(&R, &Rw, &Ews);
        if (!ec_jac_is_equal(&S, &R)) {
            printf("Failed converting to Weierstrass S\n");
            return 1;
        }

        ec_DBL(&S, &S, curve);
        ec_jac_to_ws(&Rw, &S, &Ews);
        ec_ws_to_jac(&R, &Rw, &Ews);
        if (!ec_jac_is_equal(&S, &R)) {
            printf("Failed converting to Weierstrass 2*S\n");
            return 1;
        }

        ec_jac_to_ws(&Rw, &jac_zero, &Ews);
        ec_DBL_ws(&Rw, &Rw);
        ec_ws_to_jac(&R, &Rw, &Ews);
        if (!ec_jac_is_equal(&R, &jac_zero)) {
            printf("Failed 2*0 = 0 in Weierstrass\n");
            return 1;
        }

        ec_jac_to_ws(&Rw, &S, &Ews);
        ec_DBL_ws(&Rw, &Rw);
        ec_ws_to_jac(&R, &Rw, &Ews);
        ec_DBL(&S, &S, curve);
        if (!ec_jac_is_equal(&S, &R)) {
            printf("Failed doubling in Weierstrass\n");
            return 1;
        }
    }
    printf("Jacobian tests....................................................... PASSED\n");
    return 0;
}

int
ec_test_xDBL_xADD(ec_curve_t *curve, unsigned int Ntest)
{
    ec_xz_point_t P, Q;
    ec_jac_point_t jP, jQ, jPQ;
    ec_xz_point_t xP, xQ, xPQ, xR1, xR2;
    for (unsigned int i = 0; i < Ntest; i++) {
        ec_xz_random_test(&P, curve);
        ec_xz_random_test(&Q, curve);

        ec_xz_to_jac(&jP, &P, curve);
        ec_xz_to_jac(&jQ, &Q, curve);
        ec_jac_neg(&jQ, &jQ);

        ec_ADD(&jPQ, &jP, &jQ, curve);
        ec_jac_neg(&jQ, &jQ);

        ec_copy_xz_point(&xP, &P);
        ec_copy_xz_point(&xQ, &Q);
        ec_jac_to_xz(&xPQ, &jPQ);

        // 2(P + Q) = 2P + 2Q
        ec_xADD(&xR1, &xP, &xQ, &xPQ);
        ec_xDBL(&xR1, &xR1, curve);
        ec_xDBL(&xP, &xP, curve);
        ec_xDBL(&xQ, &xQ, curve);
        ec_xDBL(&xPQ, &xPQ, curve);
        ec_xADD(&xR2, &xP, &xQ, &xPQ);
        if (!ec_xz_is_equal(&xR1, &xR2)) {
            printf("Failed 2(P + Q) = 2P + 2Q\n");
            return 1;
        }

        // (P+Q) + (P-Q) = 2P
        ec_xADD(&xR1, &xP, &xQ, &xPQ);
        ec_xDBL(&xQ, &xQ, curve);
        ec_xADD(&xR1, &xR1, &xPQ, &xQ);
        ec_xDBL(&xP, &xP, curve);
        if (!ec_xz_is_equal(&xR1, &xP)) {
            printf("Failed (P+Q) + (P-Q) = 2P\n");
            return 1;
        }
    }
    printf("xDBL and xADD tests.................................................. PASSED\n");
    return 0;
}

int
ec_test_xDBLADD(ec_curve_t *curve, unsigned int Ntest)
{
    ec_xz_point_t P, Q;
    ec_jac_point_t jP, jQ, jPQ;
    ec_xz_point_t xP, xQ, xPQ, xR1, xR2;

    for (unsigned int i = 0; i < Ntest; i++) {
        ec_xz_random_test(&P, curve);
        ec_xz_random_test(&Q, curve);

        ec_xz_to_jac(&jP, &P, curve);
        ec_xz_to_jac(&jQ, &Q, curve);
        ec_jac_neg(&jQ, &jQ);

        ec_ADD(&jPQ, &jP, &jQ, curve);
        ec_jac_neg(&jQ, &jQ);

        ec_copy_xz_point(&xP, &P);
        ec_copy_xz_point(&xQ, &Q);
        ec_jac_to_xz(&xPQ, &jPQ);

        ec_xDBLADD(&xR1, &xR2, &xP, &xQ, &xPQ, curve);
        ec_xADD(&xPQ, &xP, &xQ, &xPQ);
        if (!ec_xz_is_equal(&xR2, &xPQ)) {
            printf("Failed addition in xDBLADD\n");
            return 1;
        }
        ec_xDBL(&xP, &xP, curve);
        if (!ec_xz_is_equal(&xR1, &xP)) {
            printf("Failed doubling in xDBLADD\n");
            return 1;
        }
    }
    printf("xDBLADD tests........................................................ PASSED\n");
    return 0;
}

int
ec_test_zero_identities(ec_curve_t *curve, unsigned int Ntest)
{
    unsigned int i;

    ec_xz_point_t P, Q, R, ec_zero;

    // set ec_zero to 0
    fp2_set_one(&ec_zero.x);
    fp2_set_zero(&ec_zero.z);

    assert(ec_xz_is_zero(&ec_zero));

    for (i = 0; i < Ntest; i++) {
        ec_xz_random_test(&P, curve);

        ec_xADD(&R, &ec_zero, &ec_zero, &ec_zero);
        if (!ec_xz_is_zero(&R)) {
            printf("Failed 0 + 0 = 0\n");
            return 1;
        }

        ec_xDBL(&R, &P, curve);
        ec_xADD(&R, &P, &P, &R);
        if (!ec_xz_is_zero(&R)) {
            printf("Failed P - P = 0\n");
            return 1;
        }

        ec_xDBL(&R, &ec_zero, curve);
        if (!ec_xz_is_zero(&R)) {
            printf("Failed 2*0 = 0\n");
            return 1;
        }

        ec_xADD(&R, &P, &ec_zero, &P);
        if (!ec_xz_is_equal(&R, &P)) {
            printf("Failed P + 0 = P\n");
            return 1;
        }
        ec_xADD(&R, &ec_zero, &P, &P);
        if (!ec_xz_is_equal(&R, &P)) {
            printf("Failed P + 0 = P\n");
            return 1;
        }

        ec_xDBLADD(&R, &Q, &P, &ec_zero, &P, curve);
        if (!ec_xz_is_equal(&Q, &P)) {
            printf("Failed P + 0 = P in xDBLADD\n");
            return 1;
        }
        ec_xDBLADD(&R, &Q, &ec_zero, &P, &P, curve);
        if (!ec_xz_is_equal(&Q, &P)) {
            printf("Failed P + 0 = P in xDBLADD\n");
            return 1;
        }
        if (!ec_xz_is_zero(&R)) {
            printf("Failed 2*0 = 0 in xDBLADD\n");
            return 1;
        }
    }
    printf("xz zero identities................................................... PASSED\n");
    return 0;
}

int
ec_test_mul(ec_curve_t *curve, unsigned int Ntest)
{
    digit_t pp1[NWORDS_ORDER];
    digit_t carry = 1;

    //set pp1 to p+1
    for (int i = 0; i < NWORDS_ORDER; i++) {
        pp1[i] = CHARACTERISTIC[i] + carry;
        if ((carry != 0) && ((pp1[i]) == 0))
            carry = 1;
        else
            carry = 0;
    }
    int kbits = ec_n_bits(pp1, NWORDS_ORDER);
    int kbits_c = ec_n_bits(TORSION_ODD, 1);

    ec_xz_point_t xP, xQ, xR;
    ec_jac_point_t jP, jQ;

    for (unsigned int i = 0; i < Ntest; i++) {
        ec_xz_random_test(&xP, curve);

        ec_xMUL(&xQ, &xR, &xP, pp1, kbits, curve);

        if (!ec_xz_is_zero(&xQ)) {
            printf("Failed (p+1)*P = 0 in xMUL\n");
            return 1;
        }
        if (!ec_xz_is_equal(&xR, &xP)) {
            printf("Failed (p+2)*P = P in xMUL\n");
            return 1;
        }

        ec_xMUL(&xQ, &xR, &xP, TORSION_ODD, kbits_c, curve);
        ec_xz_to_jac(&jQ, &xQ, curve);
        ec_jac_dbl_iter(&jQ, &jQ, TORSION_EVEN_POWER, curve, 1);
        if (!ec_jac_is_zero(&jQ)) {
            printf("Failed (p+1)*P = 0 in ec_MUL and jac_dbl_iter\n");
            return 1;
        }

        ec_xz_to_jac(&jP, &xP, curve);

        ec_jac_MUL(&jQ, &jP, pp1, kbits, curve);
        if (!ec_jac_is_zero(&jQ)) {
            printf("Failed (p+1)*P = 0 in jac_MUL\n");
            return 1;
        }
    }
    printf("Multiplication tests................................................. PASSED\n");
    return 0;
}

int
ec_test_basis(ec_curve_t *curve)
{

    ec_jac_point_t P, Q;

    ec_jac_basis_2f(&P, &Q, curve, TORSION_EVEN_POWER);

    if (!ec_jac_is_on_curve(&P, curve)) {
        printf("P not on curve\n");
        return 1;
    }
    if (!ec_jac_is_on_curve(&Q, curve)) {
        printf("Q not on curve\n");
        return 1;
    }

    ec_jac_dbl_iter(&P, &P, TORSION_EVEN_POWER - 1, curve, 1);
    if (ec_jac_is_zero(&P)) {
        printf("Wrong order: 2^(e-1)*P = 0\n");
        return 1;
    }

    ec_DBL(&P, &P, curve);

    if (!ec_jac_is_zero(&P)) {
        printf("Wrong order: 2^e*P != 0\n");
        return 1;
    }

    ec_jac_dbl_iter(&Q, &Q, TORSION_EVEN_POWER - 1, curve, 1);
    if (ec_jac_is_zero(&Q)) {
        printf("Wrong order: 2^(e-1)*Q = 0\n");
        return 1;
    }

    ec_DBL(&Q, &Q, curve);
    if (!ec_jac_is_zero(&Q)) {
        printf("Wrong order: 2^e*Q != 0\n");
        return 1;
    }

    printf("Basis generation tests............................................... PASSED\n");
    return 0;
}

int
ec_test_3ptladder(ec_curve_t *curve, unsigned int Ntest)
{
    ec_jac_point_t Bj[3];
    ec_xz_point_t B[3];
    ec_xz_point_t res, cmp;
    ec_jac_point_t jcmp, tmp;
    digit_t exp[NWORDS_ORDER];
    digit_t rand;

    for (unsigned int iter = 0; iter < NWORDS_ORDER; iter++) {
        for (int k = 0; k < NWORDS_ORDER; k++) {
            randombytes((unsigned char *)&rand, sizeof(digit_t));
            exp[k] = rand;
        }

        ec_xz_random_normalized_test(&B[0], curve);
        ec_xz_random_normalized_test(&B[1], curve);
        ec_xz_to_jac(&Bj[0], &B[0], curve);
        ec_xz_to_jac(&Bj[1], &B[1], curve);
        ec_jac_neg(&jcmp, &Bj[1]);
        ec_ADD(&Bj[2], &Bj[0], &jcmp, curve);
        ec_jac_MUL(&tmp, &Bj[1], exp, ec_n_bits(exp, NWORDS_ORDER), curve);
        ec_ADD(&jcmp, &Bj[0], &tmp, curve);
        ec_jac_to_xz(&cmp, &jcmp);
        ec_jac_to_xz(&B[0], &Bj[0]);
        ec_jac_to_xz(&B[1], &Bj[1]);
        ec_jac_to_xz(&B[2], &Bj[2]);

        if (!ec_ladder3pt(&res, &B[0], &B[1], &B[2], exp, ec_n_bits(exp, NWORDS_ORDER), curve)) {
            printf("3pt ladder rejected its inputs at iteration %d\n", (int)iter);
            return 1;
        }

        if (!ec_xz_is_equal(&res, &cmp)) {
            printf("3pt ladder test failed at iteration %d\n", (int)iter);
            return 1;
        }
    }
    printf("3pt ladder tests..................................................... PASSED\n");
    return 0;
}

int
ec_test_bary_coords(ec_curve_t *curve, unsigned int Ntest)
{
    ec_xz_point_t P, Q;
    ec_jac_point_t jP, jQ, jPpQ, jPmQ;
    ec_xz_point_t xPpQ, xPmQ, xR, xS;
    ec_bary_coords_t uvw;

    for (unsigned int i = 0; i < Ntest; i++) {
        ec_xz_random_test(&P, curve);
        ec_xz_random_test(&Q, curve);

        ec_xz_to_jac(&jP, &P, curve);
        ec_xz_to_jac(&jQ, &Q, curve);

        ec_ADD(&jPpQ, &jP, &jQ, curve); // P+Q
        ec_jac_to_xz(&xPpQ, &jPpQ);
        ec_jac_neg(&jQ, &jQ);
        ec_ADD(&jPmQ, &jP, &jQ, curve); // P-Q
        ec_jac_neg(&jQ, &jQ);
        ec_jac_to_xz(&xPmQ, &jPmQ);

        ec_jac_to_bary_coords(&uvw, &jP, &jQ, curve);

        fp2_sub(&xR.x, &uvw.u, &uvw.v);
        fp2_copy(&xR.z, &uvw.w);
        fp2_add(&xS.x, &uvw.u, &uvw.v);
        fp2_copy(&xS.z, &uvw.w);

        if (!ec_xz_is_equal(&xR, &xPpQ)) {
            printf("From curve: P+Q != (u-v:*:w)\n");
            return 1;
        }
        if (!ec_xz_is_equal(&xS, &xPmQ)) {
            printf("From curve: P-Q != (u+v:*:w)\n");
            return 1;
        }

        // Test sum zeros
        ec_xz_random_test(&P, curve);
        ec_xz_to_jac(&jP, &P, curve);

        ec_jac_point_t zero;
        ec_xz_point_t xP, xP2;
        // set zero to 0
        fp2_set_zero(&zero.x);
        fp2_set_one(&zero.y);
        fp2_set_zero(&zero.z);

        ec_jac_to_xz(&xP, &jP);
        ec_jac_to_bary_coords(&uvw, &jP, &zero, curve);
        fp2_copy(&xP2.x, &uvw.u);
        fp2_copy(&xP2.z, &uvw.w);
        if (!ec_xz_is_equal(&xP, &xP2)) {
            printf("From curve: P+0 != (xP:*:zP)\n");
            return 1;
        }
        ec_jac_to_bary_coords(&uvw, &zero, &jP, curve);
        fp2_copy(&xP2.x, &uvw.u);
        fp2_copy(&xP2.z, &uvw.w);
        if (!ec_xz_is_equal(&xP, &xP2)) {
            printf("From curve: 0+P != (xP:*:zP)\n");
            return 1;
        }
    }
    printf("Barycentric coordinates tests........................................ PASSED\n");
    return 0;
}

int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int iterations = 100;
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

    printf("--------------------------------------------------------------------------------\n\n");
    printf("Testing elliptic curve arithmetic over GF(p):\n\n");
    print_seed(seed);

    randombytes_init((unsigned char *)seed, NULL, 256);

    ec_curve_t curve;
    ec_curve_init_precomputed(&curve);

    // Test on SS curve y^2 = x^3 + 6x^2 + x
    fp2_set_small(&curve.A, 6);
    ec_compute_A24(&curve);

    res |= ec_test_jacobian(&curve, iterations);
    res |= ec_test_xDBL_xADD(&curve, iterations);
    res |= ec_test_xDBLADD(&curve, iterations);
    res |= ec_test_zero_identities(&curve, iterations);
    res |= ec_test_mul(&curve, iterations);
    res |= ec_test_basis(&curve);
    res |= ec_test_3ptladder(&curve, iterations);
    res |= ec_test_bary_coords(&curve, iterations);

    if (res) {
        printf("Tests failed!\n");
    } else {
        printf("All ec arithmetic tests passed.\n");
    }

    return res;
}
