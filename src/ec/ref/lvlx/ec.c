#include <assert.h>
#include <stdio.h>
#include <ec.h>
#include <constants.h>

void
ec_curve_init_precomputed(ec_curve_t *E)
{
    // workaround since A0 is in Fp
    fp_copy(&E->A.re, &A0);
    fp_set_zero(&E->A.im);
    ec_compute_A24(E);
}

void
ec_curve_init_from_A(ec_curve_t *E, const fp2_t *A)
{
    fp2_copy(&E->A, A);
    ec_compute_A24(E);
}

void
ec_compute_ws(ec_ws_curve_t *Ews, const ec_curve_t *E)
{
    // Cost 2M + 1a
    // ao3 = A/3, a = 1-A^2/3
    fp2_t one;
    fp2_set_one(&one);
    // workaround since fp2_div3 does no exist
    fp_div3(&Ews->ao3.im, &E->A.im);
    fp_div3(&Ews->ao3.re, &E->A.re);
    fp2_mul(&Ews->a, &Ews->ao3, &E->A);
    fp2_sub(&Ews->a, &one, &Ews->a);
}

void
ec_compute_A24(ec_curve_t *E)
{
    // Cost: 2M + a
    fp2_t two;
    fp2_set_small(&two, 2);
    fp2_add(&E->A24, &E->A, &two);
    fp2_half(&E->A24, &E->A24);
    fp2_half(&E->A24, &E->A24);
}

void
ec_xz_to_jac(ec_jac_point_t *Q, const ec_xz_point_t *P, const ec_curve_t *E)
{
    fp2_t Px, Py, Pz;
    fp2_t one, z_inv, z2;
    fp2_set_one(&one);

    // Normalise the point (x:y:1)
    fp2_copy(&z_inv, &P->z);
    fp2_inv(&z_inv);
    fp2_mul(&Px, &P->x, &z_inv);
    fp2_set_one(&Pz);

    // compute y
    fp2_add(&Py, &Px, &E->A); // x+A
    fp2_mul(&Py, &Py, &Px);   // x^2+Ax
    fp2_add(&Py, &Py, &one);  // x^2 + Ax + 1
    fp2_mul(&Py, &Py, &Px);   // x^3 + Ax^2 + x
    fp2_sqrt_verify(&Py);

    // ec to jac
    fp2_sqr(&z2, &Pz);
    fp2_mul(&Q->x, &Px, &Pz);
    fp2_mul(&Q->y, &Py, &z2);
    fp2_copy(&Q->z, &Pz);

    // If Z == 0, return (0:1:0)
    uint32_t sel = fp2_is_zero(&Pz);
    fp2_select(&Q->y, &Q->y, &one, sel);
}

void
ec_jac_to_ws(ec_jac_ws_point_t *Q, const ec_jac_point_t *P, const ec_ws_curve_t *E)
{
    // X_ws = X_jac+(A*Z_jac^2)/3, Y_ws = Y_jac, Z_ws = Z_jac, T_ws = a*Z_ws^4
    // Cost of 2M + 2S + a + (precomputation ec_compute_ws(E)).

    fp2_sqr(&Q->t, &P->z);
    fp2_mul(&Q->x, &E->ao3, &Q->t);
    fp2_add(&Q->x, &Q->x, &P->x);
    fp2_sqr(&Q->t, &Q->t);
    fp2_mul(&Q->t, &Q->t, &E->a);
    fp2_copy(&Q->y, &P->y);
    fp2_copy(&Q->z, &P->z);
}

void
ec_ws_to_jac(ec_jac_point_t *Q, const ec_jac_ws_point_t *P, const ec_ws_curve_t *E)
{
    // X_jac = X_ws-(A*Z_ws^2)/3, Y_jac = Y_ws, Z_jac = Z_ws
    // Cost of 1M + 1S + a + (precomputation ec_compute_ws(E)).

    fp2_sqr(&Q->x, &P->z);
    fp2_mul(&Q->x, &E->ao3, &Q->x);
    fp2_sub(&Q->x, &P->x, &Q->x);
    fp2_copy(&Q->y, &P->y);
    fp2_copy(&Q->z, &P->z);
}

void
ec_jac_to_xz(ec_xz_point_t *Q, const ec_jac_point_t *P)
{
    // Maps infty = (0:1:0) to (1:0) and not (0:0)
    // Cost 1S

    fp2_copy(&Q->x, &P->x);
    fp2_sqr(&Q->z, &P->z);

    fp2_t one;
    uint32_t ctl = ec_jac_is_zero(P);
    fp2_set_one(&one);

    fp2_select(&Q->x, &Q->x, &one, ctl);
}

void
ec_xz_normalize(ec_xz_point_t *P)
{
    fp2_t one, zero;
    uint32_t ctl;
    fp2_set_one(&one);
    fp2_set_zero(&zero);

    ctl = fp2_is_zero(&P->z);
    fp2_inv(&P->z);
    fp2_mul(&P->x, &P->z, &P->x);

    fp2_select(&P->x, &P->x, &one, ctl);
    fp2_select(&P->z, &one, &zero, ctl);
}

uint32_t
ec_jac_is_zero(const ec_jac_point_t *P)
{
    return fp2_is_zero(&P->z);
}

uint32_t
ec_jac_is_equal(const ec_jac_point_t *P, const ec_jac_point_t *Q)
{ // Evaluate if two points in Jacobian Montgomery coordinates (X:Y:Z) are equal
  // Returns 0xFFFFFFFF (true) if P=Q, 0 (false) otherwise
    fp2_t u0, u1, t0, t1, t2, t3;

    // Check if P, Q are the points at infinity
    uint32_t l_zero = ec_jac_is_zero(P);
    uint32_t r_zero = ec_jac_is_zero(Q);

    // Check if PX * QZ^2 = QX * PZ^2 AND PY * QZ^3 = QY * PZ^3
    fp2_sqr(&u0, &Q->z);      // u0 = QZ^2
    fp2_sqr(&u1, &P->z);      // u1 = PZ^2
    fp2_mul(&t0, &P->x, &u0); // t0 = PX*QZ^2
    fp2_mul(&t1, &u1, &Q->x); // t1 = QX*PZ^2
    fp2_mul(&u0, &u0, &Q->z); // u0 = QZ^3
    fp2_mul(&u1, &u1, &P->z); // u1 = PZ^3
    fp2_mul(&t2, &P->y, &u0); // t2 = PY*QZ^3
    fp2_mul(&t3, &u1, &Q->y); // t3 = QY*PZ^3
    uint32_t lr_equal = fp2_is_equal(&t0, &t1) & fp2_is_equal(&t2, &t3);

    // Points are equal if
    // - Both are zero, or
    // - neither are zero AND PX * QZ = QX * PZ AND PY * QZ = QY * PZ
    return (l_zero & r_zero) | (~l_zero & ~r_zero & lr_equal);
}

uint32_t
ec_xz_is_equal(const ec_xz_point_t *P, const ec_xz_point_t *Q)
{ // Evaluate if two points in Jacobian Montgomery coordinates (X:Z) are equal
  // Returns 0xFFFFFFFF (true) if P=Q, 0 (false) otherwise
    fp2_t t0, t1;

    // Check if P, Q are the points at infinity
    uint32_t l_zero = ec_xz_is_zero(P);
    uint32_t r_zero = ec_xz_is_zero(Q);

    // Check if PX * QZ = QX * PZ
    fp2_mul(&t0, &P->x, &Q->z);
    fp2_mul(&t1, &P->z, &Q->x);
    uint32_t lr_equal = fp2_is_equal(&t0, &t1);

    // Points are equal if
    // - Both are zero, or
    // - neither are zero AND PX * QZ = QX * PZ
    return (l_zero & r_zero) | (~l_zero & ~r_zero & lr_equal);
}

uint32_t
ec_xz_is_zero(const ec_xz_point_t *P)
{
    return fp2_is_zero(&P->z);
}

void
ec_jac_neg(ec_jac_point_t *Q, const ec_jac_point_t *P)
{
    fp2_copy(&Q->x, &P->x);
    fp2_neg(&Q->y, &P->y);
    fp2_copy(&Q->z, &P->z);
}

void
ec_ADD(ec_jac_point_t *R, const ec_jac_point_t *P, const ec_jac_point_t *Q, const ec_curve_t *E)
{
    // Addition on a Montgomery curve, representation in Jacobian coordinates (X:Y:Z) corresponding
    // to (x,y) = (X/Z^2,Y/Z^3) This version receives the coefficient value A
    //
    // Complete routine, to handle all edge cases:
    //   if ZP == 0:            # P == inf
    //       return Q
    //   if ZQ == 0:            # Q == inf
    //       return P
    //   dy <- YQ*ZP**3 - YP*ZQ**3
    //   dx <- XQ*ZP**2 - XP*ZQ**2
    //   if dx == 0:             # x1 == x2
    //       if dy == 0:         # ... and y1 == y2: doubling case
    //           dy <- ZP*ZQ * (3*XP^2 + ZP^2 * (2*A*XP + ZP^2))
    //           dx <- 2*YP*ZP
    //       else:              # ... but y1 != y2, thus P = -Q
    //           return inf
    //   XR <- dy**2 - dx**2 * (A*ZP^2*ZQ^2 + XP*ZQ^2 + XQ*ZP^2)
    //   YR <- dy * (XP*ZQ^2 * dx^2 - XR) - YP*ZQ^3 * dx^3
    //   ZR <- dx * ZP * ZQ

    // Constant time processing:
    // - The case for P == 0 or Q == 0 is handled at the end with conditional select
    // - dy and dx are computed for both the normal and doubling cases, we switch when
    //   dx == dy == 0 for the normal case.
    // - If we have that P = -Q then dx = 0 and so ZR will be zero, giving us the point
    //   at infinity for "free".
    //
    // Cost 17M + 6S + 13a
    fp2_t t0, t1, t2, t3, u1, u2, v1, dx, dy;

    // If P is zero or Q is zero we will conditionally swap before returning.
    uint32_t ctl1 = fp2_is_zero(&P->z);
    uint32_t ctl2 = fp2_is_zero(&Q->z);

    // Precompute some values
    fp2_sqr(&t0, &P->z); // t0 = z1^2
    fp2_sqr(&t1, &Q->z); // t1 = z2^2

    // Compute dy and dx for ordinary case
    fp2_mul(&v1, &t1, &Q->z); // v1 = z2^3
    fp2_mul(&t2, &t0, &P->z); // t2 = z1^3
    fp2_mul(&v1, &v1, &P->y); // v1 = y1z2^3
    fp2_mul(&t2, &t2, &Q->y); // t2 = y2z1^3
    fp2_sub(&dy, &t2, &v1);   // dy = y2z1^3 - y1z2^3
    fp2_mul(&u2, &t0, &Q->x); // u2 = x2z1^2
    fp2_mul(&u1, &t1, &P->x); // u1 = x1z2^2
    fp2_sub(&dx, &u2, &u1);   // dx = x2z1^2 - x1z2^2

    // Compute dy and dx for doubling case
    fp2_add(&t1, &P->y, &P->y); // dx_dbl = t1 = 2y1
    fp2_add(&t2, &E->A, &E->A); // t2 = 2A
    fp2_mul(&t2, &t2, &P->x);   // t2 = 2Ax1
    fp2_add(&t2, &t2, &t0);     // t2 = 2Ax1 + z1^2
    fp2_mul(&t2, &t2, &t0);     // t2 = z1^2 * (2Ax1 + z1^2)
    fp2_sqr(&t0, &P->x);        // t0 = x1^2
    fp2_add(&t2, &t2, &t0);     // t2 = x1^2 + z1^2 * (2Ax1 + z1^2)
    fp2_add(&t2, &t2, &t0);     // t2 = 2*x1^2 + z1^2 * (2Ax1 + z1^2)
    fp2_add(&t2, &t2, &t0);     // t2 = 3*x1^2 + z1^2 * (2Ax1 + z1^2)
    fp2_mul(&t2, &t2, &Q->z);   // dy_dbl = t2 = z2 * (3*x1^2 + z1^2 * (2Ax1 + z1^2))

    // If dx is zero and dy is zero swap with double variables
    uint32_t ctl = fp2_is_zero(&dx) & fp2_is_zero(&dy);
    fp2_select(&dx, &dx, &t1, ctl);
    fp2_select(&dy, &dy, &t2, ctl);

    // Some more precomputations
    fp2_mul(&t0, &P->z, &Q->z); // t0 = z1z2
    fp2_sqr(&t1, &t0);          // t1 = (z1z2)^2
    fp2_sqr(&t2, &dx);          // t2 = dx^2
    fp2_sqr(&t3, &dy);          // t3 = dy^2

    // Compute x3 = dy**2 - dx**2 * (A*ZP^2*ZQ^2 + XP*ZQ^2 + XQ*ZP^2)
    fp2_mul(&R->x, &E->A, &t1); // x3 = A*(z1z2)^2
    fp2_add(&R->x, &R->x, &u1); // x3 = A*(z1z2)^2 + u1
    fp2_add(&R->x, &R->x, &u2); // x3 = A*(z1z2)^2 + u1 + u2
    fp2_mul(&R->x, &R->x, &t2); // x3 = dx^2 * (A*(z1z2)^2 + u1 + u2)
    fp2_sub(&R->x, &t3, &R->x); // x3 = dy^2 - dx^2 * (A*(z1z2)^2 + u1 + u2)

    // Compute y3 = dy * (XP*ZQ^2 * dx^2 - XR) - YP*ZQ^3 * dx^3
    fp2_mul(&R->y, &u1, &t2);     // y3 = u1 * dx^2
    fp2_sub(&R->y, &R->y, &R->x); // y3 = u1 * dx^2 - x3
    fp2_mul(&R->y, &R->y, &dy);   // y3 = dy * (u1 * dx^2 - x3)
    fp2_mul(&t3, &t2, &dx);       // t3 = dx^3
    fp2_mul(&t3, &t3, &v1);       // t3 = v1 * dx^3
    fp2_sub(&R->y, &R->y, &t3);   // y3 = dy * (u1 * dx^2 - x3) - v1 * dx^3

    // Compute z3 = dx * z1 * z2
    fp2_mul(&R->z, &dx, &t0);

    // Finally, we need to set R = P is Q.Z = 0 and R = Q if P.Z = 0
    ec_select_jac_point(R, R, Q, ctl1);
    ec_select_jac_point(R, R, P, ctl2);
}

void
ec_DBL(ec_jac_point_t *Q, const ec_jac_point_t *P, const ec_curve_t *E)
{
    // Cost of 6M + 6S + 14a.
    // Doubling on a Montgomery curve, representation in Jacobian coordinates (X:Y:Z) corresponding to
    // (X/Z^2,Y/Z^3) This version receives the coefficient value A
    fp2_t t0, t1, t2, t3;

    uint32_t flag = fp2_is_zero(&P->x) & fp2_is_zero(&P->z);

    fp2_sqr(&t0, &P->x); // t0 = x1^2
    fp2_add(&t1, &t0, &t0);
    fp2_add(&t0, &t0, &t1); // t0 = 3x1^2
    fp2_sqr(&t1, &P->z);    // t1 = z1^2
    fp2_mul(&t2, &P->x, &E->A);
    fp2_add(&t2, &t2, &t2); // t2 = 2Ax1
    fp2_add(&t2, &t1, &t2); // t2 = 2Ax1+z1^2
    fp2_mul(&t2, &t1, &t2); // t2 = z1^2(2Ax1+z1^2)
    fp2_add(&t2, &t0, &t2); // t2 = alpha = 3x1^2 + z1^2(2Ax1+z1^2)
    fp2_mul(&Q->z, &P->y, &P->z);
    fp2_add(&Q->z, &Q->z, &Q->z); // z2 = 2y1z1
    fp2_sqr(&t0, &Q->z);
    fp2_mul(&t0, &t0, &E->A); // t0 = 4Ay1^2z1^2
    fp2_sqr(&t1, &P->y);
    fp2_add(&t1, &t1, &t1);     // t1 = 2y1^2
    fp2_add(&t3, &P->x, &P->x); // t3 = 2x1
    fp2_mul(&t3, &t1, &t3);     // t3 = 4x1y1^2
    fp2_sqr(&Q->x, &t2);        // x2 = alpha^2
    fp2_sub(&Q->x, &Q->x, &t0); // x2 = alpha^2 - 4Ay1^2z1^2
    fp2_sub(&Q->x, &Q->x, &t3);
    fp2_sub(&Q->x, &Q->x, &t3); // x2 = alpha^2 - 4Ay1^2z1^2 - 8x1y1^2
    fp2_sub(&Q->y, &t3, &Q->x); // y2 = 4x1y1^2 - x2
    fp2_mul(&Q->y, &Q->y, &t2); // y2 = alpha(4x1y1^2 - x2)
    fp2_sqr(&t1, &t1);          // t1 = 4y1^4
    fp2_sub(&Q->y, &Q->y, &t1);
    fp2_sub(&Q->y, &Q->y, &t1); // y2 = alpha(4x1y1^2 - x2) - 8y1^4

    fp2_select(&Q->x, &Q->x, &P->x, -flag); // What kind of magic is that?
    fp2_select(&Q->z, &Q->z, &P->z, -flag);
}

void
ec_DBL_ws(ec_jac_ws_point_t *Q, const ec_jac_ws_point_t *P)
{
    // Cost of 3M + 5S + 14a.
    // Doubling on a Weierstrass curve, representation in modified Jacobian coordinates
    // (X:Y:Z:T=a*Z^4) corresponding to (X/Z^2,Y/Z^3), where a is the curve coefficient.
    // Formula from https://hyperelliptic.org/EFD/g1p/auto-shortw-modified.html

    fp2_t xx, c, cc, r, s, m;
    // XX = X^2
    fp2_sqr(&xx, &P->x);
    // A = 2*Y^2
    fp2_sqr(&c, &P->y);
    fp2_add(&c, &c, &c);
    // AA = A^2
    fp2_sqr(&cc, &c);
    // R = 2*AA
    fp2_add(&r, &cc, &cc);
    // S = (X+A)^2-XX-AA
    fp2_add(&s, &P->x, &c);
    fp2_sqr(&s, &s);
    fp2_sub(&s, &s, &xx);
    fp2_sub(&s, &s, &cc);
    // M = 3*XX+T1
    fp2_add(&m, &xx, &xx);
    fp2_add(&m, &m, &xx);
    fp2_add(&m, &m, &P->t);
    // X3 = M^2-2*S
    fp2_sqr(&Q->x, &m);
    fp2_sub(&Q->x, &Q->x, &s);
    fp2_sub(&Q->x, &Q->x, &s);
    // Z3 = 2*Y*Z
    fp2_mul(&Q->z, &P->y, &P->z);
    fp2_add(&Q->z, &Q->z, &Q->z);
    // Y3 = M*(S-X3)-R
    fp2_sub(&Q->y, &s, &Q->x);
    fp2_mul(&Q->y, &Q->y, &m);
    fp2_sub(&Q->y, &Q->y, &r);
    // T3 = 2*R*T1
    fp2_mul(&Q->t, &P->t, &r);
    fp2_add(&Q->t, &Q->t, &Q->t);
}

void
ec_jac_dbl_iter(ec_jac_point_t *Q,
                const ec_jac_point_t *P,
                const unsigned int n,
                const ec_curve_t *E,
                const int number_of_points)
{
    if (n <= 2) {
        for (int iter = 0; iter < number_of_points; iter++) {
            if (n == 0) {
                ec_copy_jac_point(&Q[iter], &P[iter]);
            } else {
                ec_DBL(&Q[iter], &P[iter], E);
                for (unsigned int i = 1; i < n; i++) {
                    ec_DBL(&Q[iter], &Q[iter], E);
                }
            }
        }
    } else {
        ec_jac_ws_point_t R;
        ec_ws_curve_t Ews;
        ec_compute_ws(&Ews, E);

        for (int iter = 0; iter < number_of_points; iter++) {
            ec_jac_to_ws(&R, &P[iter], &Ews);

            for (unsigned int i = 0; i < n; i++) {
                ec_DBL_ws(&R, &R);
            }

            ec_ws_to_jac(&Q[iter], &R, &Ews);
        }
    }
}

void
ec_jac_MUL(ec_jac_point_t *Q, const ec_jac_point_t *P, const digit_t *k, const int kbits, const ec_curve_t *curve)
{

    ec_xz_point_t P0, Q0, R0;
    ec_jac_to_xz(&P0, P);
    ec_xMUL_shifted(&Q0, &R0, &P0, k, kbits, curve, 0);
    ec_recover_y(Q, P, &Q0, &R0, curve);
}

void
ec_jac_MUL_shifted(ec_jac_point_t *Q,
                   const ec_jac_point_t *P,
                   const digit_t *k,
                   const int kbits,
                   const ec_curve_t *curve,
                   const int k_right_shift)
{

    ec_xz_point_t P0, Q0, R0;
    ec_jac_to_xz(&P0, P);
    ec_xMUL_shifted(&Q0, &R0, &P0, k, kbits, curve, k_right_shift);
    ec_recover_y(Q, P, &Q0, &R0, curve);
}

void
ec_jac_to_bary_coords(ec_bary_coords_t *bary_coords,
                      const ec_jac_point_t *P,
                      const ec_jac_point_t *Q,
                      const ec_curve_t *E)
{
    // Take P and Q in E distinct, two jac_point_t, return three components u,v and w in the base fp such
    // that the xz coordinates of P+Q are (u-v:w) and of P-Q are (u+v:w)
    // Now also works when P == 0 or Q == 0 (but not both)
    // Cost 11M + 5S + 7a

    fp2_t t0, t1, t2, t3, t4, t5, t6, xP, xQ, zP2, zQ2;
    uint32_t ctl1 = ec_jac_is_zero(P), ctl2 = ec_jac_is_zero(Q);

    fp2_copy(&xP, &P->x);
    fp2_copy(&xQ, &Q->x);

    fp2_sqr(&zP2, &P->z);               // zP2 = z1^2
    fp2_sqr(&zQ2, &Q->z);               // zQ2 = z2^2
    fp2_mul(&t2, &P->x, &zQ2);          // t2 = x1z2^2
    fp2_mul(&t3, &zP2, &Q->x);          // t3 = z1^2x2
    fp2_mul(&t4, &P->y, &Q->z);         // t4 = y1z2
    fp2_mul(&t4, &t4, &zQ2);            // t4 = y1z2^3
    fp2_mul(&t5, &P->z, &Q->y);         // t5 = z1y2
    fp2_mul(&t5, &t5, &zP2);            // t5 = z1^3y2
    fp2_mul(&t0, &zP2, &zQ2);           // t0 = (z1z2)^2
    fp2_mul(&t6, &t4, &t5);             // t6 = (z1z_2)^3y1y2
    fp2_add(&bary_coords->v, &t6, &t6); // v  = 2(z1z_2)^3y1y2
    fp2_sqr(&t4, &t4);                  // t4 = y1^2z2^6
    fp2_sqr(&t5, &t5);                  // t5 = z1^6y_2^2
    fp2_add(&t4, &t4, &t5);             // t4 = z1^6y_2^2 + y1^2z2^6
    fp2_add(&t5, &t2, &t3);             // t5 = x1z2^2 +z_1^2x2
    fp2_add(&t6, &t3, &t3);             // t6 = 2z_1^2x2
    fp2_sub(&t6, &t5, &t6);             // t6 = lambda = x1z2^2 - z_1^2x2
    fp2_sqr(&t6, &t6);                  // t6 = lambda^2 = (x1z2^2 - z_1^2x2)^2
    fp2_mul(&t1, &E->A, &t0);           // t1 = A*(z1z2)^2
    fp2_add(&t1, &t5, &t1);             // t1 = gamma =A*(z1z2)^2 + x1z2^2 +z_1^2x2
    fp2_mul(&t1, &t1, &t6);             // t1 = gamma*lambda^2
    fp2_sub(&bary_coords->u, &t4, &t1); // u  = z1^6y_2^2 + y1^2z2^6 - gamma*lambda^2
    fp2_mul(&bary_coords->w, &t6, &t0); // w  = (z1z2)^2(lambda)^2

    fp2_set_zero(&t5);
    // If P == 0, (u:v:w)=(x_Q:0:z_Q^2)
    fp2_select(&bary_coords->u, &bary_coords->u, &xQ, ctl1);
    fp2_select(&bary_coords->v, &bary_coords->v, &t5, ctl1 | ctl2);
    fp2_select(&bary_coords->w, &bary_coords->w, &zQ2, ctl1);

    // If Q == 0, (u:v:w)=(x_P:0:z_P^2)
    fp2_select(&bary_coords->u, &bary_coords->u, &xP, ctl2);
    fp2_select(&bary_coords->w, &bary_coords->w, &zP2, ctl2);
}

void
ec_xDBL(ec_xz_point_t *Q, const ec_xz_point_t *P, const ec_curve_t *E)
{
    // Doubling of a Montgomery point in projective coordinates (X:Z).
    // Input: projective Montgomery x-coordinates P = (XP:ZP), where xP=XP/ZP, and
    //        the Montgomery curve constants A24 = (A+2)/4 is (pre)computed.
    // Output: projective Montgomery x-coordinates Q <- 2*P = (XQ:ZQ) such that x(2P)=XQ/ZQ.
    // Cost 3M + 2S + 4a
    fp2_t t0, t1, t2;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sqr(&t0, &t0);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&t1, &t1);
    fp2_sub(&t2, &t0, &t1);
    fp2_mul(&Q->x, &t0, &t1);
    fp2_mul(&t0, &t2, &E->A24);
    fp2_add(&t0, &t0, &t1);
    fp2_mul(&Q->z, &t0, &t2);
}

void
ec_xADD(ec_xz_point_t *R, const ec_xz_point_t *P, const ec_xz_point_t *Q, const ec_xz_point_t *PQ)
{
    // Differential addition of Montgomery points in projective coordinates (X:Z).
    // Input: projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, and difference
    //        PQ=P-Q=(XPQ:ZPQ).
    // Output: projective Montgomery point R <- P+Q = (XR:ZR) such that x(P+Q)=XR/ZR.
    // Cost: 4M + 2S + 6a
    fp2_t t0, t1, t2, t3;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_add(&t2, &Q->x, &Q->z);
    fp2_sub(&t3, &Q->x, &Q->z);
    fp2_mul(&t0, &t0, &t3);
    fp2_mul(&t1, &t1, &t2);
    fp2_add(&t2, &t0, &t1);
    fp2_sub(&t3, &t0, &t1);
    fp2_sqr(&t2, &t2);
    fp2_sqr(&t3, &t3);
    fp2_mul(&t2, &PQ->z, &t2);
    fp2_mul(&R->z, &PQ->x, &t3);
    fp2_copy(&R->x, &t2);
}

void
ec_xDBLADD(ec_xz_point_t *R,
           ec_xz_point_t *S,
           const ec_xz_point_t *P,
           const ec_xz_point_t *Q,
           const ec_xz_point_t *PQ,
           const ec_curve_t *E)
{
    // Simultaneous doubling and differential addition.
    // Input:  projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, the difference
    //         PQ=P-Q=(XPQ:ZPQ), and the Montgomery curve constants A24 = (A+2)/4 is (pre)computed.
    // Output: projective Montgomery points R <- 2*P = (XR:ZR) such that x(2P)=XR/ZR, and S <- P+Q = (XS:ZS) such that =
    //         x(Q+P)=XS/ZS.
    // Cost: 7M + 4S + 8a (gain: it ain't much, but that's still something...)
    fp2_t t0, t1, t2;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&R->x, &t0);
    fp2_sub(&t2, &Q->x, &Q->z);
    fp2_add(&S->x, &Q->x, &Q->z);
    fp2_mul(&t0, &t0, &t2);
    fp2_sqr(&R->z, &t1);
    fp2_mul(&t1, &t1, &S->x);
    fp2_sub(&t2, &R->x, &R->z);
    fp2_mul(&R->x, &R->x, &R->z);
    fp2_mul(&S->x, &E->A24, &t2);
    fp2_sub(&S->z, &t0, &t1);
    fp2_add(&R->z, &R->z, &S->x);
    fp2_add(&S->x, &t0, &t1);
    fp2_mul(&R->z, &R->z, &t2);
    fp2_sqr(&S->z, &S->z);
    fp2_sqr(&S->x, &S->x);
    fp2_mul(&S->z, &S->z, &PQ->x);
    fp2_mul(&S->x, &S->x, &PQ->z);
}

void
ec_recover_y(ec_jac_point_t *R,
             const ec_jac_point_t *P,
             const ec_xz_point_t *Q,
             const ec_xz_point_t *PQ,
             const ec_curve_t *E)
{
    // Cost: 13M + 1S + 9a
    fp2_t v1, v2, v3, v4, v5, v6, t;

    fp2_t Px, Py, Pz;
    fp2_t Rx, Ry, Rz;
    // jac to ec on P
    fp2_t z2, one;
    fp2_sqr(&z2, &P->z);
    fp2_mul(&Pz, &z2, &P->z);
    fp2_mul(&Px, &P->x, &P->z);
    fp2_copy(&Py, &P->y);

    fp2_mul(&v1, &Px, &Q->z);  // v1 = X_P*Z_Q
    fp2_mul(&v5, &Q->x, &Pz);  // v5 = X_Q*z_P
    fp2_add(&v2, &v5, &v1);    // v2 = X_Q*z_P + X_P*Z_Q
    fp2_sub(&v3, &v5, &v1);    // v3 = X_Q*z_P - X_P*Z_Q
    fp2_sqr(&v3, &v3);         // v3 = (X_Q*z_P - X_P*Z_Q)^2
    fp2_mul(&v3, &v3, &PQ->x); // v3 = X_PQ(X_Q*z_P - X_P*Z_Q)^2
    fp2_mul(&v6, &Pz, &Q->z);  // v6 = z_P*z_Q
    fp2_mul(&v1, &E->A, &v6);  // v1 = A*z_P*_Q
    fp2_add(&v1, &v1, &v1);    // v1 = 2*A*z_P*z_Q
    fp2_add(&v2, &v1, &v2);    // v2 = X_Q*z_P + X_P*Z_Q + 2*A*z_P*z_Q
    fp2_mul(&v4, &Px, &Q->x);  // v4 = x_P*x_Q
    fp2_add(&v4, &v4, &v6);    // v4 = x_P*x_Q + z_P*z_Q
    fp2_mul(&v2, &v2, &v4);    // v2 = (x_P*x_Q + z_P*z_Q)*(X_Q*z_P + X_P*Z_Q + 2*A*z_P*z_Q)
    fp2_mul(&v1, &v1, &v6);    // v1 = 2*A*z_P^2*z_Q^2
    fp2_sub(&v2, &v2, &v1);    // v2 = (x_P*x_Q + z_P*z_Q)*(X_Q*z_P + X_P*Z_Q + 2*A*z_P*z_Q)
    // -2*A*z_P^2*z_Q^2
    fp2_mul(&v2, &v2, &PQ->z); // v2 = z_PQ*((x_P*x_Q + z_P*z_Q)*(X_Q*z_P + X_P*Z_Q + 2*A*z_P*z_Q)
    // -2*A*z_P^2*z_Q^2)
    fp2_sub(&Ry, &v2, &v3); // Y_Q = z_PQ*((x_P*x_Q + z_P*z_Q)*(X_Q*z_P + X_P*Z_Q + 2*A*z_P*z_Q)
    // -2*A*z_P^2*z_Q^2) - X_PQ(X_Q*z_P - X_P*Z_Q)^2
    fp2_add(&v1, &Py, &Py);    // v1 = 2*y_P
    fp2_mul(&v1, &v1, &v6);    // v1 = 2*y_P*z_P*z_Q
    fp2_mul(&v1, &v1, &PQ->z); // v1 = 2*y_P*z_P*z_Q*z_PQ
    fp2_mul(&Rx, &v1, &Q->x);  // v1 = 2*y_P*z_P*z_Q*z_PQ*x_Q
    fp2_mul(&Rz, &v1, &Q->z);  // v1 = 2*y_P*z_P*z_Q^2*z_PQ

    // If P == 0 or Q == 0, return R = 0
    uint32_t ctl = fp2_is_zero(&Pz) | ec_xz_is_zero(Q);
    fp2_set_one(&one);
    fp2_select(&Ry, &Ry, &one, ctl);

    // If P+Q == 0, return R = -P
    ctl = ec_xz_is_zero(PQ);
    fp2_neg(&t, &Py);
    fp2_select(&Rx, &Rx, &Px, ctl);
    fp2_select(&Ry, &Ry, &t, ctl);
    fp2_select(&Rz, &Rz, &Pz, ctl);

    // ec to jac on R
    fp2_sqr(&z2, &Rz);
    fp2_mul(&R->x, &Rx, &Rz);
    fp2_mul(&R->y, &Ry, &z2);
    fp2_copy(&R->z, &Rz);

    // If Z == 0, return (0:1:0)
    uint32_t sel = fp2_is_zero(&Rz);
    fp2_select(&R->y, &R->y, &one, sel);
}

void
ec_xMUL_shifted(ec_xz_point_t *Q,
                ec_xz_point_t *R,
                const ec_xz_point_t *P,
                const digit_t *k,
                const int kbits,
                const ec_curve_t *curve,
                const int k_right_shift)
{ // The Montgomery ladder
  // Input: projective Montgomery point P=(XP:ZP) such that xP=XP/ZP, a scalar k of bitlength kbits, and
  //        the Montgomery curve constants (A:C) (or A24 = (A+2C/4C:1) if normalized).
  // Output: projective Montgomery points Q <- k*P = (XQ:ZQ) such that x(k*P)=XQ/ZQ
    // and R <- [k+1]*P = (XR:ZR) such that x([k+1]*P)=XR/ZR
    // Uses https://eprint.iacr.org/2017/212, Algorithm 6.
    ec_xz_point_t R0, R1;
    uint32_t mask;
    unsigned int bit, prevbit = 0, swap;

    // R0 <- P, R1 <- 2*P
    ec_copy_xz_point(&R0, P);
    ec_xDBL(&R1, &R0, curve);


    // Main loop
    for (int i = kbits - 1; i >= k_right_shift; i--) {
        bit = (k[i >> LOG2RADIX] >> (i & (RADIX - 1))) & 1;
        swap = bit ^ prevbit;
        prevbit = bit;
        mask = 0 - (uint32_t)swap;

        ec_cswap_xz_points(&R0, &R1, mask);
        ec_xDBLADD(&R0, &R1, &R0, &R1, P, curve);
    }

    swap = 0 ^ prevbit;
    mask = 0 - (uint32_t)swap;

    ec_cswap_xz_points(&R0, &R1, mask);

    ec_copy_xz_point(Q, &R0);
    ec_copy_xz_point(R, &R1);
}

void
ec_xMUL(ec_xz_point_t *Q,
        ec_xz_point_t *R,
        const ec_xz_point_t *P,
        const digit_t *k,
        const int kbits,
        const ec_curve_t *curve)
{
    ec_xMUL_shifted(Q, R, P, k, kbits, curve, 0);
}

int
ec_ladder3pt(ec_xz_point_t *R,
             const ec_xz_point_t *P,
             const ec_xz_point_t *Q,
             const ec_xz_point_t *PQ,
             const digit_t *k,
             const int kbits,
             const ec_curve_t *E)
{ // The 3-point Montgomery ladder
  // Input:  projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, a scalar k of
  //         bitlength kbits, the difference PQ=P-Q=(XPQ:ZPQ), and the Montgomery curve constants A24 = (A+2C/4C:1).
  // Output: projective Montgomery point R <- P + m*Q = (XR:ZR) such that x(P + m*Q)=XR/ZR.
    // Formulas are not valid in that case
    // originally A24.z!=0
    if (fp2_is_zero(&E->A24)) {
        printf("Bad A24 in ladder3pt\n");
        return 0;
    }
    if (fp2_is_zero(&PQ->x) | fp2_is_zero(&PQ->z)) {
        printf("bad diff in ladder3pt\n");
        return 0;
    }

    ec_xz_point_t X0, X1, X2;
    ec_copy_xz_point(&X0, Q);
    ec_copy_xz_point(&X1, P);
    ec_copy_xz_point(&X2, PQ);
    uint32_t mask;
    unsigned int bit = 0;

    for (int i = 0; i < kbits + 1; i++) {
        bit = ((k[i >> LOG2RADIX] >> (i & (RADIX - 1))) & 1);
        mask = 0 - !(uint32_t)bit;
        ec_cswap_xz_points(&X1, &X2, mask);
        ec_xDBLADD(&X0, &X1, &X0, &X1, &X2, E);
        ec_cswap_xz_points(&X1, &X2, mask);
    }

    ec_copy_xz_point(R, &X1);
    return 1;
}
