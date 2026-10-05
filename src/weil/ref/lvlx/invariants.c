#include <weil.h>

static inline void
relative_invariants(mike_rel_invariants_t *rel_inv, const theta_point_weil_t *theta)
{
    // Using relations F2 = 0 and F6 = 0 from Lemma H.16 and relations from Lemma H.13,
    // optimizes the computation of I4, ..., I9, requiring a0, a3, a5, a6, a9, a15 only.
    // Actually comptes 2*I1/I0, ..., 2*I9/I0.
    // Total cost: 1I + 10M + 4S + 28a

    fp_t t0, t1, t2, t3, u0, u1, u2, u3, I0;

    // Recall: weil_to_theta[10] = {0, 1, 2, 3, 5, 6, 7, 9, 11 ,15};
    fp_add(&t0, &(*theta)[0], &(*theta)[5]); // t0 = a0 + a6
    fp_sub(&t1, &(*theta)[0], &(*theta)[5]); // t1 = a0 - a6
    fp_add(&t2, &(*theta)[7], &(*theta)[9]); // t2 = a9 + a15
    fp_sub(&t3, &(*theta)[7], &(*theta)[9]); // t3 = a9 - a15

    fp_add(&u0, &t0, &t2); // u0 = t0 + t2 = a0 + a6 + a9 + a15
    fp_sub(&u2, &t0, &t2); // u2 = t0 - t2 = a0 + a6 - a9 - a15
    fp_add(&u1, &t1, &t3); // u1 = t1 + t3 = a0 - a6 + a9 - a15
    fp_sub(&u3, &t1, &t3); // u3 = t1 - t3 = a0 - a6 - a9 + a15

    fp_sqr(&I0, &u0);            // I0 = u0^2 = (a0 + a6 + a9 + a15)^2
    fp_sqr(&(*rel_inv)[0], &u1); // I1 = u1^2 = (a0 - a6 + a9 - a15)^2
    fp_sqr(&(*rel_inv)[1], &u2); // I2 = u2^2 = (a0 + a6 - a9 - a15)^2
    fp_sqr(&(*rel_inv)[2], &u3); // I3 = u3^2 = (a0 - a6 - a9 + a15)^2

    // u0 = 16*a3*a5
    fp_mul(&u0, &(*theta)[3], &(*theta)[4]);
    fp_add(&u0, &u0, &u0);
    fp_add(&u0, &u0, &u0);
    fp_add(&u0, &u0, &u0);
    fp_add(&u0, &u0, &u0);

    fp_add(&t0, &I0, &(*rel_inv)[0]);            // t0 = I0 + I1
    fp_sub(&t1, &I0, &(*rel_inv)[0]);            // t1 = I0 - I1
    fp_add(&t2, &(*rel_inv)[1], &(*rel_inv)[2]); // t2 = I2 + I3
    fp_sub(&t3, &(*rel_inv)[1], &(*rel_inv)[2]); // t3 = I2 - I3

    fp_sub(&u2, &t0, &t2); // u2 = t0 - t2 =  I0 + I1 - I2 - I3
    fp_add(&u1, &t1, &t3); // u1 = t1 + t3 =  I0 - I1 + I2 - I3
    fp_sub(&u3, &t1, &t3); // u3 = t1 - t3 =  I0 - I1 - I2 + I3

    fp_add(&(*rel_inv)[3], &u1, &u0); // u1 + u0 = I0 - I1 + I2 - I3 + 16*a3*a5 = 2*I4
    fp_sub(&(*rel_inv)[4], &u1, &u0); // u1 - u0 = I0 - I1 + I2 - I3 - 16*a3*a5 = 2*I5
    fp_add(&(*rel_inv)[5], &u2, &u0); // u2 + u0 = I0 + I1 - I2 - I3 + 16*a3*a5 = 2*I6
    fp_sub(&(*rel_inv)[6], &u2, &u0); // u2 - u0 = I0 + I1 - I2 - I3 - 16*a3*a5 = 2*I7
    fp_add(&(*rel_inv)[7], &u3, &u0); // u3 + u0 = I0 - I1 - I2 + I3 + 16*a3*a5 = 2*I8
    fp_sub(&(*rel_inv)[8], &u3, &u0); // u3 + u0 = I0 - I1 - I2 + I3 - 16*a3*a5 = 2*I9

    fp_inv(&I0);

    for (int i = 0; i < 9; i++) {
        fp_mul(&(*rel_inv)[i], &(*rel_inv)[i], &I0);
    }

    // We have to rescale I1, I2, I3 since the following invariants are scaled to their double
    for (int i = 0; i < 3; i++) {
        fp_add(&(*rel_inv)[i], &(*rel_inv)[i], &(*rel_inv)[i]);
    }
}

void
mike_absolute_invariants(mike_abs_invariants_t *abs_inv, const theta_point_weil_t *theta)
{
    // Total cost: 1I + 28M + 22S + 64a

    mike_rel_invariants_t rel_inv;
    fp_t t[9], u[9];

    // Initialize absolute invariants to 0
    for (int i = 0; i < 4; i++) {
        fp_set_zero(&(*abs_inv)[i]);
    }

    // Compte I1/I0, ..., I9/I0
    relative_invariants(&rel_inv, theta);

    for (int i = 0; i < 9; i++) {
        fp_sqr(&t[i], &rel_inv[i]); // ti = Ii^2
    }
    for (int i = 0; i < 9; i++) {
        fp_add(&(*abs_inv)[0], &(*abs_inv)[0], &t[i]); // J2 = sum(ti)
    }

    for (int i = 0; i < 9; i++) {
        fp_mul(&u[i], &rel_inv[i], &t[i]); // ui = Ii*ti = Ii^3
    }
    for (int i = 0; i < 9; i++) {
        fp_add(&(*abs_inv)[1], &(*abs_inv)[1], &u[i]); // J3 = sum(ui)
    }

    for (int i = 0; i < 9; i++) {
        fp_sqr(&u[i], &t[i]); // ui = ti^2 = Ii^4
    }
    for (int i = 0; i < 9; i++) {
        fp_add(&(*abs_inv)[2], &(*abs_inv)[2], &u[i]); // J4 = sum(ui)
    }

    for (int i = 0; i < 9; i++) {
        fp_mul(&u[i], &rel_inv[i], &u[i]); // ui = Ii*ui = Ii^5
    }
    for (int i = 0; i < 9; i++) {
        fp_add(&(*abs_inv)[3], &(*abs_inv)[3], &u[i]); // J5 = sum(ui)
    }
}