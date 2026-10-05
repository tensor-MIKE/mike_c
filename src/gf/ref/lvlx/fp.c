#include <fp.h>
#include <encoded_sizes.h>
#include <fp_constants.h>
#include <ec_params.h>
#include <constants.h>
#include <stdio.h>

/*
 * If ctl == 0x00000000, then *d is set to a0
 * If ctl == 0xFFFFFFFF, then *d is set to a1
 * ctl MUST be either 0x00000000 or 0xFFFFFFFF.
 */
void
fp_select(fp_t *d, const fp_t *a0, const fp_t *a1, uint32_t ctl)
{
    digit_t cw = (int32_t)ctl;
    for (unsigned int i = 0; i < NWORDS_FIELD; i++) {
        (*d)[i] = (*a0)[i] ^ (cw & ((*a0)[i] ^ (*a1)[i]));
    }
}

uint32_t 
fp_is_one(const fp_t *a){
    return fp_is_equal(a, &ONE);
}

void
fp_add_one(fp_t *x, const fp_t *y){
    fp_add(x,y,&ONE);
}

uint32_t
fp_sqrt_verify(fp_t *a)
{
    fp_t t0, t1;

    fp_copy(&t0, a);
    fp_sqrt(a);
    fp_sqr(&t1, a);

    return (fp_is_equal(&t0, &t1));
}

void
fp_batched_inv(fp_t *x, const int len)
{
    fp_t t1[len], t2[len];
    fp_t inverse;

    // x = x0,...,xn
    // t1 = x0, x0*x1, ... ,x0 * x1 * ... * xn
    fp_copy(&t1[0], &x[0]);
    for (int i = 1; i < len; i++) {
        fp_mul(&t1[i], &t1[i - 1], &x[i]);
    }

    // inverse = 1/ (x0 * x1 * ... * xn)
    fp_copy(&inverse, &t1[len - 1]);
    fp_inv(&inverse);

    fp_copy(&t2[0], &inverse);
    // t2 = 1/ (x0 * x1 * ... * xn), 1/ (x0 * x1 * ... * x(n-1)) , ... , 1/xO
    for (int i = 1; i < len; i++) {
        fp_mul(&t2[i], &t2[i - 1], &x[len - i]);
    }

    fp_copy(&x[0], &t2[len - 1]);

    for (int i = 1; i < len; i++) {
        fp_mul(&x[i], &t1[i - 1], &t2[len - i - 1]);
    }
}

void
fp_proj_batched_inv(fp_t *x, const int len)
{
    // Cost: (3*len - 4)*M
    fp_t t1[len-1], t2[len];

    // x = x0,...,xn
    // t1 = x0, x0*x1, ... ,x0 * x1 * ... * x(n-1)
    fp_copy(&t1[0], &x[0]);
    for (int i = 1; i < len-1; i++) {
        fp_mul(&t1[i], &t1[i - 1], &x[i]);
    }

    // coeff = x0 * x1 * ... * xn
    //fp_copy(&coeff, &t1[len - 1]);

    fp_set_one(&t2[0]);
    // t2 = 1, xn , x(n-1) * xn, ... , x1 * ... * xn
    for (int i = 1; i < len; i++) {
        fp_mul(&t2[i], &t2[i - 1], &x[len - i]);
    }

    fp_copy(&x[0], &t2[len - 1]);

    for (int i = 1; i < len; i++) {
        fp_mul(&x[i], &t1[i - 1], &t2[len - i - 1]);
    }
}

void
fp_proj_batched_inv_with_coeff(fp_t *x, fp_t *coeff, const int len)
{
    fp_t t1[len], t2[len];

    // x = x0,...,xn
    // t1 = x0, x0*x1, ... ,x0 * x1 * ... * xn
    fp_copy(&t1[0], &x[0]);
    for (int i = 1; i < len; i++) {
        fp_mul(&t1[i], &t1[i - 1], &x[i]);
    }

    // coeff = x0 * x1 * ... * xn
    fp_copy(coeff, &t1[len - 1]);

    fp_set_one(&t2[0]);
    // t2 = 1, xn , x(n-1) * xn, ... , x1 * ... * xn
    for (int i = 1; i < len; i++) {
        fp_mul(&t2[i], &t2[i - 1], &x[len - i]);
    }

    fp_copy(&x[0], &t2[len - 1]);

    for (int i = 1; i < len; i++) {
        fp_mul(&x[i], &t1[i - 1], &t2[len - i - 1]);
    }
}

void
fp_print(const char *name, const fp_t *a)
{
    printf("%s0x", name);

    uint8_t buf[FP_ENCODED_BYTES];
    fp_encode(&buf, a); // Encoding ensures canonical rep
    for (int i = 0; i < FP_ENCODED_BYTES; i++) {
        printf("%02x", buf[FP_ENCODED_BYTES - i - 1]);
    }
    printf("\n");
}