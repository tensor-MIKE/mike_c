
#ifndef TEST_EXTRAS_H
#define TEST_EXTRAS_H

#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <constants.h>
#include <ec.h>
#include <fp.h>

#define PASSED 0
#define FAILED 1

// Generating a pseudo-random field element in [0, p-1]
void ec_fp2_random_test(fp2_t *a);

// Generating a pseudo-random element in GF(p^2)
void ec_fp2_random_test(fp2_t *a);

// Generating a random x-only projective point
void ec_xz_random_test(ec_xz_point_t *P, const ec_curve_t *curve);

// Generating a random x-only projective normalized point
void ec_xz_random_normalized_test(ec_xz_point_t *P, const ec_curve_t *curve);

#endif
