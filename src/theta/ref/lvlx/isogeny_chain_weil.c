
#include <theta.h>
#include <assert.h>

#ifndef NDEBUG
// Compute bitwise scalar product of i and j mod 2
static uint8_t
scalprod(uint8_t i, uint8_t j)
{
    uint8_t ret = 0;
    ret ^= (i & 1) & (j & 1);
    ret ^= ((i >> 1) & 1) & ((j >> 1) & 1);
    ret ^= ((i >> 2) & 1) & ((j >> 2) & 1);
    ret ^= ((i >> 3) & 1) & ((j >> 3) & 1);
    return ret;
}

// Debug function used to check at that at each step of the chain, everything is ok.
static uint32_t
debug_test_codomain(theta_point_dim4_t *T3_8, theta_point_dim4_t *T4_8, theta_struct_weil_t *codomain, int step)
{
    theta_point_dim4_t T3_4, T4_4;
    uint32_t ret = -1;
    fp_t tmp_1, tmp_2, tmp_3, tmp_4, tmp_5;

    // evaluate through the 2-isogeny
    isogeny_weil_eval(&T3_4, T3_8, codomain, step != 0);
    isogeny_weil_eval(&T4_4, T4_8, codomain, step != 0);

    theta_struct_weil_precomp_t precomp;
    theta_struct_weil_arith_precomp(&precomp, codomain);

    theta_hadamard_dim4(&T3_4, &T3_4);
    theta_hadamard_dim4(&T4_4, &T4_4);

    // Test 4 torsion
    for (uint8_t i = 0; i < 16; i++) {
        if (scalprod(i, 4) == 1) {
            ret &= fp_is_zero(&T3_4[i]);
        }
    }
    for (uint8_t i = 0; i < 16; i++) {
        if (scalprod(i, 8) == 1) {
            ret &= fp_is_zero(&T4_4[i]);
        }
    }
    theta_hadamard_dim4(&T3_4, &T3_4);
    theta_hadamard_dim4(&T4_4, &T4_4);

    // Test 2 torsion is indeed some 2 torsion
    theta_point_dim4_t null_point_3, null_point_4;

    theta_DBL_weil(&T3_4, &T3_4, &precomp);
    theta_DBL_weil(&null_point_3, &T3_4, &precomp);

    theta_hadamard_dim4(&T3_4, &T3_4);
    theta_hadamard_dim4(&null_point_3, &null_point_3);

    fp_copy(&tmp_1, &T3_4[0]);
    fp_inv(&tmp_1);
    fp_copy(&tmp_3, &null_point_3[0]);
    fp_inv(&tmp_3);
    for (uint8_t i = 0; i < 16; i++) {

        fp_mul(&tmp_2, &tmp_1, &T3_4[i]);
        fp_mul(&tmp_4, &tmp_3, &null_point_3[i]);

        if (scalprod(i, 4) == 0) {
            fp_sub(&tmp_5, &tmp_2, &tmp_4);
        } else {
            fp_add(&tmp_5, &tmp_2, &tmp_4);
        }

        ret &= fp_is_zero(&tmp_5);
    }

    theta_hadamard_dim4(&T3_4, &T3_4);
    theta_hadamard_dim4(&null_point_3, &null_point_3);

    theta_DBL_weil(&T4_4, &T4_4, &precomp);
    theta_DBL_weil(&null_point_4, &T4_4, &precomp);

    theta_hadamard_dim4(&T4_4, &T4_4);
    theta_hadamard_dim4(&null_point_4, &null_point_4);

    fp_copy(&tmp_1, &T4_4[0]);
    fp_inv(&tmp_1);
    fp_copy(&tmp_3, &null_point_4[0]);
    fp_inv(&tmp_3);
    for (uint8_t i = 0; i < 16; i++) {

        fp_mul(&tmp_2, &tmp_1, &T4_4[i]);
        fp_mul(&tmp_4, &tmp_3, &null_point_4[i]);

        if (scalprod(i, 8) == 0) {
            fp_sub(&tmp_5, &tmp_2, &tmp_4);
        } else {
            fp_add(&tmp_5, &tmp_2, &tmp_4);
        }

        ret &= fp_is_zero(&tmp_5);
    }

    theta_hadamard_dim4(&T4_4, &T4_4);
    theta_hadamard_dim4(&null_point_4, &null_point_4);

    // Test we have the same theta null point
    fp_copy(&tmp_1, &null_point_3[0]);
    fp_inv(&tmp_1);
    fp_copy(&tmp_3, &null_point_4[0]);
    fp_inv(&tmp_3);
    for (uint8_t i = 0; i < 16; i++) {
        fp_mul(&tmp_2, &tmp_1, &null_point_3[i]);
        fp_mul(&tmp_4, &tmp_3, &null_point_4[i]);

        ret &= fp_is_equal(&tmp_2, &tmp_4);
    }

    // Test it is the same theta null point as the one of the codomain !
    theta_point_weil_t dual_null_point;
    theta_invert_weil(&dual_null_point, &codomain->inv_dual_null_point);
    theta_weil_to_dim4_point(&null_point_4, &dual_null_point);
    fp_copy(&tmp_1, &null_point_3[0]);
    fp_inv(&tmp_1);
    fp_copy(&tmp_3, &null_point_4[0]);
    fp_inv(&tmp_3);
    for (uint8_t i = 0; i < 16; i++) {
        fp_mul(&tmp_2, &tmp_1, &null_point_3[i]);
        fp_mul(&tmp_4, &tmp_3, &null_point_4[i]);

        ret &= fp_is_equal(&tmp_2, &tmp_4);
    }

    return ret;
}
#endif

void
isogeny_chain_weil(theta_struct_weil_t *codomain,
                   theta_struct_weil_precomp_t *codomain_precomp,
                   const theta_point_dim4_t *T1,
                   const theta_point_dim4_t *T2,
                   const unsigned int len)
{
    // co_domain is the input domain and output codomain
    // Space to store strategy points
    int space = 1;
    int step = 0;
    int dbl_loop = 0;
    int j; 
    for (unsigned int i = 1; i < len; i *= 2) {
        ++space;
    }
    // To store point orders
    uint16_t orders[space];
    orders[0] = len;
    int current = 0;   // Index of where we are in the strategy
    uint16_t num_dbls; // Number of doublings
    // To store multiples of T1 and T2
    theta_point_dim4_t ker[space][2];
    theta_copy_theta_dim4(&ker[0][0], T1);
    theta_copy_theta_dim4(&ker[0][1], T2);

    while (current >= 0 && orders[current]) {
        assert(current < space);
        while (orders[current] != 1) {
            if ((dbl_loop == 0) & (step != 0)) // compute the needed values for doubling (except the first times)
                theta_struct_weil_arith_precomp(codomain_precomp, codomain);
            dbl_loop++;
            assert(orders[current] >= 2);
            ++current;
            assert(current < space);
            num_dbls = orders[current - 1] / 2;
            assert(num_dbls && num_dbls < orders[current - 1]);
            theta_DBL_iter_weil(&ker[current][0], &ker[current - 1][0], codomain_precomp, num_dbls);
            theta_DBL_iter_weil(&ker[current][1], &ker[current - 1][1], codomain_precomp, num_dbls);
            orders[current] = orders[current - 1] - num_dbls;
        }
        dbl_loop = 0;

        isogeny_weil_compute_codomain(codomain, ker[current], step != 0);

#ifndef NDEBUG
        // Test that the codomain is correctly computed.
        uint32_t correct = debug_test_codomain(&ker[current][0], &ker[current][1], codomain, step);
        if (!correct) {
            printf("problem during 4D chain\n");
            assert(false);
        }
#endif

        // pushing the kernel
        assert(orders[current] == 1);
        for (j = 0; j < current; ++j) {
            isogeny_weil_eval(&ker[j][0], &ker[j][0], codomain, step != 0);
            isogeny_weil_eval(&ker[j][1], &ker[j][1], codomain, step != 0);
            assert(orders[j]);
            --orders[j];
        }

        --current;
        step++;
    }
    assert(current == -1);
}
