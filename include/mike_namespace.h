
#ifndef MIKE_NAMESPACE_H
#define MIKE_NAMESPACE_H

//#define DISABLE_NAMESPACING

#if defined(_WIN32)
#define MIKE_API __declspec(dllexport)
#else
#define MIKE_API __attribute__((visibility("default")))
#endif

#define PARAM_JOIN3_(a, b, c) mike_##a##_##b##_##c
#define PARAM_JOIN3(a, b, c) PARAM_JOIN3_(a, b, c)
#define PARAM_NAME3(end, s) PARAM_JOIN3(MIKE_VARIANT, end, s)

#define PARAM_JOIN2_(a, b) mike_##a##_##b
#define PARAM_JOIN2(a, b) PARAM_JOIN2_(a, b)
#define PARAM_NAME2(end, s) PARAM_JOIN2(end, s)

#ifndef DISABLE_NAMESPACING
#define MIKE_NAMESPACE_GENERIC(s) PARAM_NAME2(gen, s)
#else
#define MIKE_NAMESPACE_GENERIC(s) s
#endif

#if defined(MIKE_VARIANT) && !defined(DISABLE_NAMESPACING)
#if defined(MIKE_BUILD_TYPE_REF)
#define MIKE_NAMESPACE(s) PARAM_NAME3(ref, s)
#elif defined(MIKE_BUILD_TYPE_OPT)
#define MIKE_NAMESPACE(s) PARAM_NAME3(opt, s)
#elif defined(MIKE_BUILD_TYPE_BROADWELL)
#define MIKE_NAMESPACE(s) PARAM_NAME3(broadwell, s)
#elif defined(MIKE_BUILD_TYPE_ARM64)
#define MIKE_NAMESPACE(s) PARAM_NAME3(arm64, s)
#elif defined(MIKE_BUILD_TYPE_ARM64CRYPTO)
#define MIKE_NAMESPACE(s) PARAM_NAME3(arm64crypto, s)
#else
#error "Build type not known"
#endif

#else
#define MIKE_NAMESPACE(s) s
#endif

// Namespacing symbols exported from bases.c:
#undef gluing_basis_compute
#undef mike_compute_basis

#define gluing_basis_compute                                            MIKE_NAMESPACE(gluing_basis_compute)
#define mike_compute_basis                                              MIKE_NAMESPACE(mike_compute_basis)

// Namespacing symbols exported from basis.c:
#undef ec_basis_2f
#undef ec_jac_basis_2f
#undef ec_n_bits
#undef ec_x_is_on_curve

#define ec_basis_2f                                                     MIKE_NAMESPACE(ec_basis_2f)
#define ec_jac_basis_2f                                                 MIKE_NAMESPACE(ec_jac_basis_2f)
#define ec_n_bits                                                       MIKE_NAMESPACE(ec_n_bits)
#define ec_x_is_on_curve                                                MIKE_NAMESPACE(ec_x_is_on_curve)

// Namespacing symbols exported from ec.c:
#undef ec_ADD
#undef ec_DBL
#undef ec_DBL_ws
#undef ec_compute_A24
#undef ec_compute_ws
#undef ec_curve_init_from_A
#undef ec_curve_init_precomputed
#undef ec_jac_MUL
#undef ec_jac_MUL_shifted
#undef ec_jac_dbl_iter
#undef ec_jac_is_equal
#undef ec_jac_is_zero
#undef ec_jac_neg
#undef ec_jac_to_bary_coords
#undef ec_jac_to_ws
#undef ec_jac_to_xz
#undef ec_ladder3pt
#undef ec_recover_y
#undef ec_ws_to_jac
#undef ec_xADD
#undef ec_xDBL
#undef ec_xDBLADD
#undef ec_xMUL
#undef ec_xMUL_shifted
#undef ec_xz_is_equal
#undef ec_xz_is_zero
#undef ec_xz_normalize
#undef ec_xz_to_jac

#define ec_ADD                                                          MIKE_NAMESPACE(ec_ADD)
#define ec_DBL                                                          MIKE_NAMESPACE(ec_DBL)
#define ec_DBL_ws                                                       MIKE_NAMESPACE(ec_DBL_ws)
#define ec_compute_A24                                                  MIKE_NAMESPACE(ec_compute_A24)
#define ec_compute_ws                                                   MIKE_NAMESPACE(ec_compute_ws)
#define ec_curve_init_from_A                                            MIKE_NAMESPACE(ec_curve_init_from_A)
#define ec_curve_init_precomputed                                       MIKE_NAMESPACE(ec_curve_init_precomputed)
#define ec_jac_MUL                                                      MIKE_NAMESPACE(ec_jac_MUL)
#define ec_jac_MUL_shifted                                              MIKE_NAMESPACE(ec_jac_MUL_shifted)
#define ec_jac_dbl_iter                                                 MIKE_NAMESPACE(ec_jac_dbl_iter)
#define ec_jac_is_equal                                                 MIKE_NAMESPACE(ec_jac_is_equal)
#define ec_jac_is_zero                                                  MIKE_NAMESPACE(ec_jac_is_zero)
#define ec_jac_neg                                                      MIKE_NAMESPACE(ec_jac_neg)
#define ec_jac_to_bary_coords                                           MIKE_NAMESPACE(ec_jac_to_bary_coords)
#define ec_jac_to_ws                                                    MIKE_NAMESPACE(ec_jac_to_ws)
#define ec_jac_to_xz                                                    MIKE_NAMESPACE(ec_jac_to_xz)
#define ec_ladder3pt                                                    MIKE_NAMESPACE(ec_ladder3pt)
#define ec_recover_y                                                    MIKE_NAMESPACE(ec_recover_y)
#define ec_ws_to_jac                                                    MIKE_NAMESPACE(ec_ws_to_jac)
#define ec_xADD                                                         MIKE_NAMESPACE(ec_xADD)
#define ec_xDBL                                                         MIKE_NAMESPACE(ec_xDBL)
#define ec_xDBLADD                                                      MIKE_NAMESPACE(ec_xDBLADD)
#define ec_xMUL                                                         MIKE_NAMESPACE(ec_xMUL)
#define ec_xMUL_shifted                                                 MIKE_NAMESPACE(ec_xMUL_shifted)
#define ec_xz_is_equal                                                  MIKE_NAMESPACE(ec_xz_is_equal)
#define ec_xz_is_zero                                                   MIKE_NAMESPACE(ec_xz_is_zero)
#define ec_xz_normalize                                                 MIKE_NAMESPACE(ec_xz_normalize)
#define ec_xz_to_jac                                                    MIKE_NAMESPACE(ec_xz_to_jac)

// Namespacing symbols exported from encode.c:
#undef hash_shared_secret
#undef public_key_from_bytes
#undef public_key_to_bytes
#undef secret_key_from_bytes
#undef secret_key_to_bytes

#define hash_shared_secret                                              MIKE_NAMESPACE(hash_shared_secret)
#define public_key_from_bytes                                           MIKE_NAMESPACE(public_key_from_bytes)
#define public_key_to_bytes                                             MIKE_NAMESPACE(public_key_to_bytes)
#define secret_key_from_bytes                                           MIKE_NAMESPACE(secret_key_from_bytes)
#define secret_key_to_bytes                                             MIKE_NAMESPACE(secret_key_to_bytes)

// Namespacing symbols exported from exchange.c:
#undef protocols_exchange

#define protocols_exchange                                              MIKE_NAMESPACE(protocols_exchange)

// Namespacing symbols exported from fp.c:
#undef fp_add_one
#undef fp_batched_inv
#undef fp_is_one
#undef fp_print
#undef fp_proj_batched_inv
#undef fp_proj_batched_inv_with_coeff
#undef fp_select
#undef fp_sqrt_verify

#define fp_add_one                                                      MIKE_NAMESPACE(fp_add_one)
#define fp_batched_inv                                                  MIKE_NAMESPACE(fp_batched_inv)
#define fp_is_one                                                       MIKE_NAMESPACE(fp_is_one)
#define fp_print                                                        MIKE_NAMESPACE(fp_print)
#define fp_proj_batched_inv                                             MIKE_NAMESPACE(fp_proj_batched_inv)
#define fp_proj_batched_inv_with_coeff                                  MIKE_NAMESPACE(fp_proj_batched_inv_with_coeff)
#define fp_select                                                       MIKE_NAMESPACE(fp_select)
#define fp_sqrt_verify                                                  MIKE_NAMESPACE(fp_sqrt_verify)

// Namespacing symbols exported from fp2.c:
#undef fp2_add
#undef fp2_add_one
#undef fp2_batched_inv
#undef fp2_copy
#undef fp2_cswap
#undef fp2_decode
#undef fp2_encode
#undef fp2_half
#undef fp2_inv
#undef fp2_is_equal
#undef fp2_is_one
#undef fp2_is_square
#undef fp2_is_zero
#undef fp2_mul
#undef fp2_mul_small
#undef fp2_neg
#undef fp2_pow_vartime
#undef fp2_print
#undef fp2_select
#undef fp2_set_one
#undef fp2_set_small
#undef fp2_set_zero
#undef fp2_sqr
#undef fp2_sqrt
#undef fp2_sqrt_verify
#undef fp2_sub

#define fp2_add                                                         MIKE_NAMESPACE(fp2_add)
#define fp2_add_one                                                     MIKE_NAMESPACE(fp2_add_one)
#define fp2_batched_inv                                                 MIKE_NAMESPACE(fp2_batched_inv)
#define fp2_copy                                                        MIKE_NAMESPACE(fp2_copy)
#define fp2_cswap                                                       MIKE_NAMESPACE(fp2_cswap)
#define fp2_decode                                                      MIKE_NAMESPACE(fp2_decode)
#define fp2_encode                                                      MIKE_NAMESPACE(fp2_encode)
#define fp2_half                                                        MIKE_NAMESPACE(fp2_half)
#define fp2_inv                                                         MIKE_NAMESPACE(fp2_inv)
#define fp2_is_equal                                                    MIKE_NAMESPACE(fp2_is_equal)
#define fp2_is_one                                                      MIKE_NAMESPACE(fp2_is_one)
#define fp2_is_square                                                   MIKE_NAMESPACE(fp2_is_square)
#define fp2_is_zero                                                     MIKE_NAMESPACE(fp2_is_zero)
#define fp2_mul                                                         MIKE_NAMESPACE(fp2_mul)
#define fp2_mul_small                                                   MIKE_NAMESPACE(fp2_mul_small)
#define fp2_neg                                                         MIKE_NAMESPACE(fp2_neg)
#define fp2_pow_vartime                                                 MIKE_NAMESPACE(fp2_pow_vartime)
#define fp2_print                                                       MIKE_NAMESPACE(fp2_print)
#define fp2_select                                                      MIKE_NAMESPACE(fp2_select)
#define fp2_set_one                                                     MIKE_NAMESPACE(fp2_set_one)
#define fp2_set_small                                                   MIKE_NAMESPACE(fp2_set_small)
#define fp2_set_zero                                                    MIKE_NAMESPACE(fp2_set_zero)
#define fp2_sqr                                                         MIKE_NAMESPACE(fp2_sqr)
#define fp2_sqrt                                                        MIKE_NAMESPACE(fp2_sqrt)
#define fp2_sqrt_verify                                                 MIKE_NAMESPACE(fp2_sqrt_verify)
#define fp2_sub                                                         MIKE_NAMESPACE(fp2_sub)

// Namespacing symbols exported from fp_p308_633.c, fp_p474_593.c, fp_p628_317.c:
#undef fp_add
#undef fp_copy
#undef fp_cswap
#undef fp_decode
#undef fp_decode_reduce
#undef fp_div3
#undef fp_encode
#undef fp_exp3div4
#undef fp_half
#undef fp_inv
#undef fp_is_equal
#undef fp_is_square
#undef fp_is_zero
#undef fp_mul
#undef fp_mul_small
#undef fp_neg
#undef fp_set_one
#undef fp_set_small
#undef fp_set_zero
#undef fp_sqr
#undef fp_sqrt
#undef fp_sub

#define fp_add                                                          MIKE_NAMESPACE(fp_add)
#define fp_copy                                                         MIKE_NAMESPACE(fp_copy)
#define fp_cswap                                                        MIKE_NAMESPACE(fp_cswap)
#define fp_decode                                                       MIKE_NAMESPACE(fp_decode)
#define fp_decode_reduce                                                MIKE_NAMESPACE(fp_decode_reduce)
#define fp_div3                                                         MIKE_NAMESPACE(fp_div3)
#define fp_encode                                                       MIKE_NAMESPACE(fp_encode)
#define fp_exp3div4                                                     MIKE_NAMESPACE(fp_exp3div4)
#define fp_half                                                         MIKE_NAMESPACE(fp_half)
#define fp_inv                                                          MIKE_NAMESPACE(fp_inv)
#define fp_is_equal                                                     MIKE_NAMESPACE(fp_is_equal)
#define fp_is_square                                                    MIKE_NAMESPACE(fp_is_square)
#define fp_is_zero                                                      MIKE_NAMESPACE(fp_is_zero)
#define fp_mul                                                          MIKE_NAMESPACE(fp_mul)
#define fp_mul_small                                                    MIKE_NAMESPACE(fp_mul_small)
#define fp_neg                                                          MIKE_NAMESPACE(fp_neg)
#define fp_set_one                                                      MIKE_NAMESPACE(fp_set_one)
#define fp_set_small                                                    MIKE_NAMESPACE(fp_set_small)
#define fp_set_zero                                                     MIKE_NAMESPACE(fp_set_zero)
#define fp_sqr                                                          MIKE_NAMESPACE(fp_sqr)
#define fp_sqrt                                                         MIKE_NAMESPACE(fp_sqrt)
#define fp_sub                                                          MIKE_NAMESPACE(fp_sub)

// Namespacing symbols exported from gluing.c:
#undef gluing_change_theta_structure_dim2_to_U_compatible_theta_struct
#undef gluing_change_theta_structure_dim2_to_V_compatible_theta_struct
#undef gluing_compute_weil_gluing
#undef gluing_diag_isogeny_compute
#undef gluing_diag_isogeny_eval
#undef gluing_evaluate_weil_gluing
#undef gluing_from_couple_dim2_to_dim4_compatible_with_isogeny
#undef gluing_mike_compute
#undef gluing_mike_eval
#undef gluing_special_inv
#undef scholten_eval
#undef scholten_gluing_eval
#undef scholten_jacobian_to_theta
#undef scholten_precompute
#undef scholten_variable_compute
#undef schoten_gluing_compute_and_verify

#define gluing_change_theta_structure_dim2_to_U_compatible_theta_struct MIKE_NAMESPACE(gluing_change_theta_structure_dim2_to_U_compatible_theta_struct)
#define gluing_change_theta_structure_dim2_to_V_compatible_theta_struct MIKE_NAMESPACE(gluing_change_theta_structure_dim2_to_V_compatible_theta_struct)
#define gluing_compute_weil_gluing                                      MIKE_NAMESPACE(gluing_compute_weil_gluing)
#define gluing_diag_isogeny_compute                                     MIKE_NAMESPACE(gluing_diag_isogeny_compute)
#define gluing_diag_isogeny_eval                                        MIKE_NAMESPACE(gluing_diag_isogeny_eval)
#define gluing_evaluate_weil_gluing                                     MIKE_NAMESPACE(gluing_evaluate_weil_gluing)
#define gluing_from_couple_dim2_to_dim4_compatible_with_isogeny         MIKE_NAMESPACE(gluing_from_couple_dim2_to_dim4_compatible_with_isogeny)
#define gluing_mike_compute                                             MIKE_NAMESPACE(gluing_mike_compute)
#define gluing_mike_eval                                                MIKE_NAMESPACE(gluing_mike_eval)
#define gluing_special_inv                                              MIKE_NAMESPACE(gluing_special_inv)
#define scholten_eval                                                   MIKE_NAMESPACE(scholten_eval)
#define scholten_gluing_eval                                            MIKE_NAMESPACE(scholten_gluing_eval)
#define scholten_jacobian_to_theta                                      MIKE_NAMESPACE(scholten_jacobian_to_theta)
#define scholten_precompute                                             MIKE_NAMESPACE(scholten_precompute)
#define scholten_variable_compute                                       MIKE_NAMESPACE(scholten_variable_compute)
#define schoten_gluing_compute_and_verify                               MIKE_NAMESPACE(schoten_gluing_compute_and_verify)

// Namespacing symbols exported from invariants.c:
#undef mike_absolute_invariants

#define mike_absolute_invariants                                        MIKE_NAMESPACE(mike_absolute_invariants)

// Namespacing symbols exported from isog_chains.c:
#undef ec_iso_isogeny_2chain
#undef ec_iso_isogeny_2chain_with_strategy
#undef ec_iso_xeval_4
#undef ec_iso_xisog_4
#undef ec_xDBL_A24
#undef ec_xz_is_four_torsion
#undef ec_xz_is_two_torsion
#undef ec_xz_order_even2

#define ec_iso_isogeny_2chain                                           MIKE_NAMESPACE(ec_iso_isogeny_2chain)
#define ec_iso_isogeny_2chain_with_strategy                             MIKE_NAMESPACE(ec_iso_isogeny_2chain_with_strategy)
#define ec_iso_xeval_4                                                  MIKE_NAMESPACE(ec_iso_xeval_4)
#define ec_iso_xisog_4                                                  MIKE_NAMESPACE(ec_iso_xisog_4)
#define ec_xDBL_A24                                                     MIKE_NAMESPACE(ec_xDBL_A24)
#define ec_xz_is_four_torsion                                           MIKE_NAMESPACE(ec_xz_is_four_torsion)
#define ec_xz_is_two_torsion                                            MIKE_NAMESPACE(ec_xz_is_two_torsion)
#define ec_xz_order_even2                                               MIKE_NAMESPACE(ec_xz_order_even2)

// Namespacing symbols exported from isogeny_chain_weil.c:
#undef isogeny_chain_weil

#define isogeny_chain_weil                                              MIKE_NAMESPACE(isogeny_chain_weil)

// Namespacing symbols exported from isogeny_dim2.c:
#undef isogeny_dim2_compute_codomain
#undef isogeny_dim2_eval

#define isogeny_dim2_compute_codomain                                   MIKE_NAMESPACE(isogeny_dim2_compute_codomain)
#define isogeny_dim2_eval                                               MIKE_NAMESPACE(isogeny_dim2_eval)

// Namespacing symbols exported from isogeny_dim4.c:
#undef isogeny_weil_compute_codomain
#undef isogeny_weil_eval

#define isogeny_weil_compute_codomain                                   MIKE_NAMESPACE(isogeny_weil_compute_codomain)
#define isogeny_weil_eval                                               MIKE_NAMESPACE(isogeny_weil_eval)

// Namespacing symbols exported from keygen.c:
#undef protocols_keygen

#define protocols_keygen                                                MIKE_NAMESPACE(protocols_keygen)

// Namespacing symbols exported from mem.c:
#undef mike_secure_clear
#undef mike_secure_free

#define mike_secure_clear                                               MIKE_NAMESPACE_GENERIC(mike_secure_clear)
#define mike_secure_free                                                MIKE_NAMESPACE_GENERIC(mike_secure_free)

// Namespacing symbols exported from mike.c:
#undef mike_exchange
#undef mike_keypair

#define mike_exchange                                                   MIKE_NAMESPACE(mike_exchange)
#define mike_keypair                                                    MIKE_NAMESPACE(mike_keypair)

// Namespacing symbols exported from mike_dim4_chain.c:
#undef mike_isogeny_chain_dim4

#define mike_isogeny_chain_dim4                                         MIKE_NAMESPACE(mike_isogeny_chain_dim4)

// Namespacing symbols exported from normalize.c:
#undef ec_compute_montgomery_coefficient
#undef ec_find_max_coefficient_in_list
#undef ec_fp2_frob
#undef ec_fp2_less_than
#undef ec_fp2_mul_by_i
#undef ec_normalize_montgomery
#undef ec_theta_to_montgomery

#define ec_compute_montgomery_coefficient                               MIKE_NAMESPACE(ec_compute_montgomery_coefficient)
#define ec_find_max_coefficient_in_list                                 MIKE_NAMESPACE(ec_find_max_coefficient_in_list)
#define ec_fp2_frob                                                     MIKE_NAMESPACE(ec_fp2_frob)
#define ec_fp2_less_than                                                MIKE_NAMESPACE(ec_fp2_less_than)
#define ec_fp2_mul_by_i                                                 MIKE_NAMESPACE(ec_fp2_mul_by_i)
#define ec_normalize_montgomery                                         MIKE_NAMESPACE(ec_normalize_montgomery)
#define ec_theta_to_montgomery                                          MIKE_NAMESPACE(ec_theta_to_montgomery)

// Namespacing symbols exported from theta_struct_dim2.c:
#undef theta_dim2_DBL
#undef theta_dim2_DBL_iter
#undef theta_dim2_diff_ADD
#undef theta_struct_dim2_arith_precomp

#define theta_dim2_DBL                                                  MIKE_NAMESPACE(theta_dim2_DBL)
#define theta_dim2_DBL_iter                                             MIKE_NAMESPACE(theta_dim2_DBL_iter)
#define theta_dim2_diff_ADD                                             MIKE_NAMESPACE(theta_dim2_diff_ADD)
#define theta_struct_dim2_arith_precomp                                 MIKE_NAMESPACE(theta_struct_dim2_arith_precomp)

// Namespacing symbols exported from theta_struct_dim4.c:
#undef theta_DBL_iter_weil
#undef theta_DBL_weil
#undef theta_struct_weil_arith_precomp

#define theta_DBL_iter_weil                                             MIKE_NAMESPACE(theta_DBL_iter_weil)
#define theta_DBL_weil                                                  MIKE_NAMESPACE(theta_DBL_weil)
#define theta_struct_weil_arith_precomp                                 MIKE_NAMESPACE(theta_struct_weil_arith_precomp)


#endif

