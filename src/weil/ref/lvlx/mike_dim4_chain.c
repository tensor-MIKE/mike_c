#include <weil.h>


uint32_t
mike_isogeny_chain_dim4(mike_abs_invariants_t *mike_abs_inv,
                            const mike_chain_basis_t *mike_basis,
                            const mike_gluing_basis_t *gluing_basis,
                            const ec_curve_t *E,
                            const unsigned int length)
{
    mike_gluing_t gluing;
    theta_point_dim4_t T3, T4;
    theta_struct_weil_t weil_codomain;
    theta_struct_weil_precomp_t weil_codomain_precomp;
    theta_point_weil_t final_codomain;
    uint32_t res, res2 = -1;

    // Compute the gluing
    res = gluing_mike_compute(&gluing, E, gluing_basis);

    if(!res)
        return 0; 
    // Push kernel point through the gluing
    gluing_mike_eval(&T3, &mike_basis->T3_1, &mike_basis->T3_2, &gluing);
    gluing_mike_eval(&T4, &mike_basis->T4_1, &mike_basis->T4_2, &gluing);

    // Copy gluing codomain
    theta_copy_weil_struct(&weil_codomain, &gluing.weil_gluing.weil_codom);
    theta_copy_weil_struct_precomp(&weil_codomain_precomp, &gluing.weil_gluing.weil_precomp );
    
    // Compute chain
    isogeny_chain_weil(&weil_codomain, &weil_codomain_precomp,  &T3, &T4, length - 3);

    // Compute dual theta null point
    theta_invert_weil(&final_codomain, &weil_codomain.inv_dual_null_point);

    // Compute invariants
    mike_absolute_invariants(mike_abs_inv, &final_codomain);

    // Check invariant are not all zeros
    for(uint8_t i =0; i < 4; i++)
        res2 &= fp_is_zero(&(*mike_abs_inv)[i]); 

    // 0 if res = 0 or res2 = -1.  
    return res & -(res2 + 1);
}