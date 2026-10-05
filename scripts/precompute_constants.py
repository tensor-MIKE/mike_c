#!/usr/bin/env python3

from sage.all import *
proof.all(False)

import os
from cformat import FpEl, MpEl, Object, ObjectFormatter

def constants_files(prefix,f,c,d_word_params):
    p = ZZ(c*2**f-1)

    # Finite field instance
    Fp = GF(p)
    Fp2 = GF((p,2), modulus=[1,0,1], name='i')

    # A0 Montgomery curve coefficient, same as in torsion_basis.EM_basis_E0
    alpha = 1  # Pick this root to send Q to be above zero
    s = ~Fp2(3 * alpha**2 - 1).sqrt()
    A0 = -3 * alpha * s
    A0 = Fp(A0) 
    #if is_square(A0-2): # to ensure A0-2 is not a square
        #A0 = -A0

    # List of non-quadratic residue mod p
    nqr_list = []
    nqr = Fp(1)
    for i in range(20):
        while is_square(nqr):
            nqr += 1
        nqr_list.append(int(nqr))
        nqr += 1

    defs = {}
    defs["CHARACTERISTIC"] = p
    defs["TORSION_EVEN"] = 2**f
    defs["TORSION_EVEN_POWER"] = f
    defs["TORSION_ODD"] = (p+1)//(2**f)
    defs["NQR_TABLE"] = [FpEl(n,p,d_word_params,True) for n in nqr_list]
    defs["A0"] = FpEl(A0,p,d_word_params,True)

    objs = ObjectFormatter([
        Object('fp_t', 'A0', defs["A0"]),
        Object('fp_t []', 'NQR_TABLE', defs["NQR_TABLE"]),
        #Object('uint16_t', "TORSION_EVEN_POWER", int(defs["TORSION_EVEN_POWER"])),
        Object('digit_t[]', "CHARACTERISTIC", MpEl(defs["CHARACTERISTIC"],p)),
        Object('digit_t[]', "TORSION_ODD", MpEl(defs["TORSION_ODD"],p)),
        Object('digit_t[]', "TORSION_EVEN", MpEl(defs["TORSION_EVEN"],p))
    ])

    with open(f'{prefix}include/constants.h','w') as hfile:
        with open(f'{prefix}constants.c','w') as cfile:
            print(f'#include <fp2.h>', file=hfile)
            print(f'#include <inttypes.h>', file=hfile)
            print(f'#include <constants.h>', file=cfile)

            objs.header(file=hfile)
            objs.implementation(file=cfile)

if __name__ == "__main__":
    from parameters import prefix,p,f,fp_word_sizes
    constants_files(prefix,p,f,fp_word_sizes)