#!/usr/bin/env python3

from sage.all import *
proof.all(False)

from cformat import FpEl, Object, ObjectFormatter
from torsion_basis import even_torsion_basis_E0, EM_basis_E0

def theta_to_montgomery(a,b):
    return 2*(a**4+b**4)/(a**4-b**4)

def get_montgomery_coeffs(P):
    x = P.x()
    z = 1
    i = x.parent().gen()
    p = x.parent().characteristic()

    a1, b1 = x, z
    a2, b2 = x+z, x-z
    a3, b3 = x+i*z, i*x+z

    A1 = theta_to_montgomery(a1,b1)
    A2 = theta_to_montgomery(a2,b2)
    A3 = theta_to_montgomery(a3,b3)

    L_mont = [A1,A2,A3]
    for i in range(3):
        L_mont.append(-L_mont[i])
    for i in range(6):
        L_mont.append(L_mont[i]**p)

    return L_mont

def e0_basis_files(prefix, f, c, sk_chain_length, d_word_params):
    p = ZZ(c*2**f-1)

    Fp2 = GF((p,2), modulus=[1,0,1], name='i')
    E0 = EllipticCurve(Fp2, [1, 0])

    P, Q = even_torsion_basis_E0(E0, f)
    E0p = EllipticCurve(Fp2, [-1, 0])
    E0M, PM, QM = EM_basis_E0(E0p, f)
    PmQM = PM-QM
    PMsk = (2**(f-sk_chain_length-2))*PM;
    QMsk = (2**(f-sk_chain_length-2))*QM;
    PmQMsk = PMsk-QMsk
    PpQMsk = PMsk+QMsk

    def Fp2_to_list(el):
        return [FpEl(int(c), p, d_word_params, True) for c in Fp2(el)]

    objs = ObjectFormatter([
        Object('fp2_t', 'BASIS_E0_PX', Fp2_to_list(P.x())),
        Object('fp2_t', 'BASIS_E0_QX', Fp2_to_list(Q.x())),
        Object('fp2_t', 'BASIS_EM_PX', Fp2_to_list(PM.x())),
        Object('fp2_t', 'BASIS_EM_QX', Fp2_to_list(QM.x())),
        Object('fp2_t', 'BASIS_EM_PMQX', Fp2_to_list(PmQM.x())),
        Object('fp2_t', 'BASIS_sk_PX', Fp2_to_list(PMsk.x())),
        Object('fp2_t', 'BASIS_sk_QX', Fp2_to_list(QMsk.x())),
        Object('fp2_t', 'BASIS_sk_PMQX', Fp2_to_list(PmQMsk.x())),
    ])

    with open(f'{prefix}include/e0_basis.h','w') as hfile:
        with open(f'{prefix}e0_basis.c','w') as cfile:
            print(f'#include <fp2.h>', file=hfile)
            print(f'#include <e0_basis.h>', file=cfile)

            objs.header(file=hfile)
            objs.implementation(file=cfile)
