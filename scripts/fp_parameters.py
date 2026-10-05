#!/usr/bin/env python3

from sage.all import *
proof.all(False)

def factor_pp1(p, factor_l=True):
    L = list(factor(p + 1))
    f = L[0][1]
    c = (p + 1) // (2**f)
    return c, f

def  compute_constants_at_fp_gen(f,c,hd_margin=2,conservative=False):
    p = ZZ(c*2**f-1)

    logp = ceil(log(p, 2))

    if conservative:
        security_bits = round(logp/(3*64))*64
        sk_length = f-hd_margin
    else:
        security_bits = round(8*logp/(21*64))*64
        sk_length = 7*security_bits//4+4

    return sk_length, security_bits

