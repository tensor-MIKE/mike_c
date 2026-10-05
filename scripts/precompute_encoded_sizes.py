#!/usr/bin/env python3

from sage.all import *
proof.all(False)

def encoded_sizes_hfile(prefix, f, c, security_bits, sk_length, hd_margin, shared_bytes):
    p = ZZ(c*2**f-1)
    logp = ceil(log(p, 2))
    defs = dict()

    TORSION_2POWER_BYTES = (f + 7) // 8
    SK_LENGTH = sk_length
    SK_LENGTH_BYTES = (SK_LENGTH + 9) // 8

    fpsz = (logp + 7)//8
    fp2sz = 2 * fpsz
    defs['SECURITY_BITS'] = security_bits
    defs['FP_ENCODED_BYTES'] = fpsz
    defs['FP2_ENCODED_BYTES'] = fp2sz
    defs['EC_CURVE_ENCODED_BYTES'] = fp2sz  # just the A

    defs['PUBLICKEY_BYTES'] = defs['EC_CURVE_ENCODED_BYTES']
    defs['SECRETKEY_BYTES'] = SK_LENGTH_BYTES
    defs['SECRETKEY_BITS'] = SK_LENGTH
    defs['HD_MARGIN'] = hd_margin
    defs['SECRETKEY_CHAIN_LENGTH'] = sk_length
    defs['SHARED_BYTES'] = shared_bytes

    size_privkey = defs['SECRETKEY_BYTES']
    size_pubkey = defs['PUBLICKEY_BYTES']
    sk_words_32 = (size_privkey+3)//4
    sk_words_64 =  (size_privkey+7)//8

    sk_words = "#if RADIX == 32\n#define SECRETKEY_WORDS "+ str(sk_words_32) + "\n#elif RADIX == 64\n#define SECRETKEY_WORDS " +str(sk_words_64)+"\n#endif\n"

    with open(f'{prefix}include/encoded_sizes.h','w') as hfile:
        for k,v in defs.items():
            v = ZZ(v)
            print(f'#define {k} {v}', file=hfile)
        print(sk_words, file=hfile)


if __name__ == "__main__":
    from parameters import prefix, p, f, security_bits, sk_length, hd_margin, shared_bytes
    encoded_sizes_hfile(prefix, p, f, security_bits, sk_length, hd_margin, shared_bytes)