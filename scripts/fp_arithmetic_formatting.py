#!/usr/bin/env python3

from sage.all import *
proof.all(False)

from re import match, sub
from textwrap import dedent

# Fixing static inline declarations that are wring in modarith generated files
def fix_syntax_field_file(input_file):
    pattern = r'(static)\s+(\w+(?:\s*\*)*)\s+(inline)'

    def replace_func(match):
        static_keyword = match.group(1)
        type_name = match.group(2)
        inline_keyword = match.group(3)
        
        # Return the corrected order: static inline <type>
        return f'{static_keyword} {inline_keyword} {type_name}'

    def fix_static_inline_declarations(content):
        # Apply the replacement
        corrected_content = sub(pattern, replace_func, content)
        return corrected_content

    # Read the input file
    with open(input_file, 'r') as f:
        content = f.read()
        
    # Fix the declarations
    corrected_content = fix_static_inline_declarations(content)
        
    # Write back to the same file
    with open(input_file, 'w') as f:
        f.write(corrected_content)

def int_to_montgemery_fp_const(x, p, Nlimbs, Radix):
    R = pow(2, Nlimbs * Radix, p)
    el = (x * R) % p
    vs = [(int(el) >> Radix * i) % 2**Radix for i in range(Nlimbs)]
    return "{" + ",\n".join(map(hex, vs)) + "}"

# Adding generic API functions to call from the header file
def write_field_file(p,write=True,name=None,wl=64):
    d_word_params = {}

    params = ["Wordlength", "Nlimbs", "Radix", "Nbits", "Nbytes"]

    if write:
        filename = "field.c"
    else:
        filename = "../src/gf/ref/" + name + "/fp_" + name + "_" + str(wl) + ".in.c"
    for line in open(filename).readlines():
        for param in params:
            if m := match(f"^#define {param} ([0-9]*)$", line):
                d_word_params[param] = int(m.group(1))

    Nlimbs = d_word_params["Nlimbs"]
    Nbytes = d_word_params["Nbytes"]
    Radix = d_word_params["Radix"]

    lines = []

    lines += [""]
    lines += ["/******************************************************************************"]
    lines += ["MIKE API functions calling generated code above"]
    lines += ["******************************************************************************/"]
    lines += [""]

    lines += [f"#include <fp.h>\n"]
    lines += [f"const digit_t ZERO[NWORDS_FIELD] = {int_to_montgemery_fp_const(0, p, Nlimbs, Radix)};"]
    lines += [f"const digit_t ONE[NWORDS_FIELD] = {int_to_montgemery_fp_const(1, p, Nlimbs, Radix)};"]
    lines += ["// Montgomery representation of 2^-1"]
    lines += [
        f"static const digit_t TWO_INV[NWORDS_FIELD] = {int_to_montgemery_fp_const(inverse_mod(2, p), p, Nlimbs, Radix)};"
    ]
    lines += ["// Montgomery representation of 3^-1"]
    lines += [
        f"static const digit_t THREE_INV[NWORDS_FIELD] = {int_to_montgemery_fp_const(inverse_mod(3, p), p, Nlimbs, Radix)};"
    ]

    lines += [
        dedent(
            """
        void
        fp_set_small(fp_t *x, const digit_t val)
        {
            modint((int)val, *x);
        }

        void
        fp_mul_small(fp_t *x, const fp_t *a, const uint32_t val)
        {
            modmli(*a, (int)val, *x);
        }

        void
        fp_set_zero(fp_t *x)
        {
            modzer(*x);
        }

        void
        fp_set_one(fp_t *x)
        {
            modone(*x);
        }

        uint32_t
        fp_is_equal(const fp_t *a, const fp_t *b)
        {
            return -(uint32_t)modcmp(*a, *b);
        }

        uint32_t
        fp_is_zero(const fp_t *a)
        {
            return -(uint32_t)modis0(*a);
        }

        void
        fp_copy(fp_t *out, const fp_t *a)
        {
            modcpy(*a, *out);
        }

        void
        fp_cswap(fp_t *a, fp_t *b, uint32_t ctl)
        {
            modcsw((int)(ctl & 0x1), *a, *b);
        }

        void
        fp_add(fp_t *out, const fp_t *a, const fp_t *b)
        {
            modadd(*a, *b, *out);
        }

        void
        fp_sub(fp_t *out, const fp_t *a, const fp_t *b)
        {
            modsub(*a, *b, *out);
        }

        void
        fp_neg(fp_t *out, const fp_t *a)
        {
            modneg(*a, *out);
        }

        void
        fp_sqr(fp_t *out, const fp_t *a)
        {
            modsqr(*a, *out);
        }

        void
        fp_mul(fp_t *out, const fp_t *a, const fp_t *b)
        {
            modmul(*a, *b, *out);
        }

        void
        fp_inv(fp_t *x)
        {
            modinv(*x, NULL, *x);
        }

        uint32_t
        fp_is_square(const fp_t *a)
        {
            return -(uint32_t)modqr(NULL, *a);
        }

        void
        fp_sqrt(fp_t *a)
        {
            modsqrt(*a, NULL, *a);
        }

        void
        fp_half(fp_t *out, const fp_t *a)
        {
            modmul(TWO_INV, *a, *out);
        }

        void
        fp_exp3div4(fp_t *out, const fp_t *a)
        {
            modpro(*a, *out);
        }

        void
        fp_div3(fp_t *out, const fp_t *a)
        {
            modmul(THREE_INV, *a, *out);
        }

        void
        fp_encode(void *dst, const fp_t *a)
        {
            // little-endian canonical encoding
            int i;
            spint c[Nlimbs]= {0};
            redc(*a, c);
            for (i = 0; i < Nbytes; i++) {
                ((char *)dst)[i] = c[0] & (spint)0xff;
                (void)modshr(8, c);
            }
        }

        uint32_t
        fp_decode(fp_t *d, const void *src)
        {
            // inverse of fp_encode; returns 0xFFFFFFFF iff the input was canonical
            int i;
            spint res;
            const unsigned char *b = src;
            for (i = 0; i < Nlimbs; i++) {
                (*d)[i] = 0;
            }
            for (i = Nbytes - 1; i >= 0; i--) {
                modshl(8, *d);
                (*d)[0] += (spint)b[i];
            }
            res = (spint)-modfsb(*d);
            nres(*d, *d);
            for (i = 0; i < Nlimbs; i++) {
                (*d)[i] &= res;
            }
            return (uint32_t)res;
        }

        void
        fp_decode_reduce(fp_t *d, const void *src, size_t len)
        {
            // Reduce a little-endian byte string of arbitrary length mod p.
            // Radix-256 Horner using only fp_add/fp_set_small: correct for ANY prime,
            // no prime-shape-specific fast reduction required.
            const unsigned char *b = src;
            size_t i;
            fp_set_zero(d);
            for (i = len; i-- > 0;) {
                fp_t t;
                int k;
                for (k = 0; k < 8; k++) {
                    fp_add(d, d, d); // d *= 2  (eight times => d *= 256)
                }
                fp_set_small(&t, (digit_t)b[i]);
                fp_add(d, d, &t);
            }
        }
    """
        )
    ]

    if write:
        with open("field.c", "a+") as file:
            file.writelines([line + "\n" for line in lines])

    return d_word_params

def write_field_wrapper_file(name):
    lines = ["#include <tutil.h>"]
    lines += ["#if RADIX == 32"]
    lines += [f"#include \"fp_{name}_32.in.c\""]
    lines += ["#else"]
    lines += [f"#include \"fp_{name}_64.in.c\""]
    lines += ["#endif"]

    with open(f"../src/gf/ref/{name}/fp_{name}.c", "w") as file:
        file.writelines([line + "\n" for line in lines])
