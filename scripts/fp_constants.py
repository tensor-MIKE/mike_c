#!/usr/bin/env python3

def fp_constants_hfile(name,d_word_params):
    with open(f"../src/precomp/ref/{name}/include/fp_constants.h", "w") as hfile:
        hfile.write("#if RADIX == 32\n")

        hfile.write("#if defined(MIKE_GF_IMPL_SAT64)\n")
        hfile.write(f"#define NWORDS_FIELD {(d_word_params[32]["Nbits"]+31)//32}\n")
        hfile.write("#else\n")
        hfile.write(f"#define NWORDS_FIELD {d_word_params[32]["Nlimbs"]}\n")
        hfile.write("#endif\n")

        hfile.write(f"#define NWORDS_ORDER {(d_word_params[32]["Nbits"]+31)//32}\n")

        hfile.write("#elif RADIX == 64\n")

        hfile.write("#if defined(MIKE_GF_IMPL_SAT64)\n")
        hfile.write(f"#define NWORDS_FIELD {(d_word_params[64]["Nbits"]+63)//64}\n")
        hfile.write("#else\n")
        hfile.write(f"#define NWORDS_FIELD {d_word_params[64]["Nlimbs"]}\n")
        hfile.write("#endif\n")

        hfile.write(f"#define NWORDS_ORDER {(d_word_params[64]["Nbits"]+63)//64}\n")

        hfile.write("#endif\n")
        # LOG2P, BITS are not used

