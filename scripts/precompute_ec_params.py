#!/usr/bin/env python3

from cformat import Object, ObjectFormatter

def ec_params_files(prefix, f, c):
    obj_cof = ObjectFormatter(
        [
            Object('digit_t[]', 'p_cofactor_for_2f', [c]),
        ]
    )

    with open(f"{prefix}include/ec_params.h", "w") as hfile:
        with open(f"{prefix}ec_params.c", "w") as cfile:
            hfile.write('#ifndef EC_PARAMS_H\n')
            hfile.write('#define EC_PARAMS_H\n')
            hfile.write('\n')

            hfile.write('#include <fp.h>\n')
            cfile.write('#include <ec_params.h>\n')
            hfile.write('\n')

            hfile.write(f'#define TORSION_EVEN_POWER {f}\n')
            hfile.write('\n')

            hfile.write('// p+1 divided by the power of 2\n')
            obj_cof.header(file=hfile)
            obj_cof.implementation(file=cfile)
            hfile.write(f'#define P_COFACTOR_FOR_2F_BITLENGTH {c.bit_length()}\n')
            hfile.write('\n')
            cfile.write('\n')

            hfile.write('#endif\n')