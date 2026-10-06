#!/usr/bin/env python3

from sage.all import *
proof.all(False)

import os
import argparse
from shutil import which

from global_constants import hd_margin, shared_bytes
from fp_arithmetic_formatting import fix_syntax_field_file, write_field_file, write_field_wrapper_file
from precompute_constants import constants_files
from fp_constants import fp_constants_hfile
from precompute_ec_params import ec_params_files
from precompute_encoded_sizes import encoded_sizes_hfile
from precompute_e0_basis import e0_basis_files
from fp_cmakefiles import write_field_cmakelists
from mike_parameters import d_mike_parameters
from fp_parameters import compute_constants_at_fp_gen

def build_from_parameters(f,c,name, write=True,conservative=False,check_primality=True):
    d_word_params = {}

    p = ZZ(c*2**f-1)
    if check_primality:
        if not is_prime(p):
            raise ValueError(f"{name} is not a prime.")

    # Writing fp_name.c files in gf/ref
    os.system("mkdir -p ../src/gf/ref/" + name)
    for wl in [32,64]:
        if write:
            os.system(f"PATH=$PATH:external/addchain/bin python external/modarith/monty.py {wl} {p} > /dev/null 2>&1")
            fix_syntax_field_file('field.c')
        d_word_params[wl] = write_field_file(p, write, name=name ,wl=wl)
        if write:
            os.system("mv -v field.c ../src/gf/ref/" + name + "/fp_" + name + "_" + str(wl) + ".in.c")
    write_field_wrapper_file(name)
    if write:
        os.system("rm -v time.c")
        os.system("rm -v ac.txt")

    # Arithmetic Generation

    prefix = f'../src/precomp/ref/{name}/'
    sk_length, security_bits = compute_constants_at_fp_gen(f, c, hd_margin, conservative=conservative)
    # Writing CMakeLists.txt wherever necessary in modules (gf, ec, theta, weil, nike)
    os.system(f"mkdir -p ../src/precomp/ref/{name}")
    os.system(f"mkdir -p ../src/precomp/ref/{name}/include")
    write_field_cmakelists(name, p, sk_length, security_bits, hd_margin)
    # Writing precomp/ref/<name>/fp_constants.c
    fp_constants_hfile(name, d_word_params)

    # Precomputations

    # Writing precomp/ref/<name>/constants.c and precomp/ref/<name>/include/constants.h
    constants_files(prefix, f, c, d_word_params)
    # Writing precomp/ref/<name>/ec_params.c and precomp/ref/<name>/include/ec_params.h
    ec_params_files(prefix, f, c)
    # Writing precomp/ref/<name>/include/encoded_sizes.h
    encoded_sizes_hfile(prefix, f, c, security_bits, sk_length, hd_margin, shared_bytes)
    # Writing precomp/ref/<name>/e0_basis.c and precomp/ref/<name>/include/e0_basis.h
    e0_basis_files(prefix, f, c, sk_length, d_word_params)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()


    #parser.add_argument("-p", "--prime")
    parser.add_argument("-f")
    parser.add_argument("-c")
    parser.add_argument("-n", "--name")
    parser.add_argument("-cons", "--conservative", action="store_true")
    parser.add_argument("-no", "--no_arith_gen", action="store_true")
    parser.add_argument("-all", "--all_params", action="store_true")

    args = parser.parse_args()

    # The no_arith_gen option (no build of field files with modarith)
    write = not args.no_arith_gen

    # `addchain` must be compiled, and so cannot be shipped as-is
    if write:
        addchain_binary = which("addchain", path=os.environ['PATH'] + ":external/addchain/bin")
        if addchain_binary:
            print(f"Found addchain binary at {addchain_binary}")
        if not addchain_binary:
            print("Cannot find `addchain` in `$PATH`. This is necessary for `modarith`")
            print("Installing now to external/addchain/bin")
            os.system("cd external/addchain && bash install.sh")

    if args.all_params:
        for name in d_mike_parameters:
            print(f'Generating files for {name}')
            f = ZZ(d_mike_parameters[name]['f'])
            c = ZZ(d_mike_parameters[name]['c'])
            conservative = d_mike_parameters[name]['conservative']
            build_from_parameters(f,c,name,write,conservative,False)
    else:
        # Bad habit to use `eval`, could be replaced with https://stackoverflow.com/a/69540962
        f = eval(args.f)
        c = eval(args.c)
        # Setting the name if not already done
        if args.name is None:
            name = "p"+str(f)+"_"+str(c)
        print(f'Generating files for {name}')
        # p = eval(args.prime)
        # Main build of files
        build_from_parameters(f,c,name,write,args.conservative)
