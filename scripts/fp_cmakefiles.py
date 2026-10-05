#!/usr/bin/env python3

from sage.all import *
proof.all(False)

import os

def write_field_cmakelists(name, p, sk_length, security_bits, hd_margin):
    # Writing gf/ref/<name>/CMakeLists.txt
    lines = ["set(SOURCE_FILES_GF_SPECIFIC"]
    lines += [f"\tfp_{name}.c"]
    lines += [")\n"]
    lines += ["include(../lvlx.cmake)"]
    with open(f"../src/gf/ref/{name}/CMakeLists.txt", "w") as file:
        file.writelines([line + "\n" for line in lines])

    # Writing gf/ref/<name>/test/CMakeLists.txt
    os.system(f"mkdir -p ../src/gf/ref/{name}/test")
    lines = ["include(../../lvlx_test.cmake)"]
    with open(f"../src/gf/ref/{name}/test/CMakeLists.txt", "w") as file:
        file.writelines([line + "\n" for line in lines])

    # Modules with tests
    for module in ['ec','theta','weil','nike']:
        os.system(f"mkdir -p ../src/{module}/ref/{name}")
        lines = ["include(../lvlx.cmake)"]
        with open(f"../src/{module}/ref/{name}/CMakeLists.txt", "w") as file:
            file.writelines([line + "\n" for line in lines])

        os.system(f"mkdir -p ../src/{module}/ref/{name}/test")
        lines = ["include(../../lvlx_test.cmake)"]
        with open(f"../src/{module}/ref/{name}/test/CMakeLists.txt", "w") as file:
            file.writelines([line + "\n" for line in lines])

    # Modules without tests
    for module in ['precomp']:
        os.system(f"mkdir -p ../src/{module}/ref/{name}")
        lines = ["include(../lvlx.cmake)"]
        with open(f"../src/{module}/ref/{name}/CMakeLists.txt", "w") as file:
            file.writelines([line + "\n" for line in lines])

    # make parameter files
    for module in ['precomp']:
        os.system(f"mkdir -p ../src/{module}/ref/{name}/include")
        lines = [f"name = {name}",f"p = {p}",f"sk_length = {sk_length}",f"security_bits = {security_bits}"]
        with open(f"../src/{module}/ref/{name}/mike_parameters.txt", "w") as file:
            file.writelines([line + "\n" for line in lines])
