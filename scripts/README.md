# Parameter generation scripts

This folder contains prime generation scripts to generate new MIKE parameters.

## Requirements

- Python 3.10
- Sagemath 3.8

For sage, one can either use a sage environment activated via `conda activate sage` or `conda activate <path to your sage installation>`, or use sage directly. 

External libraries added as submodules in the `external` folder:
- [Modarith](https://github.com/mcarrickscott/modarith), for prime field arithmetic generation.
- [Addchain](https://github.com/mmcloughlin/addchain), as a dependency of modarith.

These submodules can be downloaded from an already cloned repository with:
`git submodule init`
`git submodule update`

## Generating pure C arithmetic

All MIKE primes already in the source file (`src`) were generated using `sage` (which can be replaced by `python` in an active `sage` environment in all following commands):

`sage generate_fp.py -all [-no]`

To generate constants for a specific prime, type:

`sage generate_fp.py -f=<2 power> -c=<cofactor> [-n=<prime name>] [-cons] [-no]`

If these commands fail with `sage`, you may run them with `python` instead (just replace `sage` by `python`).

The `-n` option is the short form of the `--name` option that is set to `pf_c` by default, where $$p=c\cdot 2^f-1$$ is the defined prime.

The `-cons` option is the short form of the `--conservative` option, that stands for "conservative security estimate" according to Benjamin Wesolowski's attack ($\log_2(p)=3\lambda$ where $\lambda$ is the security level in bits). By default, this option is not active.

The `-no` option is the short form of the `--no_arith_gen` option, that can only be used when there is already arithmetic generated for a given prime to avoid generating it again but only changing some constants. This option is mainly useful for developers.  

Example: for the non-conservative NIST level 1 prime `p308_633`

`sage generate_fp.py -f=308 -c=633`

or

`python generate_fp.py -f=308 -c=633`

## Generating sat64 assembly arithmetic

The saturated 64-bit GF(p) backend (`src/gf/sat64/<pe_c>/{fp.c, broadwell/fp_asm.S, arm64/fp_asm.S}`) is generated for all primes by running `./sat64_generator/generate_all.sh` (requires only Python 3). The C in `fp.c` is shared by both architectures; only `fp_mul`, `fp_sqr`, `fp2_mul` and `fp2_sqr` are assembly, and the build picks `broadwell/` or `arm64/` from `-DMIKE_BUILD_TYPE`.

We recommend building the arm64 backend with clang as it compiles the C field additions and subtractions to a single carry chain, while GCC does not, so its `fp_add` and `fp_sub` are much slower. This could be avoided by hand-writing these functions in assembly too, but this is not done at the moment.

## Rerunning only file generation for constants in src/precomp


This can be achieved by running `sage generate_fp.py -all -no` or `python generate_fp.py -all -no` in the currect directory.



