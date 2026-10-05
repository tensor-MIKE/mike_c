#!/bin/sh
# Generate the sat64 GF(p) backend (fp.c + {broadwell,arm64}/fp_asm.S) for every MIKE prime.
set -e
cd "$(dirname "$0")"

gen() {
    python3 gen_fp_sat64.py "$2" --source-dir "../../src/gf/sat64/$1"
}

gen p308_633 "633 * 2**308 - 1"
gen p374_117 "117 * 2**374 - 1"
gen p474_593 "593 * 2**474 - 1"
gen p566_77  "77 * 2**566 - 1"
gen p628_317 "317 * 2**628 - 1"
gen p758_41  "41 * 2**758 - 1"
