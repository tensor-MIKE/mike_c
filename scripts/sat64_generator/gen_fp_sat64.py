#!/usr/bin/env python3
"""
Generator of the saturated 64-bit GF(p) backend: the portable C from gen_fp.py
with fp_mul/fp_sqr replaced by per-architecture assembly from
gen_fp_asm_broadwell.py and gen_fp_asm_arm64.py.

Writes <source-dir>/fp.c, <source-dir>/broadwell/fp_asm.S and
<source-dir>/arm64/fp_asm.S.

Usage:
    python gen_fp_sat64.py "633 * 2**308 - 1" --source-dir ../../src/gf/sat64/p308_633
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gen_fp import FieldGenerator, parse_int
import gen_fp_asm_arm64
import gen_fp_asm_broadwell

ASM_GENERATORS = {
    "broadwell": gen_fp_asm_broadwell,
    "arm64": gen_fp_asm_arm64,
}


class Sat64FieldGenerator(FieldGenerator):
    def generate_mul(self) -> str:
        return "// fp_mul: implemented in <arch>/fp_asm.S\n"

    def generate_sqr(self) -> str:
        return "// fp_sqr: implemented in <arch>/fp_asm.S\n"

    def generate(self, source_dir: Path):
        super().generate(source_dir)
        for arch, gen in ASM_GENERATORS.items():
            arch_dir = source_dir / arch
            arch_dir.mkdir(parents=True, exist_ok=True)
            (arch_dir / "fp_asm.S").write_text(gen.generate(self.p))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "prime",
        help="prime as decimal or 0x-prefixed integer, or a Python expression like '633 * 2**308 - 1'",
    )
    ap.add_argument("--source-dir", required=True)
    args = ap.parse_args()

    Sat64FieldGenerator(parse_int(args.prime)).generate(Path(args.source_dir))


if __name__ == "__main__":
    main()
