#!/usr/bin/env python3
# This file was generated using a LLM which had information of both:
# 1. The pure C implementation written for an fp2 library in C
#    https://github.com/GiacomoPope/fp2_c
# 2. The AArch64 generator written for SQIsign
#    https://github.com/SQISign/the-sqisign/blob/main/scripts/gen_fp/gen_fp_asm_arm64.py
# The instruction helpers, the multiply-accumulate macros and the Montgomery
# loop are taken from there; the routines built on them are rewritten to keep
# every value fully reduced (see below).

"""
Generator of AArch64 assembly for fp_mul, fp_sqr, fp2_mul and fp2_sqr, for
Montgomery-friendly primes p = c * 2^t - 1.

Matches gen_fp.py's fp_t representation exactly (radix
2^64, N = ceil(bits/64) limbs, R = 2^(64*N) mod p) and is fully canonical:
unlike the SQIsign AArch64 generator, which leaves outputs lazily reduced in
[0, 2^bits) and needs 4 spare top-limb bits to stay bounded, every routine
here takes inputs in [0, p) and returns outputs in [0, p).

Because t >= 64*(N-1), the low N-1 limbs of p are all 0xFFFFFFFFFFFFFFFF,
forcing mu = -p^-1 mod 2^64 = 1, and p+1 has a single non-zero 64-bit limb,
so each Montgomery reduction step is one 64x64 multiply into the top two
accumulator limbs. A Montgomery product of x, y with x*y < 4p^2 is below
4p^2/R + p, which is < 2p when R >= 4p, i.e. with at least two spare bits;
every routine relies on this and finishes each product with one masked
subtraction of p:

    fp_mul   a x b                    a, b < p
    fp2_mul  a0 x b0, a1 x b1,        a, b < p
             (a0+a1) x (b0+b1)        sums < 2p
             c0 = P0 - P1, c1 = P2 - P0 - P1 (mod p)
    fp2_sqr  c0 = (a0+a1) x (a0-a1+p) both < 2p
             c1 = a1 x 2*a0           2*a0 < 2p

fp_sqr is fp_mul with b = a. The multiplicand is register-resident (N
registers), the accumulator takes N+1 and the scalar and a temp one each:
2N+3 out of the 27 usable registers (x18 is platform-reserved, x29 is left
alone). fp2_mul/fp2_sqr also stream a scalar operand through x1, which caps
them at 2N+3 <= 26; for N = 12 they spill the output pointer x0 to the stack
and use it as the 27th register.

Usage:
    python gen_fp_asm_arm64.py "633 * 2**308 - 1" -o fp_asm.S
"""

import argparse
from pathlib import Path

MASK64 = (1 << 64) - 1

# x0/x2 stay pointers where live, x18 is platform-reserved, x29 is untouched.
# x30 (saved) and x1 (dead once the operands are loaded) close the pool.
POOL = (
    ["x%d" % i for i in range(3, 18)]
    + ["x%d" % i for i in range(19, 29)]
    + ["x30", "x1"]
)
CALLEE_SAVED = tuple("x%d" % i for i in range(19, 29)) + ("x30",)


class Params:
    def __init__(self, p: int):
        self.p = p
        self.bits = p.bit_length()
        self.N = (self.bits + 63) // 64
        self.W = 64 * self.N
        self.spare = self.W - self.bits

        if self.N < 3:
            raise ValueError("primes below 3 limbs are not supported")
        modulus = [(p >> (64 * i)) & MASK64 for i in range(self.N)]
        if any(limb != MASK64 for limb in modulus[:-1]):
            raise ValueError(
                "the low N-1 limbs of p must be 0xFFFFFFFFFFFFFFFF, i.e. "
                "p = c*2^t - 1 with t >= 64*(N-1)"
            )
        if 4 * p > 1 << self.W:
            raise ValueError(
                "p needs at least two spare bits (R >= 4p) for the Montgomery "
                "products of unreduced sums to stay below 2p; got %d" % self.spare
            )

        # The single non-zero limb of p+1, and the top limb of p.
        self.pp1_top = (p + 1) >> (64 * (self.N - 1))
        self.p_top = modulus[-1]
        assert self.pp1_top == self.p_top + 1 and self.pp1_top <= MASK64

        n = self.N
        # fp_mul: multiplicand A, accumulator Z, scalar BI, temp T0.
        # The reduction borrows BI as its second temp (dead there).
        self.A = POOL[:n]
        self.Z = POOL[n : 2 * n + 1]
        self.BI = POOL[2 * n + 1]
        self.T0 = POOL[2 * n + 2]

        # fp2_mul/fp2_sqr: x1 streams a scalar operand, so it leaves the
        # pool; when that is one register short, x0 is spilled and joins it.
        pool2 = [r for r in POOL if r != "x1"]
        self.spill_out = 2 * n + 3 > len(pool2)
        if self.spill_out:
            pool2.append("x0")
        if 2 * n + 3 > len(pool2):
            raise ValueError(
                "p is too large: %d limbs need %d registers" % (n, 2 * n + 3)
            )
        self.A2 = pool2[:n]
        self.Z2 = pool2[n : 2 * n + 1]
        self.BI2 = pool2[2 * n + 1]
        self.T02 = pool2[2 * n + 2]


def hexq(v):
    return "0x%016X" % v


def saved_regs(used):
    return [r for r in CALLEE_SAVED if r in used]


def reg_range(regs):
    """Compact display: [x14,x15,x16,x17,x19] -> 'x14..x17,x19'."""
    nums = [int(r[1:]) for r in regs]
    parts, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        if j == i:
            parts.append("x%d" % nums[i])
        elif j == i + 1:
            parts.extend(["x%d" % nums[i], "x%d" % nums[j]])
        else:
            parts.append("x%d..x%d" % (nums[i], nums[j]))
        i = j + 1
    return ",".join(parts)


# --------------------------------------------------------------------------
# Instruction formatting
# --------------------------------------------------------------------------
def ins(mn, args, ind="    "):
    return "%s%-5s %s\n" % (ind, mn, args)


def tail(line, comment):
    """Append an aligned // comment to a formatted instruction line."""
    text = line.rstrip("\n")
    pad = max(37, len(text) + 1)
    return "%s// %s\n" % (text.ljust(pad), comment)


def mov_imm(reg, v, ind="    ", comment=""):
    """Materialise a 64-bit constant with the shortest movz/movn/movk run.
    None of these touch the flags, so they may sit inside a carry chain."""
    v &= MASK64
    ch = [(v >> (16 * i)) & 0xFFFF for i in range(4)]
    nz = [i for i in range(4) if ch[i]]
    nf = [i for i in range(4) if ch[i] != 0xFFFF]
    lsl = lambda i: "" if i == 0 else ", lsl #%d" % (16 * i)
    out = []
    if not nz:
        out.append(ins("mov", "%s, xzr" % reg, ind))
    elif v > (1 << 64) - 0x10001:
        out.append(ins("mov", "%s, #-%d" % (reg, (1 << 64) - v), ind))
    elif len(nf) < len(nz):
        # the immediate is complemented: spell the value out in the comment
        if comment:
            comment += " = %s" % hexq(v)
        i = nf[0]
        out.append(ins("movn", "%s, #0x%04x%s" % (reg, ~ch[i] & 0xFFFF, lsl(i)), ind))
        for i in nf[1:]:
            out.append(ins("movk", "%s, #0x%04x%s" % (reg, ch[i], lsl(i)), ind))
    else:
        i = nz[0]
        out.append(ins("movz", "%s, #0x%04x%s" % (reg, ch[i], lsl(i)), ind))
        for i in nz[1:]:
            out.append(ins("movk", "%s, #0x%04x%s" % (reg, ch[i], lsl(i)), ind))
    if comment:
        out[-1] = tail(out[-1], comment)
    return "".join(out)


def _adr(base, i):
    return "[%s]" % base if i == 0 else "[%s, #%d]" % (base, 8 * i)


def ld_vec(regs, base, index=0, ind="    "):
    """ldp/ldr the given registers from base + 8*index."""
    s, i = [], 0
    while i < len(regs):
        if i + 1 < len(regs):
            s.append(
                ins(
                    "ldp",
                    "%s, %s, %s" % (regs[i], regs[i + 1], _adr(base, index + i)),
                    ind,
                )
            )
            i += 2
        else:
            s.append(ins("ldr", "%s, %s" % (regs[i], _adr(base, index + i)), ind))
            i += 1
    return "".join(s)


def st_vec(regs, base, index=0, ind="    "):
    s, i = [], 0
    while i < len(regs):
        if i + 1 < len(regs):
            s.append(
                ins(
                    "stp",
                    "%s, %s, %s" % (regs[i], regs[i + 1], _adr(base, index + i)),
                    ind,
                )
            )
            i += 2
        else:
            s.append(ins("str", "%s, %s" % (regs[i], _adr(base, index + i)), ind))
            i += 1
    return "".join(s)


# --------------------------------------------------------------------------
# Prologue / epilogue
# --------------------------------------------------------------------------
def _save_pairs(saved):
    pairs = [(saved[k], saved[k + 1]) for k in range(0, len(saved) - 1, 2)]
    odd = saved[-1] if len(saved) % 2 else None
    return pairs, odd


def frame_open(saved, bufbytes=0, buffer_name="", ind="    "):
    """Pre-index stp frame, or sub-sp frame when a scratch buffer is needed."""
    pairs, odd = _save_pairs(saved)
    save_bytes = 16 * ((len(saved) + 1) // 2)
    s = []
    if bufbytes:
        buf = (bufbytes + 15) & ~15
        line = ins("sub", "sp, sp, #%d" % (buf + save_bytes), ind)
        if saved:
            line = tail(
                line, "%d for %s + %d for callee saves" % (buf, buffer_name, save_bytes)
            )
        else:
            line = tail(line, "scratch for %s" % buffer_name)
        s.append(line)
        for i, (r0, r1) in enumerate(pairs):
            s.append(ins("stp", "%s, %s, [sp, #%d]" % (r0, r1, buf + 16 * i), ind))
        if odd:
            s.append(ins("str", "%s, [sp, #%d]" % (odd, buf + 16 * len(pairs)), ind))
    elif len(saved) == 1:
        s.append(ins("str", "%s, [sp, #-16]!" % saved[0], ind))
    elif saved:
        s.append(
            ins(
                "stp",
                "%s, %s, [sp, #-%d]!" % (pairs[0][0], pairs[0][1], save_bytes),
                ind,
            )
        )
        for i, (r0, r1) in enumerate(pairs[1:], 1):
            s.append(ins("stp", "%s, %s, [sp, #%d]" % (r0, r1, 16 * i), ind))
        if odd:
            s.append(ins("str", "%s, [sp, #%d]" % (odd, 16 * len(pairs)), ind))
    return "".join(s)


def frame_close(saved, bufbytes=0, ind="    "):
    pairs, odd = _save_pairs(saved)
    save_bytes = 16 * ((len(saved) + 1) // 2)
    s = []
    if bufbytes:
        buf = (bufbytes + 15) & ~15
        if odd:
            s.append(ins("ldr", "%s, [sp, #%d]" % (odd, buf + 16 * len(pairs)), ind))
        for i in range(len(pairs) - 1, -1, -1):
            s.append(
                ins(
                    "ldp",
                    "%s, %s, [sp, #%d]" % (pairs[i][0], pairs[i][1], buf + 16 * i),
                    ind,
                )
            )
        s.append(ins("add", "sp, sp, #%d" % (buf + save_bytes), ind))
    elif len(saved) == 1:
        s.append(ins("ldr", "%s, [sp], #16" % saved[0], ind))
    elif saved:
        if odd:
            s.append(ins("ldr", "%s, [sp, #%d]" % (odd, 16 * len(pairs)), ind))
        for i in range(len(pairs) - 1, 0, -1):
            s.append(
                ins(
                    "ldp", "%s, %s, [sp, #%d]" % (pairs[i][0], pairs[i][1], 16 * i), ind
                )
            )
        s.append(
            ins(
                "ldp", "%s, %s, [sp], #%d" % (pairs[0][0], pairs[0][1], save_bytes), ind
            )
        )
    return "".join(s)


# --------------------------------------------------------------------------
# Header
# --------------------------------------------------------------------------
def gen_header(P):
    return f"""// Generated by scripts/sat64_generator/gen_fp_asm_arm64.py -- do not edit by hand.
// p = 0x{P.p:x} ({P.bits} bits, {P.N} limbs, {P.spare} spare bits)
//
// AArch64 assembly for GF(p) / GF(p^2) multiplication and squaring, with
// all inputs and outputs fully reduced in [0, p). 64x64->128 products use
// MUL (low) + UMULH (high); each multiply-accumulate is two ADDS/ADCS
// chains (low halves, then high halves one limb up). Montgomery reduction
// uses mu = 1 and the single non-zero limb of p+1 (PP1_TOP =
// {hexq(P.pp1_top)}, rematerialised by MOVZ/MOVK), so each step is one
// 64x64 multiply into the top two accumulator limbs; limb 0 is dropped by
// operand rotation. Constants are materialised inline: no .rodata.
//
// Calling convention (AAPCS64): x0 = out, x1 = a, x2 = b.

// Symbol names go through MIKE's namespacing macros, so CDECL must
// expand its argument before pasting the Mach-O underscore onto it.
#include <mike_namespace.h>

#if defined(__APPLE__)
#define CDECL_(x) _##x
#define CDECL(x) CDECL_(x)
#else
#define CDECL(x) x
#endif

#if defined(__linux__) && defined(__ELF__)
.section .note.GNU-stack,"",@progbits
#endif

.text
.p2align 4
"""


# --------------------------------------------------------------------------
# Macros
# --------------------------------------------------------------------------
def _arglist(*groups):
    return ", ".join(",".join(g) for g in groups)


def gen_macros(P):
    N, W = P.N, P.W
    Z = ["Z%d" % i for i in range(N + 1)]
    A = ["A%d" % i for i in range(N)]
    bz = ["\\" + z for z in Z]
    ba = ["\\" + a for a in A]

    out = [f"""
///////////////////////////////////////////////////////////////// MACROS
// z = a x bi   (initial product, no input accumulator)
// Inputs: a in registers [{A[0]}:{A[N-1]}],
//         bi in register BI
// Output: [{Z[0]}:{Z[N]}]
// Temps:  reg T0
// Notes:  the {N} low halves land directly in [{Z[0]}:{Z[N-1]}], then the
//         high halves are summed one limb up in a single ADDS/ADCS
//         chain. One temp suffices: MUL/UMULH leave the flags untouched,
//         and the renamer gives each product a fresh physical T0.
/////////////////////////////////////////////////////////////////
.macro MUL64x{W} {_arglist(Z, A, ["BI"], ["T0"])}
"""]
    for i in range(N):
        out.append(ins("mul", "%s, %s, \\BI" % (bz[i], ba[i])))
    for i in range(N - 1):
        out.append(ins("umulh", "\\T0, %s, \\BI" % ba[i]))
        out.append(
            ins("adds" if i == 0 else "adcs", "%s, %s, \\T0" % (bz[i + 1], bz[i + 1]))
        )
    out.append(ins("umulh", "%s, %s, \\BI" % (bz[N], ba[N - 1])))
    out.append(ins("adc", "%s, %s, xzr" % (bz[N], bz[N])))
    out.append(f""".endm

/////////////////////////////////////////////////////////////////
// z = a x bi + z
// Inputs: a in registers [{A[0]}:{A[N-1]}],
//         bi in register BI,
//         accumulator z in [{Z[0]}:{Z[N]}]
// Output: [{Z[0]}:{Z[N]}]
// Temps:  reg T0
// Notes:  same single-temp two-chain schedule as MUL64x{W}.
/////////////////////////////////////////////////////////////////
.macro MULADD64x{W} {_arglist(Z, A, ["BI"], ["T0"])}
""")
    for i in range(N):
        out.append(ins("mul", "\\T0, %s, \\BI" % ba[i]))
        out.append(ins("adds" if i == 0 else "adcs", "%s, %s, \\T0" % (bz[i], bz[i])))
    out.append(ins("adc", "%s, %s, xzr" % (bz[N], bz[N])))
    for i in range(N):
        out.append(ins("umulh", "\\T0, %s, \\BI" % ba[i]))
        out.append(
            ins(
                "adds" if i == 0 else ("adc" if i == N - 1 else "adcs"),
                "%s, %s, \\T0" % (bz[i + 1], bz[i + 1]),
            )
        )
    out.append(""".endm

// Montgomery word-reduction step: m = Z0; add m x PP1_TOP into the top two
// limbs (mu = 1). PP1_TOP is rematerialised in T0; callers pass the scalar
// register as T1 (dead during the reduction, reloaded for the next column).
.macro MULADD64x64 %s, T0,T1
""" % _arglist(Z))
    out.append(mov_imm("\\T0", P.pp1_top))
    out.append(ins("mul", "\\T1, %s, \\T0" % bz[0]))
    out.append(ins("umulh", "\\T0, %s, \\T0" % bz[0]))
    out.append(ins("adds", "%s, %s, \\T1" % (bz[N - 1], bz[N - 1])))
    out.append(ins("adc", "%s, %s, \\T0" % (bz[N], bz[N])))
    out.append(".endm\n")
    return "".join(out)


def _window(names, order):
    """'[x10:x17, x9]' accumulator-window display."""
    idx = {r: i for i, r in enumerate(order)}
    segs, i = [], 0
    while i < len(names):
        j = i
        while j + 1 < len(names) and idx[names[j + 1]] == idx[names[j]] + 1:
            j += 1
        segs.append(names[i] if i == j else "%s:%s" % (names[i], names[j]))
        i = j + 1
    return "[%s]" % ", ".join(segs)


def gen_fpmul_macro(P):
    N, W = P.N, P.W
    Z = ["Z%d" % i for i in range(N + 1)]
    A = ["A%d" % i for i in range(N)]
    st = list(Z)
    out = [f"""
///////////////////////////////////////////////////////////////// MACRO
// z = a x b / R (Montgomery, raw: z < a*b/R + p)
// Inputs: scalar source pointer M0 (b; b[0] is already folded into z by
//         the caller's initial product), a in registers [{A[0]}:{A[N-1]}],
//         accumulator z in [{Z[0]}:{Z[N]}] pre-loaded with a x b[0].
// Output: [{Z[0]}:{Z[N]}] (rotated by {N} positions)
// Temps:  BI (scalar; doubles as the reduction's 2nd temp), reg T0
/////////////////////////////////////////////////////////////////
.macro FPMUL{W}x{W} M0, {_arglist(A, Z, ["BI"], ["T0"])}
"""]

    def red(state):
        s = "    // %s <- z = (z0 x p_plus_1 + z)/2^64\n" % _window(state[1:], Z)
        s += "    MULADD64x64 %s, \\T0, \\BI\n" % _arglist(["\\" + z for z in state])
        return s, state[1:] + state[:1]

    code, st = red(st)
    out.append(code)
    for i in range(1, N):
        out.append("\n    // %s <- z = a x b%d + z\n" % (_window(st, Z), i))
        out.append(ins("mov", "\\%s, xzr" % st[N]))
        out.append(ins("ldr", "\\BI, [\\M0, #%d]" % (8 * i)))
        out.append(
            "    MULADD64x%d  %s\n"
            % (
                W,
                _arglist(
                    ["\\" + z for z in st], ["\\" + a for a in A], ["\\BI"], ["\\T0"]
                ),
            )
        )
        code, st = red(st)
        out.append(code)
    out.append(".endm\n")
    return "".join(out)


# --------------------------------------------------------------------------
# Building blocks (function bodies, real registers)
# --------------------------------------------------------------------------
def montmul(P, M0, A, st, bi, t0, first_scalar):
    """st = A x [M0] / R (raw, < 2p for the operand bounds above); returns
    the rotated accumulator, whose first N registers hold the result and
    whose last one is zero. first_scalar loads [M0]'s limb 0 into bi."""
    s = ["    // %s <- z = a x b0\n" % reg_range(st)]
    s.append(ins("ldr", "%s, %s" % (bi, first_scalar)))
    s.append("    MUL64x%d  %s\n" % (P.W, _arglist(st, A, [bi], [t0])))
    s.append("    FPMUL%dx%d  %s, %s\n" % (P.W, P.W, M0, _arglist(A, st, [bi], [t0])))
    return "".join(s), st[P.N :] + st[: P.N]


def reduce_once(P, a, m, t, ind="    "):
    """a -= p, then add p back if that borrowed: maps [0, 2p) to [0, p).
    The low N-1 limbs of p are all-ones, so one register (m) holds them,
    and later the add-back mask."""
    s = ["%s// reduce %s into [0, p)\n" % (ind, reg_range(a))]
    s.append(mov_imm(m, MASK64, ind))
    s.append(mov_imm(t, P.p_top, ind, comment="p[%d]" % (P.N - 1)))
    s.append(ins("subs", "%s, %s, %s" % (a[0], a[0], m), ind))
    for r in a[1:-1]:
        s.append(ins("sbcs", "%s, %s, %s" % (r, r, m), ind))
    s.append(ins("sbcs", "%s, %s, %s" % (a[-1], a[-1], t), ind))
    s.append(_add_p_on_borrow(P, a, m, t, ind))
    return "".join(s)


def _add_p_on_borrow(P, a, m, t, ind="    "):
    """a += p & -borrow, where the borrow is the carry flag left by the
    preceding sbcs chain; t must already hold p's top limb."""
    s = [tail(ins("sbc", "%s, xzr, xzr" % m, ind), "mask = -borrow")]
    s.append(ins("and", "%s, %s, %s" % (t, t, m), ind))
    s.append(ins("adds", "%s, %s, %s" % (a[0], a[0], m), ind))
    for r in a[1:-1]:
        s.append(ins("adcs", "%s, %s, %s" % (r, r, m), ind))
    s.append(ins("adc", "%s, %s, %s" % (a[-1], a[-1], t), ind))
    return "".join(s)


def sub_mod(P, a, b, m, t, ind="    "):
    """a = a - b mod p for a, b in [0, p), both in registers."""
    s = [mov_imm(t, P.p_top, ind, comment="p[%d]" % (P.N - 1))]
    s.append(ins("subs", "%s, %s, %s" % (a[0], a[0], b[0]), ind))
    for x, y in zip(a[1:], b[1:]):
        s.append(ins("sbcs", "%s, %s, %s" % (x, x, y), ind))
    s.append(_add_p_on_borrow(P, a, m, t, ind))
    return "".join(s)


def add_raw(P, dst, base, index, t0, t1, ind="    "):
    """dst += [base + 8*index ..] without reduction (the sum must fit in N
    limbs), the second operand ldp'd pairwise through t0/t1."""
    s, n = [], P.N
    for i in range(0, n - 1, 2):
        s.append(ld_vec([t0, t1], base, index + i, ind))
        s.append(
            ins("adds" if i == 0 else "adcs", "%s, %s, %s" % (dst[i], dst[i], t0), ind)
        )
        op = "adc" if i + 1 == n - 1 else "adcs"
        s.append(ins(op, "%s, %s, %s" % (dst[i + 1], dst[i + 1], t1), ind))
    if n % 2:
        s.append(ld_vec([t0], base, index + n - 1, ind))
        s.append(ins("adc", "%s, %s, %s" % (dst[n - 1], dst[n - 1], t0), ind))
    return "".join(s)


def sub_raw(P, dst, base, index, t0, t1, ind="    "):
    """dst -= [base + 8*index ..], leaving the borrow in the carry flag."""
    s, n = [], P.N
    for i in range(0, n - 1, 2):
        s.append(ld_vec([t0, t1], base, index + i, ind))
        s.append(
            ins("subs" if i == 0 else "sbcs", "%s, %s, %s" % (dst[i], dst[i], t0), ind)
        )
        s.append(ins("sbcs", "%s, %s, %s" % (dst[i + 1], dst[i + 1], t1), ind))
    if n % 2:
        s.append(ld_vec([t0], base, index + n - 1, ind))
        s.append(ins("sbcs", "%s, %s, %s" % (dst[n - 1], dst[n - 1], t0), ind))
    return "".join(s)


# --------------------------------------------------------------------------
# fp_mul / fp_sqr
# --------------------------------------------------------------------------
def gen_fp_mul(P):
    N = P.N
    A, st, bi, t0 = P.A, list(P.Z), P.BI, P.T0
    saved = saved_regs(A + st + [bi, t0])
    s = [f"""
//***********************************************************************
//  Field multiplication in GF(p)
//  Operation: c = a x b mod p
//  Inputs: a stored in [x1], b stored in [x2]
//  Output: c stored in [x0]
//  Register allocation: a (resident) = {reg_range(A)};
//                       acc z0..z{N} = {reg_range(st)},
//                       scalar bi = {bi}, MULADD temp {t0};
//                       the reduction borrows bi (dead there).
//***********************************************************************
.global CDECL(fp_mul)
.p2align 6
CDECL(fp_mul):
"""]
    s.append(frame_open(saved))
    s.append(ld_vec(A, "x1"))
    code, st = montmul(P, "x2", A, st, bi, t0, "[x2]")
    s.append(code)
    s.append(reduce_once(P, st[:N], t0, bi))
    s.append(st_vec(st[:N], "x0"))
    s.append(frame_close(saved))
    s.append("    ret\n")
    return "".join(s)


def gen_fp_sqr(P):
    return """
//***********************************************************************
//  Field squaring in GF(p): fp_mul with b = a
//***********************************************************************
.global CDECL(fp_sqr)
.p2align 6
CDECL(fp_sqr):
    mov   x2, x1
    b     CDECL(fp_mul)
"""


# --------------------------------------------------------------------------
# fp2_mul / fp2_sqr
# --------------------------------------------------------------------------
def _out_ptr(P, scratch, slot, ind="    "):
    """The register holding the output pointer, reloading it into scratch
    when x0 was spilled."""
    if not P.spill_out:
        return "x0", ""
    return scratch, tail(
        ins("ldr", "%s, [sp, #%d]" % (scratch, 8 * slot), ind), "out pointer"
    )


def _spill_out(P, slot, ind="    "):
    if not P.spill_out:
        return ""
    return tail(
        ins("str", "x0, [sp, #%d]" % (8 * slot), ind),
        "spill the out pointer: x0 joins the pool",
    )


def gen_fp2_mul(P):
    """c = a x b in GF(p^2) with three Montgomery multiplications:
    P0 = a0 x b0, P1 = a1 x b1, P2 = (a0+a1) x (b0+b1), each reduced
    c0 = P0 - P1, c1 = P2 - P0 - P1 (mod p)"""
    N = P.N
    A, st, bi, t0 = P.A2, list(P.Z2), P.BI2, P.T02
    saved = saved_regs(A + st + [bi, t0])
    # [sp + 0 .. N): s, then P2; [sp + N .. 2N): P0; [sp + 2N]: out pointer
    slots = 2 * N + (1 if P.spill_out else 0)
    s = [f"""
//***********************************************************************
//  Multiplication in GF(p^2) (Karatsuba, 3 Montgomery multiplications)
//  Operation: c0 = a0 x b0 - a1 x b1 ; c1 = a0 x b1 + a1 x b0
//    as P0 = a0 x b0, P1 = a1 x b1, P2 = (a0+a1) x (b0+b1), each reduced
//    into [0, p), then c0 = P0 - P1 and c1 = P2 - P1 - P0 (mod p).
//    a0+a1 and b0+b1 are kept unreduced (< 2p): P2 < 4p^2/R + p < 2p.
//  Inputs: a = [a0, a1] stored in [x1]
//          b = [b0, b1] stored in [x2]
//  Output: c = [c0, c1] stored in [x0]  (written only after a and b are
//          fully read: c may alias a or b)
//  Register allocation: multiplicand (r/b0/b1 in turn) = {reg_range(A)};
//                       acc z0..z{N} = {reg_range(st)},
//                       scalar bi = {bi}, MULADD temp {t0};
//                       scalars (s/a0/a1) streamed from [sp]/[x1].
//  Stack: [sp+0 .. {8 * N}) holds s = a0+a1, then P2;
//         [sp+{8 * N} .. {16 * N}) holds P0.
//***********************************************************************
.global CDECL(fp2_mul)
.p2align 6
CDECL(fp2_mul):
"""]
    s.append(
        frame_open(saved, 8 * slots, "s/P2 + P0" + (" + out" if P.spill_out else ""))
    )
    s.append(_spill_out(P, 2 * N))

    s.append("\n    // s = a0 + a1 (unreduced, < 2p) -> [sp]\n")
    s.append(ld_vec(st[:N], "x1"))
    s.append(add_raw(P, st[:N], "x1", N, bi, t0))
    s.append(st_vec(st[:N], "sp"))
    s.append("    // r = b0 + b1 (unreduced, < 2p) -> %s\n" % reg_range(A))
    s.append(ld_vec(A, "x2"))
    s.append(add_raw(P, A, "x2", N, bi, t0))

    s.append("\n    // P2 = r x s -> [sp], overwriting s\n")
    code, st = montmul(P, "sp", A, st, bi, t0, "[sp]")
    s.append(code)
    s.append(reduce_once(P, st[:N], bi, t0))
    s.append(st_vec(st[:N], "sp"))

    s.append("\n    // P0 = b0 x a0 -> [sp+%d]\n" % (8 * N))
    s.append(ld_vec(A, "x2"))
    code, st = montmul(P, "x1", A, st, bi, t0, "[x1]")
    s.append(code)
    s.append(reduce_once(P, st[:N], bi, t0))
    s.append(st_vec(st[:N], "sp", N))

    s.append("\n    // P1 = b1 x a1 -> %s\n" % reg_range(st[:N]))
    s.append(ld_vec(A, "x2", N))
    s.append(tail(ins("add", "x1, x1, #%d" % (8 * N)), "a1; a and b are dead after P1"))
    code, st = montmul(P, "x1", A, st, bi, t0, "[x1]")
    s.append(code)
    s.append(reduce_once(P, st[:N], bi, t0))

    out, reload = _out_ptr(P, st[N], 2 * N)
    s.append("\n    // c0 = P0 - P1 -> %s\n" % reg_range(A))
    s.append(ld_vec(A, "sp", N))
    s.append(sub_mod(P, A, st[:N], bi, t0))
    s.append(reload)
    s.append(st_vec(A, out))

    s.append("\n    // c1 = (P2 - P1) - P0\n")
    s.append(ld_vec(A, "sp"))
    s.append(sub_mod(P, A, st[:N], bi, t0))
    s.append(ld_vec(st[:N], "sp", N))
    s.append(sub_mod(P, A, st[:N], bi, t0))
    s.append(st_vec(A, out, N))

    s.append(frame_close(saved, 8 * slots))
    s.append("    ret\n")
    return "".join(s)


def gen_fp2_sqr(P):
    """c = a^2 in GF(p^2) with two Montgomery multiplications:
        c0 = (a0+a1) x (a0-a1+p), c1 = a1 x 2*a0
    c0 is stored before a1 is reloaded for c1, which is safe under exact
    aliasing (c == a): the c0 store only touches a0's slot."""
    N = P.N
    A, st, bi, t0 = P.A2, list(P.Z2), P.BI2, P.T02
    saved = saved_regs(A + st + [bi, t0])
    # [sp + 0 .. N): s; [sp + N .. 2N): t; [sp + 2N]: out pointer
    slots = 2 * N + (1 if P.spill_out else 0)
    s = [f"""
//***********************************************************************
//  Squaring in GF(p^2)
//  Operation: c0 = (a0+a1) x (a0-a1) ; c1 = 2 a0 x a1
//    with a0+a1, a0-a1+p and 2*a0 kept unreduced (< 2p), so both
//    Montgomery products are < 2p before their final reduction.
//  Inputs: a = [a0, a1] stored in [x1]
//  Output: c = [c0, c1] stored in [x0]  (c may alias a)
//  Register allocation: multiplicand (d = a0-a1+p, then a1) = {reg_range(A)};
//                       acc z0..z{N} = {reg_range(st)},
//                       scalar bi = {bi}, MULADD temp {t0};
//                       scalars streamed from [sp].
//  Stack: [sp+0 .. {8 * N}) holds s = a0+a1;
//         [sp+{8 * N} .. {16 * N}) holds t = 2*a0.
//***********************************************************************
.global CDECL(fp2_sqr)
.p2align 6
CDECL(fp2_sqr):
"""]
    s.append(frame_open(saved, 8 * slots, "s + t" + (" + out" if P.spill_out else "")))
    s.append(_spill_out(P, 2 * N))

    s.append(
        "\n    // t = 2*a0 (unreduced) -> [sp+%d] ; a0 stays in %s\n"
        % (8 * N, reg_range(A))
    )
    s.append(ld_vec(A, "x1"))
    s.append(ins("adds", "%s, %s, %s" % (st[0], A[0], A[0])))
    for z, a in zip(st[1 : N - 1], A[1:-1]):
        s.append(ins("adcs", "%s, %s, %s" % (z, a, a)))
    s.append(ins("adc", "%s, %s, %s" % (st[N - 1], A[-1], A[-1])))
    s.append(st_vec(st[:N], "sp", N))

    s.append("    // s = a0 + a1 (unreduced) -> [sp]\n")
    s.append(ld_vec(st[:N], "x1", N))
    s.append(ins("adds", "%s, %s, %s" % (st[0], st[0], A[0])))
    for z, a in zip(st[1 : N - 1], A[1:-1]):
        s.append(ins("adcs", "%s, %s, %s" % (z, z, a)))
    s.append(ins("adc", "%s, %s, %s" % (st[N - 1], st[N - 1], A[-1])))
    s.append(st_vec(st[:N], "sp"))

    s.append("    // d = a0 - a1 + p (unreduced, in (0, 2p)) -> %s\n" % reg_range(A))
    s.append(sub_raw(P, A, "x1", N, st[0], st[1]))
    s.append(mov_imm(st[0], MASK64))
    s.append(mov_imm(st[1], P.p_top, comment="p[%d]" % (N - 1)))
    s.append(ins("adds", "%s, %s, %s" % (A[0], A[0], st[0])))
    for r in A[1:-1]:
        s.append(ins("adcs", "%s, %s, %s" % (r, r, st[0])))
    s.append(ins("adc", "%s, %s, %s" % (A[-1], A[-1], st[1])))

    s.append("\n    // c0 = d x s -> [out] (a0's slot; a1 is still intact if c == a)\n")
    code, st = montmul(P, "sp", A, st, bi, t0, "[sp]")
    s.append(code)
    s.append(reduce_once(P, st[:N], bi, t0))
    out, reload = _out_ptr(P, bi, 2 * N)
    s.append(reload)
    s.append(st_vec(st[:N], out))

    s.append("\n    // c1 = a1 x t -> [out+%d]\n" % (8 * N))
    s.append(ld_vec(A, "x1", N))
    s.append(
        tail(
            ins("add", "x1, sp, #%d" % (8 * N)),
            "x1 is dead after the a1 load: reuse it as the t pointer",
        )
    )
    code, st = montmul(P, "x1", A, st, bi, t0, "[x1]")
    s.append(code)
    s.append(reduce_once(P, st[:N], bi, t0))
    out, reload = _out_ptr(P, bi, 2 * N)
    s.append(reload)
    s.append(st_vec(st[:N], out, N))

    s.append(frame_close(saved, 8 * slots))
    s.append("    ret\n")
    return "".join(s)


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------
def generate(p: int) -> str:
    P = Params(p)
    return "".join(
        [
            gen_header(P),
            gen_macros(P),
            gen_fpmul_macro(P),
            gen_fp_mul(P),
            gen_fp_sqr(P),
            gen_fp2_mul(P),
            gen_fp2_sqr(P),
        ]
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "prime",
        help="prime as decimal or 0x-prefixed integer, or a Python expression like '633 * 2**308 - 1'",
    )
    ap.add_argument("-o", "--output", default="fp_asm.S")
    args = ap.parse_args()

    from gen_fp import parse_int

    Path(args.output).write_text(generate(parse_int(args.prime)))


if __name__ == "__main__":
    main()
