#!/usr/bin/env python3
# This file was generated using a LLM which had information of both:
# 1. The pure C implementation written for an fp2 library in C
#    https://github.com/GiacomoPope/fp2_c
# 2. The ASM generation script written for SQISign
#    https://github.com/SQISign/the-sqisign/blob/main/scripts/gen_fp/gen_fp_asm_broadwell.py

"""
Generator of x86-64 Broadwell (MULX/ADCX/ADOX) assembly for fp_mul,
fp_sqr, fp2_mul and fp2_sqr, for Montgomery-friendly primes
p = c * 2^t - 1.

Matches gen_fp.py's fp_t representation exactly (radix 2^64,
N = ceil(bits/64) limbs, R = 2^(64*N) mod p) and is fully canonical:
unlike the SQIsign broadwell generator, which leaves outputs lazily
reduced in [0, 2^bits) and needs 4 spare top-limb bits to stay bounded,
every routine here finishes with one conditional subtraction of p. That
needs the raw result < 2p, which Montgomery reduction guarantees for
fp_mul/fp_sqr/fp2_mul with one spare bit and for fp2_sqr (unreduced
operands < 2p) with two; AsmParams checks this.

Because these primes have p = c*2^t - 1 with t >= 64*(N-1), the low N-1
limbs of p are all 0xFFFFFFFFFFFFFFFF, forcing mu = -p^-1 mod 2^64 = 1
(a row's own low limb *is* the reduction multiplier) and making p+1 have
a single non-zero 64-bit limb (folding in "m * p" collapses from an
N-limb multiply-add to one 64x64 multiply). fp_mul is row-interleaved
CIOS: one schoolbook row (paired MULX/ADCX/ADOX so the two carry chains
run independently) immediately followed by that one-multiply reduction
fold, keeping the accumulator to N+1 limbs throughout. fp2_mul and
fp2_sqr reuse that loop with fp_mul's register plan, one fused pass per
component; fp_sqr has its own (see emit_sqr_routine).

Usage:
    python gen_fp_asm_broadwell.py --prime "633 * 2**308 - 1" -o fp_asm.S
"""

import argparse
from pathlib import Path

MASK64 = (1 << 64) - 1


class AsmParams:
    def __init__(self, p: int):
        self.p = p
        self.bits = p.bit_length()
        self.n = (self.bits + 63) // 64
        if self.n < 2:
            raise ValueError("primes below 2 limbs are not supported")

        self.modulus = [(p >> (64 * i)) & MASK64 for i in range(self.n)]
        if self.modulus[0] != MASK64:
            raise ValueError(
                "this generator requires mu = -p^-1 mod 2^64 == 1, i.e. the "
                "low limb of p must be 0xFFFFFFFFFFFFFFFF (true for every "
                "p = c*2^t - 1 with t >= 64*(n-1)); got low limb "
                "0x%016x" % self.modulus[0]
            )

        # p + 1 = c * 2^t must have exactly one non-zero 64-bit limb.
        p_plus_1 = p + 1
        top_idx = self.n - 1
        if any(((p_plus_1 >> (64 * i)) & MASK64) != 0 for i in range(top_idx)):
            raise ValueError(
                "p + 1 must have a single non-zero 64-bit limb "
                "(the top one): this generator only targets "
                "p = c*2^t - 1 primes with t >= 64*(n-1)"
            )
        p_plus_1_top = p_plus_1 >> (64 * top_idx)
        if p_plus_1_top.bit_length() > 64:
            raise ValueError(
                "p+1's non-zero limb is %d bits: the reduction fold's single "
                "mulx needs it to fit in 64 bits (got 0x%x)"
                % (p_plus_1_top.bit_length(), p_plus_1_top)
            )
        self.top_idx = top_idx
        self.p_plus_1_top = p_plus_1_top
        # fp2_sqr multiplies two unreduced operands < 2p; its single final
        # subtraction needs the Montgomery output (< 4p^2/R + p) below 2p.
        if 4 * p >= 1 << (64 * self.n):
            raise ValueError(
                "fp2_sqr needs p < 2^(64n)/4 (two spare bits in the top limb)"
            )


def hexq(v: int) -> str:
    return "0x%016X" % (v & MASK64)


# --------------------------------------------------------------------------
# Row-interleaved CIOS multiply
# --------------------------------------------------------------------------
def emit_row(window, b_ptr, n, t0, t1, clear_top=True):
    """window[0..n-1] += a_i * b[0..n-1] (rdx = a_i already), carry into
    window[n]. Paired ADCX/ADOX: term j's low half joins the ADCX chain
    into window[j] (the very first term has nothing to chain off, so it
    goes via ADOX instead) and its high half joins the ADOX chain into
    window[j+1]; one final ADC merges the two chains.

    clear_top: window[n] is the freshly-rotated-in slot and gets zeroed
    here (also clearing CF/OF) before use. Pass False to instead
    accumulate onto an existing window[n] -- e.g. calling this twice into
    the same window for fused sum-of-products -- while still clearing
    CF/OF via a scratch register for this call's own j=0 term.
    """
    if not clear_top:
        lines = [f"xor    {t0}, {t0}"]
    elif is_mem(window[n]):
        # xor is what clears CF/OF, so it needs a register; t0 is free
        # until the first mulx overwrites it.
        lines = [f"xor    {t0}, {t0}", f"mov    {window[n]}, {t0}"]
    else:
        lines = [f"xor    {window[n]}, {window[n]}"]
    for j in range(n):
        b_off = "[%s]" % b_ptr if j == 0 else "[%s+%d]" % (b_ptr, 8 * j)
        lines.append(f"mulx   {t0}, {t1}, {b_off}")
        if j == 0:
            lines.extend(emit_adx("adox", window[0], t1))
            lines.extend(emit_adx("adox", window[1], t0))
        else:
            lines.extend(emit_adx("adcx", window[j], t1))
            lines.extend(emit_adx("adox", window[j + 1], t0))
    lines.append(f"adc    {window[n]}, 0")
    return lines


def is_mem(loc: str) -> bool:
    return "[" in loc


def emit_adx(op, dst, src):
    """dst += src (+ CF for adcx / OF for adox). ADCX/ADOX only take a
    register destination, so a memory-resident window slot is added into
    the (about to be dead) mulx temporary instead and stored back; mov
    leaves both flags alone, so neither carry chain notices."""
    if not is_mem(dst):
        return [f"{op}   {dst}, {src}"]
    return [f"{op}   {src}, {dst}", f"mov    {dst}, {src}"]


def emit_reduction_fold(window, n, top_idx, t0, t1):
    """m = window[0] (mu=1); fold m * (p+1)'s one non-zero limb into
    window[n-1..n] (always the top two slots, since top_idx == n-1)."""
    assert top_idx == n - 1
    return [
        f"mov    rdx, {window[0]}",
        f"mulx   {t0}, {t1}, [rip + FP_P_PLUS_1_TOP]",
        f"add    {window[n - 1]}, {t1}",
        f"adc    {window[n]}, {t0}",
    ]


# --------------------------------------------------------------------------
# Final canonicalization: one conditional subtraction of p
# --------------------------------------------------------------------------
def emit_final_reduce(result_regs, n, mask, tmp, frame_base):
    """result_regs holds a value < 2p (standard Montgomery CIOS bound).
    Subtract p; if that borrows, add it back (branch-free), mirroring
    gen_fp.py's own fp_sub/fp_add correction. The n masked-p limbs are
    staged at frame_base first: AND clobbers flags, so every masked limb
    must be ready before the add/adc chain starts -- interleaving
    "mov+and" between chain steps silently eats the next adc's carry.
    """
    lines = []
    for i in range(n):
        op = "sub" if i == 0 else "sbb"
        lines.append(f"{op}    {result_regs[i]}, [rip + FP_MODULUS + {8 * i}]")
    lines.append(
        f"sbb    {mask}, {mask}"
    )  # mask = 0 - CF: all-ones if we must add p back
    for i in range(n):
        lines.append(f"mov    {tmp}, [rip + FP_MODULUS + {8 * i}]")
        lines.append(f"and    {tmp}, {mask}")
        lines.append(f"mov    [{frame_base}+{8 * i}], {tmp}")
    for i in range(n):
        op = "add" if i == 0 else "adc"
        lines.append(f"{op}    {result_regs[i]}, [{frame_base}+{8 * i}]")
    return lines


# --------------------------------------------------------------------------
# Register pools
# --------------------------------------------------------------------------
# fp_mul needs n+3 registers live at once (n+1-limb window + 2 mul temps).
CALLER_SAVED_POOL = ["r8", "r9", "r10", "r11", "rax"]
CALLEE_SAVED_POOL = ["rbx", "r12", "r13", "r14", "r15", "rbp"]

# Past 11 registers, free up pointer registers one at a time: rcx first
# (copy b to a stack buffer once, so its row reads become [rsp+...]
# instead of [rcx+...] and rcx itself becomes another accumulator), then
# rdi (park the output pointer in a stack slot, reload only for the final
# store), then rsi (park a's pointer too; each row reloads it once to
# fetch a[i] into rdx).
SPILLABLE_POINTER_REGS = ["rcx", "rdi", "rsi"]

MAX_LIMBS_REGS = len(CALLER_SAVED_POOL) + len(CALLEE_SAVED_POOL) - 3
MAX_LIMBS_SPILL = MAX_LIMBS_REGS + len(SPILLABLE_POINTER_REGS)

# With every pointer spilled, rdx (mulx's implicit operand) is the only
# general-purpose register left out, so one more limb is supported by
# keeping a single window slot in a stack slot instead of a register.
MAX_LIMBS_MEM = MAX_LIMBS_SPILL + 1


def choose_registers(n: int):
    """Register-resident if n+3 fits in 11 registers, otherwise spill
    pointer registers one at a time (SPILLABLE_POINTER_REGS) until it
    does, and past that keep one window slot in memory. Returns
    (window_regs, t0, t1, to_save, spill, mem_window): when mem_window is
    set, window_regs holds only n registers and the caller appends the
    stack slot as the window's last entry."""
    base_pool = CALLER_SAVED_POOL + CALLEE_SAVED_POOL
    needed = n + 3
    mem_window = False
    if needed <= len(base_pool):
        spill = []
    else:
        extra = needed - len(base_pool)
        if extra > len(SPILLABLE_POINTER_REGS) + 1:
            raise ValueError(
                "p needs %d limbs: this generator tops out at %d limbs "
                "(register-resident up to %d limbs, pointer-spilling up "
                "to %d, one memory-resident window slot up to %d)"
                % (n, MAX_LIMBS_MEM, MAX_LIMBS_REGS, MAX_LIMBS_SPILL, MAX_LIMBS_MEM)
            )
        if extra > len(SPILLABLE_POINTER_REGS):
            mem_window = True
            needed -= 1
            extra -= 1
        spill = SPILLABLE_POINTER_REGS[:extra]

    pool = CALLER_SAVED_POOL + spill + CALLEE_SAVED_POOL
    regs = pool[:needed]
    n_window_regs = n if mem_window else n + 1
    window_regs = regs[:n_window_regs]
    t0, t1 = regs[n_window_regs], regs[n_window_regs + 1]
    to_save = [r for r in CALLEE_SAVED_POOL if r in regs]
    return window_regs, t0, t1, to_save, spill, mem_window


# --------------------------------------------------------------------------
# Dedicated squaring: r = a^2, about n(n+1)/2 + n multiplies against
# fp_mul's n^2 + n.
#
#   1. cross products a_i*a_j (i < j), one row per a_i with fp_mul's dual
#      carry-chain emit_row on a shrinking window: row i touches positions
#      [2i+1, i+n] only, so positions 2i+1 and 2i+2 are final after it.
#      Final positions below n go to a stack buffer and free their
#      register; positions >= n stay in registers.
#   2. one pass doubles the cross products and adds the squares a_k^2:
#      "adcx x, x" is x = 2x + CF (the bit shifted out carries on through
#      CF) while "adox x, a_k^2 half" adds the squares on the OF chain.
#   3. Montgomery reduction. With mu = 1 and p + 1 = c * 2^(64(n-1)), round
#      i takes m_i = t[i], cancels it exactly, and folds m_i * c in at
#      positions i+n-1 and i+n. Only round 0's fold lands below n (at
#      n-1), so m_i = t[i] untouched for 0 < i < n-1, m_{n-1} = t[n-1] +
#      lo(m_0 * c): all n folds can then be added in one pass with two
#      carry chains instead of rippling each one to the top. Every step
#      only adds to the high half, which ends < 2p < 2^(64n), so nothing
#      overflows along the way.
#   4. one conditional subtraction of p.
#
# At most n + 2 registers are live at once besides the pointers (r, a)
# and rdx: past that r is parked on the stack, then a is copied to the
# stack, which covers 12 limbs without a memory-resident window slot.
# --------------------------------------------------------------------------
SQR_POOL_ORDER = [
    "r8",
    "r9",
    "r10",
    "r11",
    "rax",
    "rcx",
    "rbx",
    "r12",
    "r13",
    "r14",
    "r15",
    "rbp",
]


def emit_sqr_routine(P: "AsmParams") -> list:
    n = P.n
    # 12 always-available registers (with rcx: no b pointer here) plus rdi
    # and rsi once spilled; rdx is mulx's.
    spill_rdi = n + 2 > len(SQR_POOL_ORDER)
    spill_rsi = n + 2 > len(SQR_POOL_ORDER) + 1
    if n + 2 > len(SQR_POOL_ORDER) + 2:
        raise ValueError("fp_sqr: %d limbs needs more than 14 registers" % n)
    pool = (
        SQR_POOL_ORDER[:6]
        + (["rdi"] if spill_rdi else [])
        + (["rsi"] if spill_rsi else [])
        + SQR_POOL_ORDER[6:]
    )

    frame_reduce = ((n * 8 + 15) // 16) * 16
    low_off = frame_reduce  # t[0..n-1]
    a_off = low_off + 8 * n  # copy of a, if rsi is spilled
    rdi_slot_off = a_off + (8 * n if spill_rsi else 0)
    frame_used = rdi_slot_off + (8 if spill_rdi else 0)
    frame_size = ((frame_used + 15) // 16) * 16
    a_base = f"rsp+{a_off}" if spill_rsi else "rsi"

    free = list(pool)
    used = set()

    def alloc():
        r = free.pop(0)
        used.add(r)
        return r

    def release(r):
        free.insert(0, r)

    def low(pos):
        return f"qword ptr [rsp+{low_off + 8 * pos}]"

    def a_mem(k):
        return f"[{a_base}+{8 * k}]" if k else f"[{a_base}]"

    hi, lo = alloc(), alloc()
    body = []
    if spill_rdi:
        body.append(f"mov    [rsp+{rdi_slot_off}], rdi  // park the output pointer")
    if spill_rsi:
        body.append("// copy a to the stack so rsi can hold a limb")
        for k in range(n):
            body.append(f"mov    {hi}, [rsi+{8 * k}]")
            body.append(f"mov    [rsp+{a_off + 8 * k}], {hi}")
    body.append("")

    # ---- 1. cross products
    loc = {}
    body.append("// 1. cross products a_i*a_j, i < j")
    for pos in range(1, n):
        loc[pos] = alloc()
        body.append(f"xor    {loc[pos]}, {loc[pos]}")
    for i in range(n - 1):
        top = i + n
        loc[top] = alloc()  # emit_row zeroes it (clear_top)
        row = n - 1 - i
        window = [loc[pos] for pos in range(2 * i + 1, top + 1)]
        body.append(
            f"// row {i}: positions [{2 * i + 1}, {top}] += a[{i}] * a[{i + 1}..{n - 1}]"
        )
        body.append(f"mov    rdx, {a_mem(i)}")
        body.extend(emit_row(window, f"{a_base}+{8 * (i + 1)}", row, hi, lo))
        for pos in (2 * i + 1, 2 * i + 2):
            if pos < n and pos in loc:
                body.append(f"mov    {low(pos)}, {loc[pos]}")
                release(loc.pop(pos))
    body.append("")

    # ---- 2. double + squares
    body.append(
        "// 2. t = 2*t + sum a_k^2 * 2^(128k): CF chain doubles, OF chain adds squares"
    )
    x = alloc()
    body.append(f"xor    {x}, {x}  // clears CF and OF")
    for k in range(n):
        body.append(f"mov    rdx, {a_mem(k)}")
        body.append(f"mulx   {hi}, {lo}, rdx")
        for pos, add in ((2 * k, lo), (2 * k + 1, hi)):
            if pos == 0:
                body.append(
                    f"mov    {low(0)}, {lo}  // t[0] = 2*0 + lo, flags untouched"
                )
            elif pos == 2 * n - 1:
                loc[pos] = x
                body.append(
                    f"mov    {x}, 0  // top limb starts at 0 (mov keeps the flags)"
                )
                body.append(f"adcx   {x}, {x}")
                body.append(f"adox   {x}, {add}")
            elif pos in loc:
                body.append(f"adcx   {loc[pos]}, {loc[pos]}")
                body.append(f"adox   {loc[pos]}, {add}")
            else:
                body.append(f"mov    {x}, {low(pos)}")
                body.append(f"adcx   {x}, {x}")
                body.append(f"adox   {x}, {add}")
                body.append(f"mov    {low(pos)}, {x}")
    body.append("")

    # ---- 3. reduction folds
    body.append("// 3. fold m_i * c in at positions i+n-1, i+n (m_i = t[i])")
    R = [loc[pos] for pos in range(n, 2 * n)]
    body.append(f"xor    {lo}, {lo}  // clears CF and OF")
    body.append(f"mov    rdx, {low(0)}")
    body.append(f"mulx   {hi}, {lo}, [rip + FP_P_PLUS_1_TOP]")
    body.append(f"adcx   {lo}, {low(n - 1)}")
    body.append(f"mov    {low(n - 1)}, {lo}  // m_(n-1) = t[n-1] + lo(m_0 * c)")
    body.append(f"adox   {R[0]}, {hi}")
    for i in range(1, n):
        body.append(f"mov    rdx, {low(i)}")
        body.append(f"mulx   {hi}, {lo}, [rip + FP_P_PLUS_1_TOP]")
        body.append(f"adcx   {R[i - 1]}, {lo}")
        body.append(f"adox   {R[i]}, {hi}")
    body.append(f"mov    {lo}, 0")
    body.append(f"adcx   {R[n - 1]}, {lo}")
    body.append("")

    # ---- 4. canonicalize and store
    body.append("// 4. result is < 2p; subtract p once if needed")
    body.extend(emit_final_reduce(R, n, hi, lo, "rsp"))
    if spill_rdi:
        body.append(f"mov    {hi}, [rsp+{rdi_slot_off}]  // reload the output pointer")
        ptr = hi
    else:
        ptr = "rdi"
    for i in range(n):
        body.append(
            f"mov    [{ptr} + {8 * i}], {R[i]}" if i else f"mov    [{ptr}], {R[i]}"
        )

    to_save = [r for r in CALLEE_SAVED_POOL if r in used]
    lines = [
        ".global CDECL(fp_sqr)",
        "CDECL(fp_sqr):",
        "// void fp_sqr(fp_t *r, const fp_t *a);",
    ]
    lines += [f"push   {r}" for r in to_save]
    lines.append(f"sub    rsp, {frame_size}")
    lines += body
    lines.append(f"add    rsp, {frame_size}")
    lines += [f"pop    {r}" for r in reversed(to_save)]
    lines += ["ret", ""]
    return lines


# --------------------------------------------------------------------------
# Fused GF(p^2) multiplication: r = a * b with
#     r.re = a0*b0 - a1*b1 = a0*b0 + a1*(p - b1)
#     r.im = a0*b1 + a1*b0
# Each half is one fused CIOS pass (two schoolbook multiply-adds per row
# into the same window, then a single reduction fold), as in
# fp_sum_of_products. Unlike that routine, this one takes the fp2 operands
# by pointer (r, a, b), so register pressure is exactly fp_mul's: it
# reuses fp_mul's register plan, including pointer spilling and the
# memory-resident window slot, and covers every limb count fp_mul does.
#
# b0, b1 and p - b1 are copied to the stack up front (b is then never read
# through its pointer again, so rcx is never needed). p - b1 lies in
# [1, p], so both sums are < 2p^2 and, with p < R/2 (at least one spare
# bit), the Montgomery output is < 2p^2/R + p < 2p: one conditional
# subtraction of p leaves it canonical. No lazy-reduction headroom is
# assumed beyond that.
#
# r.re is parked on the stack while r.im is computed, and nothing is
# written to r before both are done, so r may alias a or b.
# --------------------------------------------------------------------------
def emit_fp2_fused_routine(
    P: "AsmParams", name, signature, buffers, prep, passes
) -> list:
    """Shared skeleton for the GF(p^2) routines: r = (re, im), each half one
    fused CIOS pass over operands that prep() stages in n-limb stack
    buffers (the b-side operands always come from the stack, so no b
    pointer -- and no rcx -- is needed after prep).

    buffers: names of the n-limb stack buffers prep() fills.
    prep(lines, t0, buf): emits the staging code; buf[name] is a stack
        offset; rsi = a, rdx = the third argument (if any), both still live.
    passes: [(label, products)] for re then im, products a list of
        (a_src, b_buf): a_src is ("a", byte offset) to read a row
        multiplier through a's pointer, or ("buf", name) from the stack.
    """
    n = P.n
    window_regs, t0, t1, to_save, spill, mem_window = choose_registers(n)
    spill_rdi = "rdi" in spill
    spill_rsi = "rsi" in spill

    frame_reduce = ((n * 8 + 15) // 16) * 16
    buf = {}
    off = frame_reduce
    for b in buffers:
        buf[b] = off
        off += 8 * n
    re_off = off
    rdi_slot_off = re_off + 8 * n
    rsi_slot_off = rdi_slot_off + (8 if spill_rdi else 0)
    mem_slot_off = rsi_slot_off + (8 if spill_rsi else 0)
    frame_used = mem_slot_off + (8 if mem_window else 0)
    frame_size = ((frame_used + 15) // 16) * 16
    if mem_window:
        window_regs = window_regs + [f"qword ptr [rsp+{mem_slot_off}]"]

    lines = [f".global CDECL({name})", f"CDECL({name}):"]
    lines.extend(signature)
    for r in to_save:
        lines.append(f"push   {r}")
    lines.append(f"sub    rsp, {frame_size}")
    lines.append("")
    prep(lines, t0, buf)
    if spill_rdi:
        lines.append(f"mov    [rsp+{rdi_slot_off}], rdi  // park the output pointer")
    if spill_rsi:
        lines.append(
            f"mov    [rsp+{rsi_slot_off}], rsi  // park a's pointer; each row reloads it"
        )
    lines.append("")

    def load_rdx(a_src, i):
        kind, where = a_src
        if kind == "buf":
            return [f"mov    rdx, [rsp+{buf[where] + 8 * i}]"]
        off = where + 8 * i
        if spill_rsi:
            return [f"mov    rdx, [rsp+{rsi_slot_off}]", f"mov    rdx, [rdx + {off}]"]
        return [f"mov    rdx, [rsi + {off}]"]

    def fused_pass(label, products):
        out = [f"// {label}", "// zero the window"]
        for r in window_regs:
            if is_mem(r):
                out.append(f"mov    {r}, {window_regs[0]}")
            else:
                out.append(f"xor    {r}, {r}")
        for i in range(n):
            window = [window_regs[(i + k) % (n + 1)] for k in range(n + 1)]
            out.append(f"// row {i}")
            for k, (a_src, b_buf) in enumerate(products):
                out.extend(load_rdx(a_src, i))
                out.extend(
                    emit_row(window, f"rsp+{buf[b_buf]}", n, t0, t1, clear_top=(k == 0))
                )
            out.extend(emit_reduction_fold(window, n, P.top_idx, t0, t1))
        final_window = [window_regs[(n - 1 + k) % (n + 1)] for k in range(n + 1)]
        result_regs = final_window[1:]
        dropped_reg = final_window[0]
        mem_results = [k for k, r in enumerate(result_regs) if is_mem(r)]
        if mem_results:
            k = mem_results[0]
            out.append(
                f"mov    {dropped_reg}, {result_regs[k]}  // memory-resident result limb into the dead slot"
            )
            result_regs = result_regs[:k] + [dropped_reg] + result_regs[k + 1 :]
        out.append("// canonicalize: result is < 2p; subtract p once if needed")
        out.extend(emit_final_reduce(result_regs, n, t0, t1, "rsp"))
        return out, result_regs

    (re_label, re_products), (im_label, im_products) = passes
    re_lines, re_regs = fused_pass(re_label, re_products)
    lines.extend(re_lines)
    for i in range(n):
        lines.append(f"mov    [rsp+{re_off + 8 * i}], {re_regs[i]}")
    lines.append("")

    im_lines, im_regs = fused_pass(im_label, im_products)
    lines.extend(im_lines)
    lines.append("")

    # t0/t1 are free once the final reduction has staged its masked p.
    if spill_rdi:
        lines.append(f"mov    {t1}, [rsp+{rdi_slot_off}]  // reload the output pointer")
        ptr = t1
    else:
        ptr = "rdi"
    for i in range(n):
        lines.append(f"mov    [{ptr} + {8 * (n + i)}], {im_regs[i]}")
    for i in range(n):
        lines.append(f"mov    {t0}, [rsp+{re_off + 8 * i}]")
        lines.append(f"mov    [{ptr} + {8 * i}], {t0}")
    lines.append("")
    lines.append(f"add    rsp, {frame_size}")
    for r in reversed(to_save):
        lines.append(f"pop    {r}")
    lines.append("ret")
    lines.append("")
    return lines


def emit_fp2_mul_routine(P: "AsmParams") -> list:
    n = P.n

    def prep(lines, t0, buf):
        lines.append("// copy b0, b1 to the stack, and p - b1 next to them")
        for j in range(2 * n):
            lines.append(f"mov    {t0}, [rdx+{8 * j}]")
            lines.append(f"mov    [rsp+{buf['b0'] + 8 * j}], {t0}")
        for j in range(n):
            # mov leaves CF alone, so the sub/sbb chain survives the loads/stores
            op = "sub" if j == 0 else "sbb"
            lines.append(f"mov    {t0}, [rip + FP_MODULUS + {8 * j}]")
            lines.append(f"{op}    {t0}, [rdx+{8 * (n + j)}]")
            lines.append(f"mov    [rsp+{buf['nb1'] + 8 * j}], {t0}")

    return emit_fp2_fused_routine(
        P,
        "fp2_mul",
        [
            "// void fp2_mul(fp2_t *r, const fp2_t *a, const fp2_t *b);",
            "// SysV: rdi = r, rsi = a, rdx = b.",
        ],
        ["b0", "b1", "nb1"],
        prep,
        [
            ("r.re = a0*b0 + a1*(p - b1)", [(("a", 0), "b0"), (("a", 8 * n), "nb1")]),
            ("r.im = a0*b1 + a1*b0", [(("a", 0), "b1"), (("a", 8 * n), "b0")]),
        ],
    )


# --------------------------------------------------------------------------
# GF(p^2) squaring: r = a^2 with
#     r.re = a0^2 - a1^2 = (a0 + a1) * (a0 + p - a1)
#     r.im = 2*a0*a1     = a0 * (2*a1)
# Two plain CIOS passes on operands staged once on the stack, none of
# them reduced: a0 + a1 < 2p, a0 + p - a1 < 2p and 2*a1 < 2p all fit in n
# limbs. The worst product is < 4p^2, so with p < R/4 (two spare bits,
# checked in AsmParams) the Montgomery output is < 4p^2/R + p < 2p and one
# conditional subtraction leaves it canonical; the CIOS window stays
# below 3p, so n+1 limbs are enough. Saves fp2_sqr's three modular
# add/subs, its two fp_mul calls and the second read of a.
# --------------------------------------------------------------------------
def emit_fp2_sqr_routine(P: "AsmParams") -> list:
    n = P.n

    def prep(lines, t0, buf):
        lines.append("// s = a0 + a1")
        for j in range(n):
            op = "add" if j == 0 else "adc"
            lines.append(f"mov    {t0}, [rsi+{8 * j}]")
            lines.append(f"{op}    {t0}, [rsi+{8 * (n + j)}]")
            lines.append(f"mov    [rsp+{buf['s'] + 8 * j}], {t0}")
        lines.append("// d = p - a1, then d += a0")
        for j in range(n):
            op = "sub" if j == 0 else "sbb"
            lines.append(f"mov    {t0}, [rip + FP_MODULUS + {8 * j}]")
            lines.append(f"{op}    {t0}, [rsi+{8 * (n + j)}]")
            lines.append(f"mov    [rsp+{buf['d'] + 8 * j}], {t0}")
        for j in range(n):
            op = "add" if j == 0 else "adc"
            lines.append(f"mov    {t0}, [rsi+{8 * j}]")
            lines.append(f"{op}    [rsp+{buf['d'] + 8 * j}], {t0}")
        lines.append("// t = 2 * a1")
        for j in range(n):
            op = "add" if j == 0 else "adc"
            lines.append(f"mov    {t0}, [rsi+{8 * (n + j)}]")
            lines.append(f"{op}    {t0}, {t0}")
            lines.append(f"mov    [rsp+{buf['t'] + 8 * j}], {t0}")

    return emit_fp2_fused_routine(
        P,
        "fp2_sqr",
        ["// void fp2_sqr(fp2_t *r, const fp2_t *a);", "// SysV: rdi = r, rsi = a."],
        ["s", "d", "t"],
        prep,
        [
            ("r.re = (a0 + a1) * (a0 + p - a1)", [(("buf", "s"), "d")]),
            ("r.im = a0 * (2 * a1)", [(("a", 0), "t")]),
        ],
    )


def generate(p: int) -> str:
    P = AsmParams(p)
    n = P.n
    window_regs, t0, t1, to_save, spill, mem_window = choose_registers(n)
    spill_rcx = "rcx" in spill
    spill_rdi = "rdi" in spill
    spill_rsi = "rsi" in spill

    lines = []
    lines.append(
        "// Generated by scripts/sat64_generator/gen_fp_asm_broadwell.py -- do not edit by hand."
    )
    lines.append("// p = 0x%x (%d bits, %d limbs)" % (p, P.bits, n))
    lines.append("")
    lines.append("// Symbol names go through MIKE's namespacing macros, so CDECL must")
    lines.append("// expand its argument before pasting the Mach-O underscore onto it.")
    lines.append("#include <mike_namespace.h>")
    lines.append("")
    lines.append("#if defined(__APPLE__)")
    lines.append("#define CDECL_(x) _##x")
    lines.append("#define CDECL(x) CDECL_(x)")
    lines.append("#else")
    lines.append("#define CDECL(x) x")
    lines.append("#endif")
    lines.append("")
    lines.append(".intel_syntax noprefix")
    lines.append("")
    lines.append("#if defined(__APPLE__)")
    lines.append(".section __TEXT,__const")
    lines.append("#else")
    lines.append(".section .rodata")
    lines.append("#endif")
    lines.append(".p2align 3")
    lines.append("FP_MODULUS:")
    for i in range(n):
        lines.append(f"    .quad {hexq(P.modulus[i])}")
    lines.append("FP_P_PLUS_1_TOP:")
    lines.append(f"    .quad {hexq(P.p_plus_1_top)}")
    lines.append("")
    lines.append(".text")
    lines.append(".p2align 4")
    lines.append(".global CDECL(fp_mul)")
    lines.append("CDECL(fp_mul):")
    lines.append("// void fp_mul(fp_t *r, const fp_t *a, const fp_t *b);")
    lines.append("// SysV: rdi = r, rsi = a, rdx = b.")
    # Stack layout when spilling: final-reduce staging, then (if rcx was
    # spilled) a buffer holding all of b, then (if rdi was spilled) one
    # slot for the output pointer, then (if rsi was spilled) one for a's
    # pointer, then (if needed) the memory-resident window slot.
    frame_reduce = ((n * 8 + 15) // 16) * 16
    b_buf_off = frame_reduce
    b_buf_size = n * 8 if spill_rcx else 0
    rdi_slot_off = frame_reduce + b_buf_size
    rsi_slot_off = rdi_slot_off + (8 if spill_rdi else 0)
    mem_slot_off = rsi_slot_off + (8 if spill_rsi else 0)
    frame_used = mem_slot_off + (8 if mem_window else 0)
    frame_size = ((frame_used + 15) // 16) * 16
    if mem_window:
        window_regs = window_regs + [f"qword ptr [rsp+{mem_slot_off}]"]

    for r in to_save:
        lines.append(f"push   {r}")
    if frame_size:
        note = ["final-reduce staging"]
        if spill_rcx:
            note.append("b buffer")
        if spill_rdi:
            note.append("saved output pointer")
        if spill_rsi:
            note.append("saved a pointer")
        if mem_window:
            note.append("one window limb")
        lines.append(f"sub    rsp, {frame_size}  // {', '.join(note)}")

    if spill_rcx:
        lines.append(
            f"// p needs {n} limbs: b doesn't fit in its own pointer register here,"
        )
        lines.append(
            "// so copy it to the stack once and read every row from [rsp+...]."
        )
        for j in range(n):
            src = "[rdx]" if j == 0 else f"[rdx+{8 * j}]"
            lines.append(f"mov    {t0}, {src}")
            lines.append(f"mov    [rsp+{b_buf_off + 8 * j}], {t0}")
        b_ptr = f"rsp+{b_buf_off}"
    else:
        lines.append(
            "mov    rcx, rdx  // free up rdx for mulx's implicit operand; rcx = b"
        )
        b_ptr = "rcx"

    if spill_rdi:
        lines.append(
            f"mov    [rsp+{rdi_slot_off}], rdi  // output pointer doesn't fit either; park it"
        )
    if spill_rsi:
        lines.append(
            f"mov    [rsp+{rsi_slot_off}], rsi  // and neither does a's; each row reloads it"
        )

    lines.append("// zero the window before row 0: each row's own \"xor window[n],")
    lines.append('// window[n]" only clears the freshly-rotated-in slot, not the whole')
    lines.append("// window, so row 0 needs this done explicitly first.")
    for r in window_regs:
        if is_mem(r):
            lines.append(f"mov    {r}, {window_regs[0]}")
        else:
            lines.append(f"xor    {r}, {r}")
    lines.append("")
    for i in range(n):
        window = [window_regs[(i + k) % (n + 1)] for k in range(n + 1)]
        lines.append(
            f"// row {i}: window[0..{n}] += a[{i}] * b[0..{n - 1}], then reduce"
        )
        if spill_rsi:
            lines.append(f"mov    rdx, [rsp+{rsi_slot_off}]")
            lines.append(f"mov    rdx, [rdx + {8 * i}]" if i else "mov    rdx, [rdx]")
        else:
            lines.append(f"mov    rdx, [rsi + {8 * i}]" if i else "mov    rdx, [rsi]")
        lines.extend(emit_row(window, b_ptr, n, t0, t1))
        lines.extend(emit_reduction_fold(window, n, P.top_idx, t0, t1))
        lines.append("")
    # window[0] of the last row is stale (CIOS guarantees it's
    # mathematically zero) and simply dropped.
    final_window = [window_regs[(n - 1 + k) % (n + 1)] for k in range(n + 1)]
    result_regs = final_window[1:]
    dropped_reg = final_window[0]
    store_ptr_reg = dropped_reg
    mem_results = [k for k, r in enumerate(result_regs) if is_mem(r)]
    if mem_results:
        # The final sub/sbb chain needs every result limb in a register:
        # bring the memory-resident one into the dropped (dead) slot, and
        # use t1 -- free once the final reduction has staged its masked p
        # -- for the output pointer instead.
        k = mem_results[0]
        lines.append(
            f"mov    {dropped_reg}, {result_regs[k]}  // memory-resident result limb into the dead slot"
        )
        result_regs = result_regs[:k] + [dropped_reg] + result_regs[k + 1 :]
        store_ptr_reg = t1
    lines.append("// canonicalize: result is < 2p; subtract p once if needed")
    lines.extend(emit_final_reduce(result_regs, n, t0, t1, "rsp"))
    lines.append("")
    if spill_rdi:
        store_ptr = store_ptr_reg
        lines.append(
            f"mov    {store_ptr}, [rsp+{rdi_slot_off}]  // reload the output pointer"
        )
    else:
        store_ptr = "rdi"
    for i in range(n):
        lines.append(
            f"mov    [{store_ptr} + {8 * i}], {result_regs[i]}"
            if i
            else f"mov    [{store_ptr}], {result_regs[i]}"
        )
    lines.append("")
    if frame_size:
        lines.append(f"add    rsp, {frame_size}")
    for r in reversed(to_save):
        lines.append(f"pop    {r}")
    lines.append("ret")
    lines.append("")

    lines.extend(emit_sqr_routine(P))

    lines.extend(emit_fp2_mul_routine(P))
    lines.extend(emit_fp2_sqr_routine(P))

    lines.append("#if defined(__linux__) && defined(__ELF__)")
    lines.append('.section .note.GNU-stack,"",%progbits')
    lines.append("#endif")
    lines.append("")
    return "\n".join(lines)


def parse_prime(text: str) -> int:
    text = text.strip()
    try:
        return int(eval(text, {"__builtins__": {}}, {}))
    except Exception:
        return int(text, 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "prime",
        nargs="?",
        help="prime as decimal, 0x-prefixed, or a 'c * 2**t - 1' expression",
    )
    ap.add_argument("--prime", dest="prime_opt", help="same as the positional form")
    ap.add_argument("-o", "--output", default="fp_asm.S")
    args = ap.parse_args()

    text = args.prime_opt or args.prime
    if not text:
        ap.error("provide a prime")
    p = parse_prime(text)

    out = generate(p)
    Path(args.output).write_text(out)
    P = AsmParams(p)
    print("wrote %s: p = %d bits, %d limbs" % (args.output, P.bits, P.n))


if __name__ == "__main__":
    main()
