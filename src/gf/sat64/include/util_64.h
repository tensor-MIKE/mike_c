#ifndef UTIL_64_H
#define UTIL_64_H

#include <stdint.h>

#ifdef __has_builtin
#define UTIL_64_HAS_BUILTIN(x) __has_builtin(x)
#else
#define UTIL_64_HAS_BUILTIN(x) 0
#endif

/*
 * x86-64 GCC/Clang/MSVC:
 *
 * Use the native carry/borrow intrinsics. These expose the CPU's carry
 * flag directly instead of expressing the operation through __int128.
 */
#if ((defined(__x86_64__) || defined(_M_X64)) &&                               \
     (defined(__GNUC__) || defined(__clang__) || defined(_MSC_VER)))

#include <immintrin.h>

static inline unsigned char addcarry_u64(unsigned char c, uint64_t a,
                                         uint64_t b, uint64_t *d) {
  return _addcarry_u64(c, (unsigned long long)a, (unsigned long long)b,
                       (unsigned long long *)(void *)d);
}

static inline unsigned char subborrow_u64(unsigned char c, uint64_t a,
                                          uint64_t b, uint64_t *d) {
  return _subborrow_u64(c, (unsigned long long)a, (unsigned long long)b,
                        (unsigned long long *)(void *)d);
}

#elif UTIL_64_HAS_BUILTIN(__builtin_addcll) &&                                  \
    UTIL_64_HAS_BUILTIN(__builtin_subcll)

/*
 * Clang and GCC >= 14 (e.g. AArch64):
 *
 * The add/sub-with-carry builtins compile to a single ADCS/SBCS chain,
 * where the __int128 fallback below rebuilds the carry with CSET/CINC on
 * every limb.
 */
static inline unsigned char addcarry_u64(unsigned char c, uint64_t a,
                                         uint64_t b, uint64_t *d) {
  unsigned long long carry;

  *d = __builtin_addcll(a, b, c, &carry);
  return (unsigned char)carry;
}

static inline unsigned char subborrow_u64(unsigned char c, uint64_t a,
                                          uint64_t b, uint64_t *d) {
  unsigned long long borrow;

  *d = __builtin_subcll(a, b, c, &borrow);
  return (unsigned char)borrow;
}

#else

/*
 * Portable fallback.
 */
static inline unsigned char addcarry_u64(unsigned char c, uint64_t a,
                                         uint64_t b, uint64_t *d) {
  unsigned __int128 t =
      (unsigned __int128)a + (unsigned __int128)b + (unsigned __int128)c;

  *d = (uint64_t)t;
  return (unsigned char)(t >> 64);
}

static inline unsigned char subborrow_u64(unsigned char c, uint64_t a,
                                          uint64_t b, uint64_t *d) {
  unsigned __int128 t =
      (unsigned __int128)a - (unsigned __int128)b - (unsigned __int128)c;

  *d = (uint64_t)t;
  return (unsigned char)(t >> 127);
}

#endif

static inline void mul_u64(uint64_t *lo, uint64_t *hi, uint64_t x, uint64_t y) {
  unsigned __int128 t = (unsigned __int128)x * (unsigned __int128)y;

  *lo = (uint64_t)t;
  *hi = (uint64_t)(t >> 64);
}

static inline void mul_add(uint64_t *lo, uint64_t *hi, uint64_t x, uint64_t y,
                           uint64_t z) {
  unsigned __int128 t =
      (unsigned __int128)x * (unsigned __int128)y + (unsigned __int128)z;

  *lo = (uint64_t)t;
  *hi = (uint64_t)(t >> 64);
}

static inline void mul_add2(uint64_t *lo, uint64_t *hi, uint64_t x, uint64_t y,
                            uint64_t z1, uint64_t z2) {
  unsigned __int128 t = (unsigned __int128)x * (unsigned __int128)y +
                        (unsigned __int128)z1 + (unsigned __int128)z2;

  *lo = (uint64_t)t;
  *hi = (uint64_t)(t >> 64);
}

static inline void mul_x2(uint64_t *lo, uint64_t *hi, uint64_t x1, uint64_t y1,
                          uint64_t x2, uint64_t y2) {
  uint64_t lo1, hi1;
  uint64_t lo2, hi2;

  mul_u64(&lo1, &hi1, x1, y1);
  mul_u64(&lo2, &hi2, x2, y2);

  unsigned char c = addcarry_u64(0, lo1, lo2, &lo1);

  (void)addcarry_u64(c, hi1, hi2, &hi1);

  *lo = lo1;
  *hi = hi1;
}

static inline void mul_x2_add(uint64_t *lo, uint64_t *hi, uint64_t x1,
                              uint64_t y1, uint64_t x2, uint64_t y2,
                              uint64_t z) {
  uint64_t tlo, thi;

  mul_x2(&tlo, &thi, x1, y1, x2, y2);

  unsigned char c = addcarry_u64(0, tlo, z, &tlo);

  (void)addcarry_u64(c, thi, 0, &thi);

  *lo = tlo;
  *hi = thi;
}

static inline uint64_t sign_word(uint64_t x) {
  return (uint64_t)(((int64_t)x) >> 63);
}

#endif /* UTIL_64_H */
