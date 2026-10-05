#if RADIX == 32
#if defined(MIKE_GF_IMPL_SAT64)
#define NWORDS_FIELD 12
#else
#define NWORDS_FIELD 14
#endif
#define NWORDS_ORDER 12
#elif RADIX == 64
#if defined(MIKE_GF_IMPL_SAT64)
#define NWORDS_FIELD 6
#else
#define NWORDS_FIELD 7
#endif
#define NWORDS_ORDER 6
#endif
