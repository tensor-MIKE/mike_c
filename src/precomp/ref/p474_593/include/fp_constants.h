#if RADIX == 32
#if defined(MIKE_GF_IMPL_SAT64)
#define NWORDS_FIELD 16
#else
#define NWORDS_FIELD 17
#endif
#define NWORDS_ORDER 16
#elif RADIX == 64
#if defined(MIKE_GF_IMPL_SAT64)
#define NWORDS_FIELD 8
#else
#define NWORDS_FIELD 8
#endif
#define NWORDS_ORDER 8
#endif
