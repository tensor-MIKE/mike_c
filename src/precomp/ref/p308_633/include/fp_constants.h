#if RADIX == 32
#if defined(MIKE_GF_IMPL_SAT64)
#define NWORDS_FIELD 10
#else
#define NWORDS_FIELD 12
#endif
#define NWORDS_ORDER 10
#elif RADIX == 64
#if defined(MIKE_GF_IMPL_SAT64)
#define NWORDS_FIELD 5
#else
#define NWORDS_FIELD 6
#endif
#define NWORDS_ORDER 5
#endif
