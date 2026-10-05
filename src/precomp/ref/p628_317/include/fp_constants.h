#if RADIX == 32
#if defined(MIKE_GF_IMPL_SAT64)
#define NWORDS_FIELD 20
#else
#define NWORDS_FIELD 23
#endif
#define NWORDS_ORDER 20
#elif RADIX == 64
#if defined(MIKE_GF_IMPL_SAT64)
#define NWORDS_FIELD 10
#else
#define NWORDS_FIELD 11
#endif
#define NWORDS_ORDER 10
#endif
