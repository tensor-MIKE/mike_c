#if RADIX == 32
#if defined(MIKE_GF_IMPL_SAT64)
#define NWORDS_FIELD 18
#else
#define NWORDS_FIELD 20
#endif
#define NWORDS_ORDER 18
#elif RADIX == 64
#if defined(MIKE_GF_IMPL_SAT64)
#define NWORDS_FIELD 9
#else
#define NWORDS_FIELD 10
#endif
#define NWORDS_ORDER 9
#endif
