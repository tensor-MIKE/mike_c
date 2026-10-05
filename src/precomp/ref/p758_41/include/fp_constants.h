#if RADIX == 32
#if defined(MIKE_GF_IMPL_SAT64)
#define NWORDS_FIELD 24
#else
#define NWORDS_FIELD 27
#endif
#define NWORDS_ORDER 24
#elif RADIX == 64
#if defined(MIKE_GF_IMPL_SAT64)
#define NWORDS_FIELD 12
#else
#define NWORDS_FIELD 13
#endif
#define NWORDS_ORDER 12
#endif
