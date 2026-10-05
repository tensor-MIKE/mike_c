#include <ex.h>
#include <encoded_sizes.h>
#include <stdio.h>
#include <rng.h>
#include <string.h>

// assumes ascii
static unsigned char
cli_from_hex(char ipt)
{
    unsigned char x;
    if (ipt < 61)
        x = ipt - 48;
    else {
        if (ipt < 97) {
            x = ipt - 55;
        } else {
            x = ipt - 87;
        }
    }
    assert(x < 16);
    return (x);
}

static unsigned char
cli_from_hex_byte(const char *ipt)
{
    unsigned char first, second;
    first = cli_from_hex(ipt[0]);
    second = cli_from_hex(ipt[1]);
    if ((first > 16) || (second > 16))
        return 0;
    return (first * 16 + second);
}

// issue hex encoding of bytes sometimes only 1 char

int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int help = 0;
    int keygen = 0;
    unsigned char *pk;
    unsigned char *sk;
    unsigned char *shared;

    for (int i = 1; i < argc; i++) {
        if (!help && strcmp(argv[i], "--help") == 0) {
            help = 1;
            continue;
        }

        if (!keygen && strcmp(argv[i], "--keygen") == 0) {
            keygen = 1;
            continue;
        }
    }

    if (argc < 1 || help || (!keygen && (argc != 3))) {
        printf("Keygen:   %s [--keygen]\n", argv[0]);
        printf("Exchange: %s <sk> <pk>\n", argv[0]);
        return 1;
    }

    randombytes_select((unsigned char *)seed, sizeof(seed));

    if (keygen) {
        pk = calloc(PUBLICKEY_BYTES, sizeof(char));
        sk = calloc(SECRETKEY_BYTES, sizeof(char));
        mike_keypair(pk, sk);
        printf("MIKE sk 0x");
        for (int i = 0; i < SECRETKEY_BYTES; i++)
            printf("%hhx%hhx", (unsigned char)(sk[i] / 16), sk[i] & 15);
        printf("\n");
        printf("MIKE pk 0x");
        for (int i = 0; i < PUBLICKEY_BYTES; i++)
            printf("%hhx%hhx", (unsigned char)(pk[i] / 16), pk[i] & 15);
        printf("\n");
        free(pk);
        free(sk);
        return (0);
    }
    int sklen = strlen(argv[1]) / 2;
    int pklen = strlen(argv[2]) / 2;
    sk = calloc(sklen, sizeof(char));
    pk = calloc(pklen, sizeof(char));
    shared = calloc(SHARED_BYTES, sizeof(char));
    for (int i = 1; i < pklen; i++)
        pk[i - 1] = cli_from_hex_byte(&argv[2][2 * i]);
    for (int i = 1; i < sklen; i++)
        sk[i - 1] = cli_from_hex_byte(&argv[1][2 * i]);
    sk[sklen - 1] = 0;
    pk[pklen - 1] = 0;
    int ret = mike_exchange(shared, pk, sk);
    printf("MIKE sk 0x");
    for (int i = 0; i < SECRETKEY_BYTES; i++)
        printf("%hhx%hhx", (unsigned char)(sk[i] / 16), sk[i] & 15);
    printf("\n");
    printf("MIKE pk 0x");
    for (int i = 0; i < PUBLICKEY_BYTES; i++)
        printf("%hhx%hhx", (unsigned char)(pk[i] / 16), pk[i] & 15);
    printf("\n");
    printf("MIKE shared 0x");
    for (int i = 0; i < SHARED_BYTES; i++)
        printf("%hhx%hhx", (unsigned char)(shared[i] / 16), shared[i] & 15);
    printf("\n");
    free(shared);
    free(sk);
    free(pk);
    return (ret);
}
