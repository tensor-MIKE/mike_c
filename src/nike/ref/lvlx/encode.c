
#include <nike.h>

void
secret_key_to_bytes(unsigned char *enc, const secret_key_t *sk)
{
    digit_t mask = -1;
    digit_t val;

    for (int i = 0; i < SECRETKEY_WORDS - 1; i++) {
        for (int k = 0; k < RADIX / 8; k++) {
            val = sk->x[i] & (mask << (8 * k));
            val = val >> (8 * k);
            (enc)[(RADIX / 8) * i + k] = (uint8_t)val;
        }
    }
    for (int k = 0; k < SECRETKEY_BYTES % (RADIX / 8); k++) {
        val = sk->x[SECRETKEY_WORDS - 1] & (mask << (8 * k));
        val = val >> (8 * k);
        (enc)[(RADIX / 8) * (SECRETKEY_WORDS - 1) + k] = (uint8_t)val;
    }
}

void
public_key_to_bytes(unsigned char *enc, const public_key_t *pk)
{
    fp2_encode(enc, &pk->curveA);
}

/**
 * @brief Decodes a secret key (and public key) from a byte array
 *
 * @param sk : Structure to decode the secret key in
 * @param enc : Byte array to decode
 */
void
secret_key_from_bytes(secret_key_t *sk, const unsigned char *enc)
{
    digit_t tmp;
    for (int i = 0; i < SECRETKEY_WORDS; i++) {
        sk->x[i] = 0;
        if (i < SECRETKEY_WORDS - 1) {
            for (int k = 0; k < RADIX / 8; k++) {
                tmp = (enc)[(RADIX / 8) * i + k];
                sk->x[i] = sk->x[i] | (tmp << (8 * k));
            }
        } else {
            for (int k = 0; k < SECRETKEY_BYTES % (RADIX / 8); k++) {
                tmp = (enc)[(RADIX / 8) * i + k];
                sk->x[i] = sk->x[i] | (tmp << (8 * k));
            }
        }
    }
}

int 
public_key_from_bytes(public_key_t *pk, const unsigned char *enc)
{
   return fp2_decode(&pk->curveA, enc);
}

void hash_shared_secret(uint8_t *out, const shared_t  *shared_secret)
{
    uint8_t buf[FP_ENCODED_BYTES];

    sha3_256incctx ctx;
    sha3_256_inc_init(&ctx);

    for (int i = 0; i < 4; i++) {
        fp_encode(buf, &shared_secret->invariant[i]);
        sha3_256_inc_absorb(&ctx, buf, sizeof(buf));
    }

    sha3_256_inc_finalize(out, &ctx);
}