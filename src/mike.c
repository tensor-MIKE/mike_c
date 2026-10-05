#include <ex.h>


int
mike_keypair(unsigned char *pk, unsigned char *sk)
{
    public_key_t public_key;
    secret_key_t secret_key;
    int ret;

    // Sample key pair
    ret = protocols_keygen(&public_key, &secret_key);

    // Encode the keys
    public_key_to_bytes(pk, &public_key);
    secret_key_to_bytes(sk, &secret_key);

    return ret; 
}

int
mike_exchange(unsigned char *shared, const unsigned char *pk, const unsigned char *sk)
{
    public_key_t public_key;
    secret_key_t secret_key;
    shared_t shared_key; 
    int ret;

    // Decode secret and public keys
    secret_key_from_bytes(&secret_key, sk); 
    ret = public_key_from_bytes(&public_key, pk); 

    if(!ret)
        return 0; 

    // compute shared secret
    ret &= protocols_exchange(&shared_key, &public_key, &secret_key); 

    // hash shared secret
    hash_shared_secret(shared, &shared_key);
    
    return ret; 
}
