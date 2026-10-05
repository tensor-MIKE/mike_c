/** @file
 *
 * @brief The key generation and exchange protocols
 */

#ifndef NIKE_H
#define NIKE_H

#include <mike_namespace.h>

#include <constants.h>
#include <encoded_sizes.h>
#include <fips202.h>

#include <weil.h>
#include <ec.h>

#include "assert.h"

/** @defgroup Key Key generation and Key agreement protocols
 * @{
 */
/** @defgroup Type Types for MIKE key generation and signature protocols
 * @{
 */

/** @brief Type for the secret keys
 *
 * @typedef secret_key_t
 *
 * @struct secret_key
 *
 */
typedef struct secret_key
{
    digit_t x[SECRETKEY_WORDS];
} secret_key_t;

/** @brief Type for the public keys
 *
 * @typedef public_key_t
 *
 * @struct public_key
 *
 */
typedef struct public_key
{
    fp2_t curveA; /// the public curve A coefficient
} public_key_t;

/** @brief Type for the shared secrets
 *
 * @typedef shared_t
 *
 * @struct shared
 *
 */
typedef struct shared
{
    mike_abs_invariants_t invariant; // the absolute invariant of the scheme
} shared_t;

/** @}
 */

/*************************** Functions *****************************/

/**
 * @brief Key generation
 *
 * @param pk Output: will contain the public key
 * @param sk Output: will contain the secret key
 * @returns 1 if success, 0 otherwise
 */
int protocols_keygen(public_key_t *pk, secret_key_t *sk);

/**
 * @brief Signature computation
 *
 * @param out Output: will contain the shared secret
 * @param sk secret key
 * @param pk public key
 * @returns 1 if success, 0 otherwise
 */
int protocols_exchange(shared_t *out, const public_key_t *pk, const secret_key_t *sk);

/*************************** Encoding *****************************/

/** @defgroup encoding Encoding and decoding functions
 * @{
 */

/**
 * @brief Encodes a secret key as a byte array
 *
 * @param enc : Byte array to encode the secret key (including public key) in
 * @param sk : Secret key to encode
 */
void secret_key_to_bytes(unsigned char *enc, const secret_key_t *sk);

/**
 * @brief Encodes a public key as a byte array
 *
 * @param enc : Byte array to encode only a public key in
 * @param pk : Public key to encode
 */
void public_key_to_bytes(unsigned char *enc, const public_key_t *pk);

/**
 * @brief Decodes a secret key (and public key) from a byte array
 *
 * @param sk : Structure to decode the secret key in
 * @param enc : Byte array to decode
 */
void secret_key_from_bytes(secret_key_t *sk, const unsigned char *enc);

/**
 * @brief Decodes a public key from a byte array
 *
 * @param pk : Structure to decode the public key in
 * @param enc : Byte array to decode
 * 
 * @return -1 if valid, 0 if not valid representation
 */
int public_key_from_bytes(public_key_t *pk, const unsigned char *enc);


/**
 * @brief Hash the shared_key to a byte array
 *
 * @param out : Byte array of the 32B hash of the shared secret 
 * @param shared_secret : shared secret
 */
void hash_shared_secret(uint8_t *out, const shared_t  *shared_secret);


/** @}
 */


/** @}
 */

#endif
