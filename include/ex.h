// SPDX-License-Identifier: Apache-2.0

#ifndef MIKE_H
#define MIKE_H

#include <stdint.h>
#include <mike_namespace.h>
#include <nike.h>

/**
 * MIKE keypair generation.
 *
 * The caller is responsible to allocate sufficient memory to hold pk and sk.
 *
 * @param[out] pk MIKE public key
 * @param[out] sk MIKE secret key
 * @return int status code
 */
MIKE_API 
int mike_keypair(unsigned char *pk, unsigned char *sk);

/**
 * MIKE key exchange.
 *
 * The caller is responsible to allocate sufficient memory to hold shared.
 *
 * @param[out] shared Shared secret
 * @param[in] pk MIKE public key
 * @param[in] sk MIKE secret key
 * @return int status code
 */
MIKE_API 
int mike_exchange(unsigned char *shared,
                 const unsigned char *pk,
                 const unsigned char *sk);



#endif
