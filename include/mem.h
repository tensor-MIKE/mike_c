// SPDX-License-Identifier: Apache-2.0

#ifndef MEM_H
#define MEM_H
#include <stddef.h>
#include <mike_namespace.h>

/**
 * Clears and frees allocated memory.
 *
 * @param[out] mem Memory to be cleared and freed.
 * @param size Size of memory to be cleared and freed.
 */
void mike_secure_free(void *mem, size_t size);

/**
 * Clears memory.
 *
 * @param[out] mem Memory to be cleared.
 * @param size Size of memory to be cleared.
 */
void mike_secure_clear(void *mem, size_t size);

#endif
