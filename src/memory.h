/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef MEMORY_H_DEFINED
#define MEMORY_H_DEFINED

#include <stdlib.h>
#include <sys/types.h>

void *xmalloc( const size_t size );

void *xcalloc( const size_t elements, const size_t size );

void *xrealloc( void *ptr, const size_t size );

#endif

