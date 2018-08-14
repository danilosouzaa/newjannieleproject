/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <stdio.h>
#include "memory.h"

void *xmalloc( const size_t size )
{
    void *result = malloc( size );
    if (!result) {
        fprintf(stderr, "No more memory available. Trying to allocate %zu bytes.", size);
        exit(1);
    }

    return result;
}

void *xcalloc( const size_t elements, const size_t size )
{
    void *result = calloc( elements, size );
    if (!result) {
        fprintf(stderr, "No more memory available. Trying to calloc %zu bytes.", size);
        exit(1);
    }

    return result;
}

void *xrealloc( void *ptr, const size_t size )
{
    void *result = realloc( ptr, size );
    if (!result) {
        fprintf(stderr, "No more memory available. Trying to allocate %zu bytes.", size);
        exit(1);
    }
    return result;
}

