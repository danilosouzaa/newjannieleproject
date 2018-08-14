/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "vec_char.h"
#include "macros.h"

#define INI_CAP 128

struct _VecChar {
    int size;
    int capacity;
    char *vector;
};

VecChar *VChar_createIniCap( const int iniCap )
{
    VecChar *v;
    ALLOCATE( v, VecChar );

    v->capacity = iniCap;
    v->size = 0;
    ALLOCATE_VECTOR( v->vector, char, v->capacity );

    return v;
}

VecChar *VChar_create()
{
    return VChar_createIniCap(INI_CAP);
}

void VChar_pushBack( VecChar *vec, const char value )
{
    if ( vec->size+1 >= vec->capacity ) {
        vec->capacity *=  2;
        REALLOCATE_VECTOR( vec->vector, char, vec->capacity );
    }

    vec->vector[vec->size] = value;
    vec->size++;
}

char VChar_get( const VecChar *vec, const int pos )
{
#ifdef DEBUG
    if ( (pos<0) || (pos>=vec->size) ) {
        fprintf( stderr, "ERROR: trying to get content from position %d in a vector with size %d\n", pos, vec->size );
        abort();
    }
#endif /* DEBUG */
    return vec->vector[pos];
}

void VChar_clear( VecChar *v )
{
    v->size = 0;
}

char *VChar_getPtr( VecChar *v )
{
    return v->vector;
}

void VChar_set( VecChar *vec, const int pos, const char value )
{
#ifdef DEBUG
    if ( (pos<0) || (pos>=vec->size) ) {
        fprintf( stderr, "ERROR: trying to change content of position %d in a vector with size %d\n", pos, vec->size );
        abort();
    }
#endif /* DEBUG */
    vec->vector[pos] = value;
}

int VChar_find( const VecChar *v, const int value )
{
    int i;
    for ( i=0; (i<v->size); ++i )
        if (v->vector[i]==value)
            return i;
    return INT_MAX;
}

/* returns the number of elements in this vector */
int VChar_size( const VecChar *v )
{
    return v->size;
}

/* returns the current capacity of this vector (i.e. including pre-allocated space) */
int VChar_capacity( const VecChar *v )
{
    return v->capacity;
}

/* creates a vector and fills with some initial value */
VecChar *VChar_createFill( const int size, const char value )
{
    VecChar *v = VChar_createIniCap(size);

    FILL( v->vector, 0, v->size, value );

    v->size = size;

    return v;
}

void VChar_free( VecChar **_v )
{
    VecChar *v = *_v;

    free( v->vector );
    free( v );

    *_v = NULL;
}
