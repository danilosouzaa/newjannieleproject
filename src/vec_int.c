/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "vec_int.h"
#include "macros.h"

#define INI_CAP 128

struct _VecInt {
    int size;
    int capacity;
    int *vector;
};

VecInt *VInt_createIniCap( const int iniCap )
{
    VecInt *v;
    ALLOCATE( v, VecInt );

    v->capacity = iniCap;
    v->size = 0;
    ALLOCATE_VECTOR( v->vector, int, v->capacity );

    return v;
}

VecInt *VInt_create()
{
    return VInt_createIniCap(INI_CAP);
}

void VInt_pushBack( VecInt *vec, int value )
{
    if ( vec->size+1 >= vec->capacity ) {
        vec->capacity *=  2;
        REALLOCATE_VECTOR( vec->vector, int, vec->capacity );
    }

    vec->vector[vec->size] = value;
    vec->size++;
}


int VInt_get( const VecInt *vec, const int pos )
{
#ifdef DEBUG
    if ( (pos<0) || (pos>=vec->size) ) {
        fprintf( stderr, "ERROR: trying to get content from position %d in a vector with size %d\n", pos, vec->size );
        abort();
    }
#endif /* DEBUG */
    return vec->vector[pos];
}

void VInt_clear( VecInt *v )
{
    v->size = 0;
}

int *VInt_getPtr( VecInt *v )
{
    return v->vector;
}

void VInt_set( VecInt *vec, const int pos, const int value )
{
#ifdef DEBUG
    if ( (pos<0) || (pos>=vec->size) ) {
        fprintf( stderr, "ERROR: trying to change content of position %d in a vector with size %d\n", pos, vec->size );
        abort();
    }
#endif /* DEBUG */
    vec->vector[pos] = value;
}

int VInt_find( const VecInt *v, const int value )
{
    int i;
    for ( i=0; (i<v->size); ++i )
        if (v->vector[i]==value)
            return i;
    return INT_MAX;
}

/* returns the number of elements in this vector */
int VInt_size( const VecInt *v )
{
    return v->size;
}

/* returns the current capacity of this vector (i.e. including pre-allocated space) */
int VInt_capacity( const VecInt *v )
{
    return v->capacity;
}

/* creates a vector and fills with some initial value */
VecInt *VInt_createFill( const int size, const int value )
{
    VecInt *v = VInt_createIniCap(size);
    int i;
    for ( i=0; (i<v->size); ++i )
        v->vector[i] = value;
    v->size = size;
    return v;
}

void VInt_free( VecInt **_v )
{
    VecInt *v = *_v;

    free( v->vector );
    free( v );

    *_v = NULL;
}
