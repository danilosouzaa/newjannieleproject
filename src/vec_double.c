/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "vec_double.h"
#include "macros.h"

#define INI_CAP 128

struct _VecDbl {
    int size;
    int capacity;
    double *vector;
};

VecDbl *VDbl_createIniCap( const int iniCap )
{
    VecDbl *v;
    ALLOCATE( v, VecDbl );

    v->capacity = iniCap;
    v->size = 0;
    ALLOCATE_VECTOR( v->vector, double, v->capacity );

    return v;
}

VecDbl *VDbl_create()
{
    return VDbl_createIniCap(INI_CAP);
}

void VDbl_pushBack( VecDbl *vec, double value )
{
    if ( vec->size+1 >= vec->capacity ) {
        vec->capacity *=  2;
        REALLOCATE_VECTOR( vec->vector, double, vec->capacity );
    }

    vec->vector[vec->size] = value;
    vec->size++;
}

double VDbl_get( const VecDbl *vec, const int pos )
{
#ifdef DEBUG
    if ( (pos<0) || (pos>=vec->size) ) {
        fprintf( stderr, "ERROR: trying to get content from position %d in a vector with size %d\n", pos, vec->size );
        abort();
    }
#endif /* DEBUG */
    return vec->vector[pos];
}

void VDbl_clear( VecDbl *v )
{
    v->size = 0;
}

double *VDbl_getPtr( VecDbl *v )
{
    return v->vector;
}

void VDbl_set( VecDbl *vec, const int pos, const double value )
{
#ifdef DEBUG
    if ( (pos<0) || (pos>=vec->size) ) {
        fprintf( stderr, "ERROR: trying to change content of position %d in a vector with size %d\n", pos, vec->size );
        abort();
    }
#endif /* DEBUG */
    vec->vector[pos] = value;
}

int VDbl_find( const VecDbl *v, const double value )
{
    int i;
    for ( i=0; (i<v->size); ++i )
        if (v->vector[i]==value)
            return i;
    return INT_MAX;
}

/* returns the number of elements in this vector */
int VDbl_size( const VecDbl *v )
{
    return v->size;
}

/* returns the current capacity of this vector (i.e. including pre-allocated space) */
int VDbl_capacity( const VecDbl *v )
{
    return v->capacity;
}

/* creates a vector and fills with some initial value */
VecDbl *VDbl_createFill( const int size, const double value )
{
    VecDbl *v = VDbl_createIniCap(size);
    int i;
    for ( i=0; (i<v->size); ++i )
        v->vector[i] = value;
    v->size = size;
    return v;
}

void VDbl_free( VecDbl **_v )
{
    VecDbl *v = *_v;

    free( v->vector );
    free( v );

    *_v = NULL;
}
