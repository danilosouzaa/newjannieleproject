
/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef MACROS_H_INCLUDED
#define MACROS_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef long int Cost;

#define  ROUND( x ) ( floor(x+0.5) )
#define  ROUNDUP( x ) ( ceil(x+0.5) )

#define INT_MAX_M +32767

#define INT_RANDOM( n ) \
   ((int)( ((double)n) * (((double)rand())/(((double)RAND_MAX)+((double)1.0))) ))

#define INT_RANDOM_LU( l,u ) \
   ((int)( l + ((double)u-l+1) * (double)rand()/(((double)RAND_MAX)+((double)1.0))))

#define DOUBLE_RANDOM( l,u ) \
   ( l+ ( (u-l) * (((double)rand())/((double)RAND_MAX + 1)) ))

#define MAX( a, b ) ( (a)>(b) ? (a) : (b) )
#define MIN( a, b ) ( (a)<(b) ? (a) : (b) )

#define STR_SIZE       256
#define FILE_NAME_SIZE 1024
#define LINE_SIZE      2048

#define ALLOCATE( ptr, type ) {\
    ptr = malloc( sizeof(type) ); \
    if (!ptr) { \
       fprintf( stderr, "ERROR: no more memory available. at: %s:%d\n", __FILE__, __LINE__ ); \
       abort(); \
    }; }

#define ALLOCATE_INI( ptr, type ) {\
    ptr = calloc( 1, sizeof(type) ); \
    if (!ptr) { \
       fprintf( stderr, "ERROR: no more memory available. at: %s:%d\n", __FILE__, __LINE__ ); \
       abort(); \
    }; }

#define ALLOCATE_VECTOR( ptr, type, nElements ) {\
    ptr = malloc( sizeof(type)*(nElements) ); \
    if (!ptr) { \
       fprintf( stderr, "ERROR: no more memory available. at: %s:%d\n", __FILE__, __LINE__ ); \
       abort(); \
    }; }

/* allocate filling with zeros */
#define ALLOCATE_VECTOR_INI( ptr, type, nElements ) {\
    ptr = calloc( nElements, sizeof(type) ); \
    if (!ptr) { \
       fprintf( stderr, "ERROR: no more memory available. at: %s:%d\n", __FILE__, __LINE__ ); \
       abort(); \
    }; }

#define REALLOCATE_VECTOR( ptr, type, newSize ) { \
    type *newV = (type*) realloc( ptr, sizeof(type)*(newSize) ); \
    if (!ptr) { \
       fprintf( stderr, "ERROR: no more memory available. at: %s:%d\n", __FILE__, __LINE__ ); \
       abort(); \
    } \
    ptr = newV; \
}

#define OPEN_FILE( f, fileName, mode )  \
    f = fopen( fileName, mode ); \
    if (!f) { \
        fflush(stdout); \
        fprintf( stderr, "ERROR: could not open file %s with mode %s. At %s:%d\n", fileName, mode, __FILE__, __LINE__ ); \
        fflush(stderr);  abort(); exit(EXIT_FAILURE); }

#define True  1
#define False 0

/* fills from start until the last element before end */
#define FILL( vector, start, end, value ) { \
    int i; \
    for ( i=start ; (i<end) ; ++i ) vector[i] = value; \
} \

/* copies elements of a vector of a given type */
#define COPY_VECTOR( target, source, type, nElements ) memcpy( target, source, sizeof(type)*(nElements) );

/* combine elements of two vector of a given type */
#define COMBINE_VECTOR( type, target, sizet, source, sizes ) \
    REALLOCATE_VECTOR(target,type,(sizet + sizes)); \
    if (target != NULL)\
        mempcpy (mempcpy (target, target, sizet), source, sizes);\



#define CLEAR_VECTOR( vector, type, nElements ) \
 memset( vector, 0, sizeof( type )*(nElements) );

#define EPS 1e-5

#define E 2.71828

#define DBL_EQUAL( v1, v2 ) ( fabs(v1-v2)<=EPS )

#endif

