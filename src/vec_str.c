
/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "vec_str.h"
#include "memory.h"
#include "str_utils.h"

#define DEF_INI_CAP      256
#define LINE_SIZE       1024
#define FILE_NAME_SIZE  1024

struct _VecStr {
    int capacity;
    int size;

    char **sv;
    char *s;

    int strSize;
};

#ifdef DEBUG
void VStr_check_valid_strv( const VecStr *strv );
#endif

void VStr_increase_capacity_to( VecStr *strv, int newCapacity );

VecStr *VStr_create( int strSize )
{
#ifdef DEBUG
    assert( strSize > 0 );
#endif
    VecStr *res = xmalloc( sizeof(VecStr) );

    res->capacity = DEF_INI_CAP;
    res->size = 0;
    res->strSize = strSize;
    res->sv = xmalloc( sizeof(char*)*res->capacity );

    res->s = xcalloc( res->capacity*strSize, sizeof(char) );

    res->sv[0] = res->s;
    int i;
    for ( i=1; (i<res->capacity); i++ )
        res->sv[i] = res->sv[i-1]+strSize;

#ifdef DEBUG
    VStr_check_valid_strv( res );
#endif

    return res;
}

void VStr_clear( VecStr *strv )
{
    VStr_resize( strv, 0 );
}

int VStr_find( const VecStr *strv, const char *str )
{
    int i;
    for ( i=0; (i<VStr_size(strv)); i++ )
        if ( strcmp( str, VStr_get(strv,i) ) == 0 )
            return i;

    return INT_MAX;
}

int VStr_findSubStr( const VecStr *strv, const char *str )
{
    int i;
    for ( i=0; (i<VStr_size(strv)); i++ )
        if ( strstr( VStr_get(strv,i), str ) == 0 )
            return i;

    return INT_MAX;
}

void VStr_writeTo( VecStr *strv, const char *fileName )
{
    FILE *f = fopen( fileName, "w" );
    int i;
    for ( i=0; (i<VStr_size(strv)); i++ )
        fprintf( f, "%s\n", VStr_get( strv,i) );
    fclose( f );
}

void VStr_readFrom( VecStr *strv, const char *fileName, const char ignoreEmptyLines )
{
    char line[LINE_SIZE], *s;

    FILE *f = fopen( fileName, "r" );
    if (!f) {
        fprintf( stderr, "could not open %s.\n", fileName );
        // abort();
        exit( EXIT_FAILURE );
    }

    while ( (s=fgets(line,FILE_NAME_SIZE,f)) ) {
        strRemoveSpsEol( s );
        if ( ignoreEmptyLines && (!strlen(s)) )
            continue;

        VStr_pushBack( strv, s );
    }

    fclose(f);
}

int VStr_size( const VecStr *strv )
{
    return strv->size;
}

void VStr_pushBack( VecStr *strv, const char *str )
{
#ifdef DEBUG
    VStr_check_valid_strv( strv );
    assert( str );
#endif
    if ( strv->size+1 > strv->capacity )
        VStr_increase_capacity_to( strv, strv->capacity*2 );

    strncpy( strv->sv[strv->size], str, strv->strSize );
    strv->size++;
}


const char *VStr_get( const VecStr *strv, int pos )
{
#ifdef DEBUG
    VStr_check_valid_strv( strv );
    assert( pos >= 0 );
    assert( pos < strv->size );
#endif

    return strv->sv[pos];
}

void VStr_resize( VecStr *strv, int newSize )
{
#ifdef DEBUG
    VStr_check_valid_strv( strv );
    assert( newSize>=0 );
#endif
    if ( newSize > strv->capacity )
        VStr_increase_capacity_to( strv, newSize );
    strv->size = newSize;
}

#ifdef DEBUG
void VStr_check_valid_strv( const VecStr *strv )
{
    assert( strv );
    assert( strv->s );
    assert( strv->sv );
    assert( strv->sv[0] == strv->s );
    assert( strv->capacity > 0 );
    assert( strv->size >= 0 );
    assert( strv->strSize > 0 );
    assert( strv->capacity >= strv->size );
}
#endif

void VStr_increase_capacity_to( VecStr *strv, int newCapacity )
{
#ifdef DEBUG
    VStr_check_valid_strv( strv );
    assert( newCapacity > strv->capacity );
#endif

    char **sv = xmalloc( sizeof(char*)*newCapacity );
    char *s = xcalloc( newCapacity*strv->strSize, sizeof(char) );
    sv[0] = s;
    int i;
    for ( i=1; (i<newCapacity); i++ )
        sv[i] = sv[i-1] + strv->strSize;

    memcpy( s, strv->s, sizeof(char)*strv->size*strv->strSize );

    free( strv->sv );
    free( strv->s );
    strv->sv = sv;
    strv->s = s;

    strv->capacity = newCapacity;
}

void VStr_set( VecStr *strv, int pos, const char *str )
{
#ifdef DEBUG
    VStr_check_valid_strv( strv );
    assert( pos>=0 );
    assert( pos<strv->size );
    assert( str );
#endif
    strncpy( strv->sv[pos], str, strv->strSize );
}

char **VStr_ptr( VecStr *strv )
{
#ifdef DEBUG
    VStr_check_valid_strv( strv );
#endif
    return(strv->sv);
}

void VStr_free( VecStr **_strv )
{
    VecStr *strv = *_strv;

#ifdef DEBUG
    VStr_check_valid_strv( strv );
#endif

    free( strv->s);
    free( strv->sv);
    free( strv );

    *_strv = NULL;

}

