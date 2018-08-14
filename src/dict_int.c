#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "dict_int.h"
#include "macros.h"
#include "vec_str.h"

static const unsigned int hashval[] = { 11, 269, 3, 7, 31, 37, 131, 13, 17, 647, 653, 89, 97, 101, 39, 149, 151, 157, 821, 257, 263, 389, 397, \
                                        457, 461, 463, 331, 337, 347, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 9, 53, 59
                                      };

static const unsigned int nHashvalues = sizeof(hashval)/sizeof(int);

int cmp_string( const void *p1, const void *p2 );

static unsigned int str_hash( const char *str, const unsigned int hashSize )
{
    const unsigned int len = (unsigned int) MIN( nHashvalues, strlen(str) );

    unsigned int sum = 0, i;

    for ( i=0; (i<len) ; i++ )
        sum += (unsigned int) str[i]*hashval[i];

    return sum%hashSize;
}


typedef struct {
    char str[256];
    int value;
    int keyPos;
} Dict_Bucket_int;

struct _Dict_int {
    Dict_Bucket_int **cont;
    int *rowSize;
    int *rowCap;
    int defValue;
    unsigned int hashSize;
    VecStr *keys;
};

DictInt *DInt_create( unsigned int hashSize, int defaultValue )
{
    DictInt *res;
    {
        res = malloc( sizeof(DictInt) );
        if (!res) {
            fprintf( stderr, "ERROR: no more memory available. at: %s:%d\n", "containers.c", 48 );
            abort();
        };
    };
    {
        res->cont = calloc( hashSize, sizeof(Dict_Bucket_int *) );
        if (!res->cont) {
            fprintf( stderr, "ERROR: no more memory available. at: %s:%d\n", "containers.c", 48 );
            abort();
        };
    };
    {
        res->rowSize = calloc( hashSize, sizeof(int) );
        if (!res->rowSize) {
            fprintf( stderr, "ERROR: no more memory available. at: %s:%d\n", "containers.c", 48 );
            abort();
        };
    };
    {
        res->rowCap = calloc( hashSize, sizeof(int) );
        if (!res->rowCap) {
            fprintf( stderr, "ERROR: no more memory available. at: %s:%d\n", "containers.c", 48 );
            abort();
        };
    };
    res->defValue = defaultValue;
    res->keys = VStr_create( 256 );
    res->hashSize = hashSize;
    return res;
}

void DInt_set( DictInt *dict, const char *key, const int value )
{
    unsigned int hashPos = str_hash( key, dict->hashSize );
    int i;
    char found = 0;
    for ( i=0 ; (i<dict->rowSize[hashPos]) ; i++ ) {
        if ( strcmp( key, dict->cont[hashPos][i].str) == 0 ) {
            found = 1;
            break;
        }
    }
    if (found) {
        dict->cont[hashPos][i].value = value;
        return;
    }
    if ( dict->rowSize[hashPos]+1 > dict->rowCap[hashPos] ) {
        if (dict->rowSize[hashPos]==0) {
            {
                dict->cont[hashPos] = malloc( sizeof(Dict_Bucket_int)*(64) );
                if (!dict->cont[hashPos]) {
                    fprintf( stderr, "ERROR: no more memory available. at: %s:%d\n", "containers.c", 48 );
                    abort();
                };
            };
            dict->rowCap[hashPos] = 64;
        } else {
            dict->rowCap[hashPos]*=2;
            Dict_Bucket_int *bigger = realloc( dict->cont[hashPos], sizeof(Dict_Bucket_int)*dict->rowCap[hashPos] );
            if (!bigger) {
                fprintf( stderr, "ERROR: no more memory available" );
                abort();
                //exit(1);
            }
            dict->cont[hashPos] = bigger;
        }
    }
    strncpy( dict->cont[hashPos][dict->rowSize[hashPos]].str, key, 256 );
    dict->cont[hashPos][dict->rowSize[hashPos]].value = value;
    dict->cont[hashPos][dict->rowSize[hashPos]].keyPos = VStr_size( dict->keys );
    dict->rowSize[hashPos]++;
    VStr_pushBack( dict->keys, key );
}

void DInt_remove( DictInt *dict, const char *key )
{
    unsigned int hashPos = str_hash( key, dict->hashSize );
    int i, pBucket = -1;
    Dict_Bucket_int *bucket = ((void*)0);
    for ( i=0 ; (i<dict->rowSize[hashPos]) ; i++ ) {
        if ( strcmp( key, dict->cont[hashPos][i].str) == 0 ) {
            pBucket = i;
            bucket = &dict->cont[hashPos][i];
            break;
        }
    }
    if ( pBucket == -1 ) return;
    if ( bucket->keyPos != (VStr_size(dict->keys)-1) ) {
        const char *keyToUpdate = VStr_get(dict->keys, VStr_size(dict->keys)-1);
        VStr_set( dict->keys, bucket->keyPos, keyToUpdate );
        {
            unsigned int hashPosK = str_hash( keyToUpdate, dict->hashSize );
            Dict_Bucket_int *bucketK = ((void*)0);
            int ik;
            for ( ik=0 ; (ik<dict->rowSize[hashPosK]) ; ik++ ) {
                if ( strcmp( keyToUpdate, dict->cont[hashPosK][ik].str) == 0 ) {
                    bucketK = dict->cont[hashPosK] + ik;
                    break;
                }
            }
            bucketK->keyPos = bucket->keyPos;
        }
    }
    VStr_resize( dict->keys, VStr_size(dict->keys)-1 );
    if ( pBucket != dict->rowSize[hashPos]-1 )
        dict->cont[hashPos][pBucket] = dict->cont[hashPos][dict->rowSize[hashPos]-1];
    dict->rowSize[hashPos]--;
}

void DInt_iset( DictInt *dict, const int key, const int value )
{
    char str[32];
    sprintf( str, "%d", key );
    DInt_set( dict, str, value );
}

int DInt_iget( const DictInt *dict, const int key )
{
    char str[32];
    sprintf( str, "%d", key );
    return DInt_get( dict, str);
}

int DInt_get( const DictInt *dict, const char *str )
{
    unsigned int hashPos = str_hash( str, dict->hashSize );
    int i;
    for ( i=0 ; (i<dict->rowSize[hashPos]) ; i++ ) {
        if ( strcmp( str, dict->cont[hashPos][i].str) == 0 )
            return dict->cont[hashPos][i].value;
    }
    return dict->defValue;
}

void DInt_iremove( DictInt *dict, int key )
{
    char strKey[256];
    sprintf( strKey, "%d", key );
    DInt_remove(dict, strKey);
}

int DInt_size( const DictInt *dict )
{
    return VStr_size( dict->keys );
}

const char *DInt_key( const DictInt *dict, int pos )
{
    return VStr_get( dict->keys, pos );
}

int DInt_ikey( const DictInt *dict, int pos )
{
    return atoi(VStr_get(dict->keys, pos));
}

void DInt_clear( DictInt *dict )
{
    memset( dict->rowSize, 0, dict->hashSize*sizeof(int) );
    VStr_resize( dict->keys, 0 );
}

void DInt_cpy( DictInt *target, const DictInt *source )
{
    DInt_clear( target );
    int i;
    for ( i=0 ; (i<DInt_size(source)) ; ++i ) DInt_set( target, DInt_key(source,i), DInt_get(source, DInt_key(source,i)) );
}

void DInt_free( DictInt **dict )
{
    DictInt *d = *dict;
    unsigned int i;
    for (i=0; (i<d->hashSize); ++i) if (d->rowCap[i]) free( d->cont[i] );
    free( (*dict)->cont );
    free( (*dict)->rowSize );
    free( (*dict)->rowCap );
    VStr_free( &(*dict)->keys );
    free( *dict );
}

