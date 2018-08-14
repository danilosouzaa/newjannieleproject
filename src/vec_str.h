/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef VEC_STR_H
#define VEC_STR_H

typedef struct _VecStr VecStr;

/* creates a vector of strings, each line has strSize of space */
VecStr *VStr_create( int strSize );

/* increases/decreases number of strings */
void VStr_resize( VecStr *strv, int newSize );

/* adds another string at the end of the vector */
void VStr_pushBack( VecStr *strv, const char *str );

/* gets the string at position pos  */
const char *VStr_get( const VecStr *strv, int pos );

/* sets the contents of position pos */
void VStr_set( VecStr *strv, int pos, const char *str );

/* returns a pointer to the string matrix */
char **VStr_ptr( VecStr *strv );

/* returns the current size of the string vector */
int VStr_size( const VecStr *strv );

/* clears the string vector */
void VStr_clear( VecStr *strv );

/* finds for a substring in the strings vector */
int VStr_find( const VecStr *strv, const char *str );
int VStr_findSubStr( const VecStr *strv, const char *str );
void VStr_readFrom( VecStr *strv, const char *fileName, const char ignoreEmptyLines );
void VStr_writeTo( VecStr *strv, const char *fileName );
void VStr_free( VecStr **_strv );

#endif

