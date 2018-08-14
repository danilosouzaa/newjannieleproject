/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef VEC_CHAR_HEADER_DEFINED
#define VEC_CHAR_HEADER_DEFINED

/* opaque data structure */
typedef struct _VecChar VecChar;

/* creates an empty vector */
VecChar *VChar_create();

/* creates an empty vector with an initial capacity */
VecChar *VChar_createIniCap( const int iniCap );

/* creates a vector and fills with some initial value */
VecChar *VChar_createFill( const int size, const char value );

/* inserts a new element at the end of the vector */
void VChar_pushBack( VecChar *vec, const char value );

/* makes this vector empty */
void VChar_clear( VecChar *v );

/* returns the number of elements in this vector */
int VChar_size( const VecChar *v );

/* returns the current capacity of this vector (i.e. including pre-allocated space) */
int VChar_capacity( const VecChar *v );

/* gets element at position pos */
char VChar_get( const VecChar *v, const int pos );

/* gets a pointer to the start of the vector */
char *VChar_getPtr( VecChar *v );

/* sets element at position pos */
void VChar_set( VecChar *v, const int pos, const char value );

/* searches for "value" in this vector, returning its position,  if not found returns INT_MAX. */
int VChar_find( const VecChar *v, const int value );

/* frees all memory for memory */
void VChar_free( VecChar **v );

#endif /* VEC_CHAR_HEADER_DEFINED */

