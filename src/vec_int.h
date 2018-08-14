/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef VEC_INT_HEADER
#define VEC_INT_HEADER

/* opaque data structure */
typedef struct _VecInt VecInt;

/* creates an empty vector */
VecInt *VInt_create();

/* creates an empty vector with an initial capacity */
VecInt *VInt_createIniCap( const int iniCap );

/* creates a vector and fills with some initial value */
VecInt *VInt_createFill( const int size, const int value );

/* inserts a new element at the end of the vector */
void VInt_pushBack( VecInt *v, int value );

/* makes this vector empty */
void VInt_clear( VecInt *v );

/* returns the number of elements in this vector */
int VInt_size( const VecInt *v );

/* returns the current capacity of this vector (i.e. including pre-allocated space) */
int VInt_capacity( const VecInt *v );

/* gets element at position pos */
int VInt_get( const VecInt *v, const int pos );

/* gets a pointer to the start of the vector */
int *VInt_getPtr( VecInt *v );

/* sets element at position pos */
void VInt_set( VecInt *v, const int pos, const int value );

/* searches for "value" in this vector, returning its position,  if not found returns INT_MAX. */
int VInt_find( const VecInt *v, const int value );

/* frees all memory for memory */
void VInt_free( VecInt **_v );



#endif /* VEC_INT_HEADER */

