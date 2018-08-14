/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef VEC_DBL_HEADER
#define VEC_DBL_HEADER

/* opaque data structure */
typedef struct _VecDbl VecDbl;

/* creates an empty vector */
VecDbl *VDbl_create();

/* creates an empty vector with an initial capacity */
VecDbl *VDbl_createIniCap( const int iniCap );

/* creates a vector and fills with some initial value */
VecDbl *VDbl_createFill( const int size, const double value );

/* inserts a new element at the end of the vector */
void VDbl_pushBack( VecDbl *v, double value );

/* makes this vector empty */
void VDbl_clear( VecDbl *v );

/* returns the number of elements in this vector */
int VDbl_size( const VecDbl *v );

/* returns the current capacity of this vector (i.e. including pre-allocated space) */
int VDbl_capacity( const VecDbl *v );

/* gets element at position pos */
double VDbl_get( const VecDbl *v, const int pos );

/* gets a pointer to the start of the vector */
double *VDbl_getPtr( VecDbl *v );

/* sets element at position pos */
void VDbl_set( VecDbl *v, const int pos, const double value );

/* searches for "value" in this vector, returning its position,  if not found returns INT_MAX. */
int VDbl_find( const VecDbl *v, const double value );

/* frees all memory for memory */
void VDbl_free( VecDbl **_v );

#endif /* VEC_INT_HEADER */

