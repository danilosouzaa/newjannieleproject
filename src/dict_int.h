/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef DICT_INT_H
#define DICT_INT_H

typedef struct _Dict_int DictInt;

DictInt *DInt_create( unsigned int hashSize, int defaultValue );
void DInt_set( DictInt *dict, const char *key, const int value );
void DInt_iset( DictInt *dict, const int key, const int value );
void DInt_remove( DictInt *dict, const char *key );
void DInt_iremove( DictInt *dict, int key );
int DInt_get( const DictInt *dict, const char *str );
int DInt_size( const DictInt *dict );
const char *DInt_key( const DictInt *dict, int pos );
int DInt_ikey( const DictInt *dict, int pos );
int DInt_iget( const DictInt *dict, const int key );
void DInt_clear( DictInt *dict );
void DInt_cpy( DictInt *target, const DictInt *source );
void DInt_free( DictInt **dict );

#endif /* !DInt_H */

