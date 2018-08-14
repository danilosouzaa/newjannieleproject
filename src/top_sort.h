/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef TOP_SORT_H_INCLUDED
#define TOP_SORT_H_INCLUDED

#include "instance.h"

typedef struct _TopSort TopSort;

TopSort *TopSort_create( const Instance *inst );

void TopSort_computeProjectTS( TopSort *ts, int idxProject );

int TopSort_nTS( const TopSort *ts );

int *TopSort_ts( const TopSort *ts );

void TopSort_free( TopSort **_topSort );

#endif
