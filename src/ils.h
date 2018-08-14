/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef ILS_H
#define ILS_H

#include "instance.h"
#include "neighborhood.h"
#include "solution.h"

typedef struct _ILS ILS;

ILS *ILS_create( const Instance *inst, Solution* sol);

void ILS_run(ILS *ils);

void ILS_free( ILS **_ils );

#endif
