/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef PROJ_DECOMP_H
#define PROJ_DECOMP_H

#include "instance.h"

typedef struct _ProjDecomp ProjDecomp;

ProjDecomp *prjdc_create( const Instance *instance, int tpdSum );

void prjdc_free( ProjDecomp **_prjdc );

#endif /* !PROJ_DECOMP_H */

