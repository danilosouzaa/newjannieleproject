/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef LONG_COMPL_PATH_HPP
#define LONG_COMPL_PATH_HPP

#include "instance.h"

typedef struct _LongestComplPath LongestComplPath;

LongestComplPath *LongCP_create( const Instance *inst );

void LongCP_solve( LongestComplPath *lcp );

const int *LongCP_get( const LongestComplPath *lcp );

void LongCP_free( LongestComplPath **_lcp );

#endif
