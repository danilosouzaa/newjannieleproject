/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef MS_SOLVER_LAHC_H
#define MS_SOLVER_LAHC_H

#include "instance.h"
#include "neighborhood.h"
#include "solution.h"

typedef struct _MSSolverLAHC MSSolverLAHC;

/*creates a solver to solve the modes for each project*/
MSSolverLAHC *MS_createByProj( const Instance *inst, int sizeList, int firstJob, int lastJob );
MSSolverLAHC *MS_create( const Instance *inst, int sizeList);

/*run the solver to instance by project*/
void MS_runByProj(const Instance *inst, Solution *sol, int LFA, int iterat);
/*run the solver to instance*/
void MS_run(MSSolverLAHC *msLahc, int iterat);

/*returns the best ModeSet*/
ModeSet *MS_bestModes( const MSSolverLAHC *solver );

/*returns the neighborhood*/
Neighborhood *MS_getNeighborhood( const MSSolverLAHC *solver );

/* frees memory used by MSSolverLAHC */
void MS_free( MSSolverLAHC **_msLahc );

#endif

