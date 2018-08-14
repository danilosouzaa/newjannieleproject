/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef MS_SOLVER_MIP_H
#define MS_SOLVER_MIP_H

#include "instance.h"
#include "solution.h"
#include "mode_set.h"
#include "lp.h"

typedef struct _MSM_Solver MSM_Solver;

/* created mode min solver for jobs in firstJob...lastJob */
MSM_Solver *MSM_create( const Instance *inst, int firstJob, int lastJob, double timeLeft );

/* solves selecting an initial set of modes */
int MSM_solve( MSM_Solver *msm, double timeLeft);

int MSM_getSIdx( const MSM_Solver *msm, int j, int m);
LinearProgram *MSM_lp(const MSM_Solver *msm);
/* "controled" change in modes for diversification
 * residency[j][m] must inform how many times job j
 * was in mode m, returns true if a feasible
 * and sufficiently different solution was found */
char MSM_changeModes( MSM_Solver *msm, const ModeSet *current, int minChanges,  int maxChanges, const int **residency );
/* solution for solve or changeModes */
const ModeSet *MSM_modes( const MSM_Solver *msm );

void MSM_roundingHeuristc(MSM_Solver *msm, Solution *sol);

void MSM_roundVar(MSM_Solver *msm, Solution *sol);

/* frees memory */
void MSM_free( MSM_Solver **_msm );

#endif /* !MS_SOLVER_MIP_H */

