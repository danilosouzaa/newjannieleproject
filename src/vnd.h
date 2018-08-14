/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef VND_H
#define VND_H

#include "instance.h"
#include "neighborhood.h"
#include "solution.h"
#include "test.h"

typedef struct _VND VND;

/*creates a solver to make the allocation of jobs*/
VND *VND_create( const Instance *inst, Solution *sol, char **argv, int argc);

/*run the solver*/
void VND_runStoc(VND *down, double timeRem, int firstImprovement);
//void VND_runDet(VND *vnd, Neighborhood *neighborhood, double timeRem, int nNeighbor, int firstImprovement, int *nChangesModes, int **nTimesJobOnModes, Test *test);
void VND_runDet(VND *vnd, Neighborhood *neighborhood, double timeRem, int nNeighbor, int firstImprovement, Test *test);

void VND_checkArgs(VND *vnd,  char **argv, int argc);

Neighborhood *VND_getNeighborhood(VND *vnd);

Solution *VND_getSolThread(VND *vnd, int idx);

int *VND_getNChangesModes(VND *vnd);
int **VND_getNTimesJobOnModes(VND *vnd);
/* frees memory used by VND */
void VND_free( VND **_vnd );

#endif

