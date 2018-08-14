/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef VNS_H
#define VNS_H

#include "instance.h"
#include "neighborhood.h"
#include "solution.h"
#include "test.h"

typedef struct _VNS VNS;

/*creates a solver to make the allocation of jobs*/
VNS *VNS_create( const Instance *inst, Solution *sol,  Neighborhood *neighborhood,  char **argv, int argc );

/*run the solver*/
void VNS_run_general_vnd(VNS *vns, Neighborhood* neighborhood, double timeRem, int nNeighbor, int firstImprovement, Test *test, char **argv, int argc );
void VNS_rna(VNS *vns,  Solution *current, Neighborhood* neighborhood, double timeRem, Test *test);
void VNS_run_RNA(VNS *vns, Neighborhood* neighborhood, double timeRem, int nNeighbor, Test *test, char **argv, int argc );
void VNS_run_LAHC_shake2(VNS *vns, Neighborhood* neighborhood, double timeRem, int nNeighbor,  Test *test, char **argv, int argc );
void VNS_run_LAHC_smartshake(VNS *vns, Neighborhood* neighborhood, double timeRem, int nNeighbor,  char **argv, int argc );
void VNS_run_reduced(VNS *vns, Neighborhood* neighborhood, double timeRem, int nNeighbor, int firstImprovement, Test *test, char **argv, int argc );
void VNS_run_RNA_shake2(VNS *vns, Neighborhood* neighborhood, double timeRem, int nNeighbor,  Test *test, char **argv, int argc );
void VNS_run_RNA_shake(VNS *vns, Neighborhood* neighborhood, double timeRem, int nNeighbor,  Test *test, char **argv, int argc );

int VNS_getItRNA(VNS *vns);
void VNS_checkArgs(VNS *vns, char **argv, int argc);
int VNS_getNMoves(VNS *vns);
int VNS_getAllNeighbor(VNS *vns);
void VNS_setNChanges(VNS *vns, int j);
int VNS_getNChanges(VNS *vns, int j);
double VNS_getPercRS(VNS *vns);
double VNS_getPerc(VNS *vns);
int VNS_getLfa(VNS *vns);
int VNS_getDivTM(VNS *vns);
int VNS_getDivRJ(VNS *vns);
int VNS_getDivRM(VNS *vns);
int VNS_getDivTJ(VNS *vns);
int VNS_getItLAHC(VNS *vns);
int VNS_getNSizeSamplingShake(VNS *vns);

void VNS_increasingResidencyJobInMode(VNS* vns, Solution* current);
void VNS_increasingResidencyJobInSequence(VNS* vns, Solution* current);
Cost VNS_penaltyModes(VNS* vns, Neighborhood* neighborhood, Cost foCurrent);
Cost VNS_penaltyJobs(VNS* vns, Neighborhood* neighborhood, Cost foCurrent);
Cost VNS_penaltyTransModes(VNS* vns, Neighborhood* neighborhood, Cost foCurrent);
Cost VNS_penaltyTransJobs(VNS* vns, Neighborhood* neighborhood, Cost foCurrent);
void VNS_increasingTransitivityOfModes(VNS* vns, Neighborhood* neighborhood);
void VNS_increasingTransitivityOfSequence(VNS* vns, Neighborhood* neighborhood);
/* frees memory used by VNS */
void VNS_free( VNS **_vns );

#endif

