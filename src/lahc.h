/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef LAHC_H
#define LAHC_H

#include "instance.h"
#include "neighborhood.h"
#include "solution.h"
#include "vns.h"

typedef struct _LAHC LAHC;

/*creates a solver to make the allocation of jobs*/
LAHC *LAHC_create( const Instance *inst, Solution* sol, Neighborhood* neighborhood, char **argv, int argc );

/* -lfa*/
void LAHC_checkArgs(LAHC *lahc, char **argv,int argc);

void LAHC_setF(LAHC *lahc, int idxF, Cost value);
int LAHC_getF(LAHC *lahc, int idxF);
int LAHC_getLfa(LAHC *lahc);
int LAHC_getNChangesModes(LAHC *lahc, int idx);

void Test_writeResultLAHC(Test *test, char **argv, Cost initFO, Neighborhood *neigh, Solution* sol, int nSol, int lfa, int seed, int nCostList, int nDiversification, int nStayDiversification, int nWOImprove, double perc, int sw, double lr, int it, float psw, int tm, int tj, int rm, int rj);
int LAHC_getNDiversification(LAHC *lahc);
int LAHC_getNStayDiversification(LAHC *lahc);
int LAHC_getSW(LAHC *lahc);
int LAHC_getNWOImprove(LAHC *lahc);
int LAHC_getNCostList(LAHC *lahc);
int LAHC_getNThread(LAHC *lahc);

float LAHC_getPSW(LAHC *lahc);
int LAHC_getRJ(LAHC *lahc);
int LAHC_getRM(LAHC *lahc);
int LAHC_getTJ(LAHC *lahc);
int LAHC_getTM(LAHC *lahc);

double LAHC_getPerc(LAHC *lahc);

LearningAutomata *LAHC_getLA(LAHC *lahc);
int LAHC_getItUp(LAHC *lahc);
/*run the solver*/
void LAHC_run(LAHC *lahc, Neighborhood* neighborhood, double timeRem, Test *test, char * nameInst);
void LAHC_run_parallel(LAHC *lahc, Neighborhood* neighborhood, double timeRem, Solution *globalBest, Test *test);
void LAHC_run_nIt(LAHC *lahc, Neighborhood* neighborhood, double timeRem, int nIterations, Test *test, char* nameInst);
void LAHC_run_vns(LAHC *lahc, VNS* vns, Neighborhood* neighborhood, double timeRem,  int nIterations,  char* nameInst);
void LAHC_increasingResidency(LAHC* lahc, Solution* current);
void LAHC_increasingTransitivity(LAHC* lahc, Neighborhood* neighborhood);
void LAHC_updateF(LAHC *lahc);
Cost LAHC_penalty(LAHC* lahc, Neighborhood* neighborhood, Cost foCurrent);

/* frees memory used by SA */
void LAHC_free( LAHC **_lahc );

#endif

