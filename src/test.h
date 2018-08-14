#ifndef TEST_H_INCLUDED
#define TEST_H_INCLUDED

#include <time.h>
#include "solution.h"
typedef struct _Test Test;

Test *Test_create(int nNeigh, const Instance* inst);
void Test_readAnalisysNeigh_offline(Test *test,  char **argv);
void Test_writeAnalisysNeigh(Test *test,  char **argv);

//void Test_writeResultVND(Test *test, char **argv, Cost initFO, int *assortment, Solution* sol);
void Test_writeResultVND(Test *test, char **argv, Cost initFO, int *assortment, Solution* sol, int fi);
//void Test_writeResultVNS(Test *test, char **argv, Cost initFO, int *assortment, int first, int type, Solution* sol, int itRNA, int nMoves, int tm, int tj, int rm, int rj, int nSol);
void Test_writeResultVNS(Test *test, char **argv, Cost initFO, int *assortment, int first, int type, Solution* sol,  int lfa, int itRNA, int itLAHC, int nMoves, int nSizeSamplingShake, double perc, double percRS, int tm, int tj, int rm, int rj, int nSol);
//void Test_writeResultLAHC(Test *test, char **argv, Cost initFO, Neighborhood *neigh, Solution* sol, int lfa);

void Test_writeResultSA(Test *test, char **argv, Cost initFO, int *assortment, Solution* sol);

void Test_setImproveNeigh(Test *test, int idxVet, int value);

void Test_incrementImproveNeigh(Test *test, int idxVet, int value);

void Test_setVisitNeigh(Test *test, int idxVet, int value);

void Test_incrementVisitNeigh(Test *test, int idxVet, int value);

void Test_setTimeNeigh(Test *test, int idxVet, double value);

void Test_incrementTimeNeigh(Test *test, int idxVet, double value);

void Test_setTotalTime(Test *test, double value);

void Test_setImproveFO(Test *test, int idxVet, long int value);

void Test_setCurrentFO(Test *test, Cost value);

void Test_setCurrentTime(Test *test, double _time);

void Test_setT(Test *test, int t);

void Test_setSAmax(Test *test, int samax);

void Test_setAlpha(Test *test, double alpha);

void Test_incrementImproveFO(Test *test, int idxVet, Cost value);

double Test_getAlpha(Test *test);

int Test_getT(Test *test);

int Test_getSAmax(Test *test);

int Test_getImproveNeigh(Test *test, int idxVet);

int Test_getVisitNeigh(Test *test, int idxVet);

double Test_getTimeNeigh(Test *test, int idxVet);

double Test_getTotalTime(Test *test);

Cost Test_getImproveFO(Test *test, int idxVet);

Cost Test_getCurrentFO(Test * test);

double Test_getCurrentTime(Test *test);

int Test_getCurrentNeigh(Test *test);

void Test_setCurrentNeigh(Test *test, int nNeigh);

void Test_callTest(Test *test, int idNeigh, double _timeNeigh, Cost bestFO);

void Test_free( Test **_test );

int Test_fileExists(char fileName[]);

//void Test_improvNeigh(Instance *inst, Neighborhood *neigh, double timeRem, int nNeighbor, int firstImprovement, Test *test, char **argv);

#endif // TEST_H_INCLUDED
