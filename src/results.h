/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * mip_compact: routines for the creation of compact MIP formulations
 *              for the whole problem or for just one project
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef RESULTS
#define  RESULTS
#include "macros.h"
#include "vec_str.h"
#include "vec_int.h"
#include "vec_double.h"
#include "dict_int.h"
#include "tokenizer.h"
#include "vec_str.h"
#include "str_utils.h"
#include "parameters.h"
#include "solution.h"
#include "lp.h"

typedef struct _Results Results;

Results *Res_create();
void Res_free( Results **res );
void Res_writeResults(Parameters *par, LinearProgram *mip, const Solution* sol, Results *res, double bo, FILE *fp, char **argv, int argc, double startT, int sumtpd, int  sumtms);
double Res_getMaxViol(Results *res, int nround, int cut);
double Res_getMinViol(Results *res, int nround, int cut);
int Res_getNMaxElementsConf(Results *res, int cut);
int Res_getMaxViolCutA(Results *res, int nround, int cut, int i);
double Res_getMaxViolCutB(Results *res, int nround, int cut, int i);
int Res_getMinViolCutA(Results *res, int nround, int cut, int i);
double Res_getMinViolCutB(Results *res, int nround, int cut, int i);
void Res_setMaxViolCut(Results *res, int nround, int cut, int i, int a, double b);
void Res_setMinViolCut(Results *res, int nround, int cut, int i, int a, double b);
int Res_getNElementsCuts(Results *res, int nround, int cut);
int Res_getNElementsMinViol(Results *res, int nround, int cut);
int Res_getNElementsMaxViol(Results *res, int nround, int cut);
int Res_getNMaxElementsCuts(Results *res, int nround, int cut);
int Res_getNCutsTotal(Results *res, int nround, int cut);
double Res_getTCutsTotal(Results *res, int nround, int cut);
int Res_getNMinElementsCuts(Results *res, int nround, int cut);
void Res_setTCutsTotal(Results *res, int nround,int cut, double value);
void Res_setNElementsMaxViol(Results *res, int nround, int cut, int value);
void Res_setNElementsMinViol(Results *res, int nround, int cut, int value);
void Res_setNCutsTotal(Results *res, int nround, int cut, int value);
void Res_setNElementsCuts(Results *res, int nround, int cut, int value);
void Res_setNMaxElementsCuts(Results *res, int nround, int cut, int value);
void Res_setNMinElementsCuts(Results *res, int nround, int cut, int value);
void Res_setMaxViol(Results *res, int nround, int cut, double value);
void Res_setMinViol(Results *res, int nround, int cut, double value);
void Res_setSumViol(Results *res, int nround, int cut, double value);
void Res_setNSumAllVarWithConf(Results *res, int cut, int value);
void Res_setNSumAllElementsConf(Results *res, int cut, int value);
void Res_setNMaxElementsConf(Results *res, int cut, int value);
void Res_setIdxStartCutsMaxViol(Results *res, int cut);
void Res_setIdxStartCutsMinViol(Results *res, int cut);
void Res_initRound(Results *res);
int Res_getRound(Results *res);
void Res_setRound(Results *res, int value);
double Res_getSumViol(Results *res, int nround, int cut);

void Res_setUsoUR(Results *res, int nround, int value);
void Res_setUsoURN(Results *res, int nround, int value);
void Res_setUsoUM(Results *res, int nround, int value);
void Res_setUsoUCI(Results *res, int nround, int value);
void Res_setUsoUCP(Results *res, int nround, int value);

void Res_setTotalUsoUR(Results *res, int value);
void Res_setTotalUsoURN(Results *res, int value);
void Res_setTotalUsoUM(Results *res, int value);
void Res_setTotalUsoUCI(Results *res, int value);
void Res_setTotalUsoUCP(Results *res, int value);
void Res_setTotalTimeSep(Results *res, double value);
double Res_getTotalTimeSep(Results *res);
void Res_setTotalTimeRemAdd(Results *res, double value);
double Res_getTotalTimeRemAdd(Results *res);
void Res_setTotalTimeOpt(Results *res, double value);
void Res_setTotalTimeCGraph(Results *res, double value);
int Res_getUsoUR(Results *res, int nround);
int Res_getUsoURN(Results *res, int nround);
int Res_getUsoUM(Results *res, int nround);
int Res_getUsoUCI(Results *res, int nround);
int Res_getUsoUCP(Results *res, int nround);
int Res_getTotalUsoUR(Results *res);
int Res_getTotalUsoURN(Results *res);
int Res_getTotalUsoUM(Results *res);
int Res_getTotalUsoUCI(Results *res);
int Res_getTotalUsoUCP(Results *res);
void Res_setNJumps(Results *res, int value);
int Res_getNJumps(Results *res);



int Res_getTW(Results *res);
double Res_getTWP(Results *res);
int Res_getS(Results *res);
int Res_getE(Results *res);


int Res_getMinTW(Results *res);
double Res_getMinTWP(Results *res);
int Res_getMinS(Results *res);
int Res_getMinE(Results *res);

void Res_setTW(Results *res, int nMaxTW, double nMaxTWP, int nS, int nE);
void Res_setMinTW(Results *res, int nMinTW, double nMinTWP, int nmS, int nmE);
#endif // RESULTS


