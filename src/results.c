
/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <omp.h>
#include <time.h>
#include "results.h"
#include "cut_pool.h"

#define VERBOSE  2

struct _Results {
    int nRounds; //número de rounds para encontrar cortes
    VecDbl ** tCutsTotal;
    //double **tCutsTotal; //total de tempo para gerar os cortes cortes em cada round para cada tipo de corte
    VecInt ** nCutsTotal;
    //int **nCutsTotal; //para cada round e para cada tipo o numero total de cortes
    VecInt ** nElementsCuts;
    //int **nElementsCuts;     // para cada round para cada tipo o número total de elementos (soma dos elementos)
    VecInt ** nMaxElementsCuts;
    VecInt ** nMinElementsCuts;
    //int **nMaxElementsCuts; //para cada round e para cada tipo o maior corte (número de elementos)
    //int **nMinElementsCuts; //para cada round e para cada tipo o menor corte (número de elementos)
    VecInt ** idxStartCutsMaxViol;
    VecInt ** maxViolCutIdx;
    VecDbl ** maxViolCutCoef;
    //IntDblPair ***maxViolCut; //para cada round e para cada tipo o corte com maior violação (elementos)
    VecInt ** nElementsMaxViol;
    //int **nElementsMaxViol; // para cada round e para cada tipo o número de elementos da maior violação.
    VecInt ** idxStartCutsMinViol;
    VecInt ** minViolCutIdx;
    VecDbl ** minViolCutCoef;
    //IntDblPair ***minViolCut; //para cada round e para cada tipo o corte com maior violação (elementos)
    VecInt ** nElementsMinViol;
    //int **nElementsMinViol; // para cada round e para cada tipo o número de elementos da maior violação.

    VecDbl ** maxViol;
    VecDbl ** minViol;
    VecDbl ** sumViol;
    //double **maxViol; //para cada round e para cada tipo de corte a maior violação;
    //double **minViol; //para cada round e para cada tipo de corte a menor violação;
    //double **sumViol; //para cada round e para cada tipo de corte a soma das violações;
    VecInt *nSumAllVarWithConf;
    VecInt *nSumAllElementsConf;
    VecInt *nMaxElementsConf;
    //int *nSumAllVarWithConf; // número total de variáveis com conflitos em todas as rodadas
    //int *nSumAllElementsConf; //número total de elementos dos conflitos em todas as rodadas
    //int *nMaxElementsConf; //tamanho do maior conflito encontrado, ou seja, maior número de elementos.

    VecInt *usoUR;
    VecInt *usoURN;
    VecInt *usoUM;
    VecInt *usoUCI;
    VecInt *usoUCP;

    int totalUsoUR;
    int totalUsoURN;
    int totalUsoUM;
    int totalUsoUCI;
    int totalUsoUCP;
    double totalTimeSep;
    double totalTimeOpt;
    double totalTimeRemAdd;
    double totalTimeCGraph;
    int nMaxWindowWithCutEnergeticViolated;
    double nMaxWindowWithCutEnergeticViolatedPerc;
    int nS, nE;

    int nMinWindowWithCutEnergeticViolated;
    double nMinWindowWithCutEnergeticViolatedPerc;
    int nmS, nmE;

    int nJumps;
    //    int *maxShiftLift;
};


Results *Res_create( )
{
    Results *res;
    ALLOCATE_INI( res, Results );

    res->nRounds = 0;
    int nCuts = LP_CUT_NEW;
    ALLOCATE_VECTOR(res->tCutsTotal,VecDbl*,nCuts);
    ALLOCATE_VECTOR(res->nCutsTotal,VecInt*,nCuts);
    ALLOCATE_VECTOR(res->nElementsCuts,VecInt*,nCuts);
    ALLOCATE_VECTOR(res->nMaxElementsCuts,VecInt*,nCuts);
    ALLOCATE_VECTOR(res->nMinElementsCuts,VecInt*,nCuts);
    ALLOCATE_VECTOR(res->idxStartCutsMaxViol,VecInt*,nCuts);
    ALLOCATE_VECTOR(res->maxViolCutIdx,VecInt*,nCuts);
    ALLOCATE_VECTOR(res->maxViolCutCoef,VecDbl*,nCuts);
    ALLOCATE_VECTOR(res->nElementsMaxViol,VecInt*,nCuts);
    ALLOCATE_VECTOR(res->idxStartCutsMinViol,VecInt*,nCuts);
    ALLOCATE_VECTOR(res->minViolCutIdx,VecInt*,nCuts);
    ALLOCATE_VECTOR(res->minViolCutCoef,VecDbl*,nCuts);
    ALLOCATE_VECTOR(res->nElementsMinViol,VecInt*,nCuts);
    ALLOCATE_VECTOR(res->maxViol,VecDbl*,nCuts);
    ALLOCATE_VECTOR(res->minViol,VecDbl*,nCuts);
    ALLOCATE_VECTOR(res->sumViol,VecDbl*,nCuts);
    for(int i = 0; i< nCuts; i++) {
        res->tCutsTotal[i] = VDbl_create();
        res->nCutsTotal[i] = VInt_create();
        res->nElementsCuts[i] = VInt_create();
        res->nMaxElementsCuts[i] = VInt_create();
        res->nMinElementsCuts[i] = VInt_create();
        res->idxStartCutsMaxViol[i] = VInt_create();
        res->maxViolCutIdx[i] = VInt_create();
        res->maxViolCutCoef[i] = VDbl_create();
        res->nElementsMaxViol[i] = VInt_create();
        res->idxStartCutsMinViol[i] = VInt_create();
        res->minViolCutIdx[i] = VInt_create();
        res->minViolCutCoef[i] = VDbl_create();
        res->nElementsMinViol[i] = VInt_create();
        res->maxViol[i] = VDbl_create();
        res->minViol[i] = VDbl_create();
        res->sumViol[i] = VDbl_create();
    }

    res->nSumAllVarWithConf = VInt_create();
    res->nSumAllElementsConf = VInt_create();
    res->nMaxElementsConf = VInt_create();
    res->usoUR= VInt_create();
    res->usoURN= VInt_create();
    res->usoUM= VInt_create();
    res->usoUCI= VInt_create();
    res->usoUCP= VInt_create();

    res->totalUsoUR = 0;
    res->totalUsoURN = 0;
    res->totalUsoUM = 0;
    res->totalUsoUCI = 0;
    res->totalUsoUCP = 0;
    res->totalTimeSep = 0;
    res->totalTimeCGraph =0;
    res->nJumps = 0;
    res->nMaxWindowWithCutEnergeticViolated = 0;
    res->nMaxWindowWithCutEnergeticViolatedPerc = 0.0;
    res->nS = 0;
    res->nE = 0;
    res->nMinWindowWithCutEnergeticViolated = INT_MAX;
    res->nMinWindowWithCutEnergeticViolatedPerc = 0.0;
    res->nmS = 0;
    res->nmE = 0;
    return res;
}

void Res_initRound(Results *res)
{

    for(int i = 0; i< LP_CUT_NEW; i++) {
        VDbl_pushBack(res->tCutsTotal[i],0.0);
        VInt_pushBack(res->nCutsTotal[i],0);
        VInt_pushBack(res->nElementsCuts[i],0);
        VInt_pushBack(res->nMaxElementsCuts[i],0);
        VInt_pushBack(res->nMinElementsCuts[i],INT_MAX);
        VInt_pushBack(res->idxStartCutsMaxViol[i],0);
        VInt_pushBack(res->nElementsMaxViol[i],0);
        VInt_pushBack(res->idxStartCutsMinViol[i],0);
        VInt_pushBack(res->nElementsMinViol[i],0);
        VDbl_pushBack(res->maxViol[i],0.0);
        VDbl_pushBack(res->minViol[i],INT_MAX);
        VDbl_pushBack(res->sumViol[i],0.0);
    }
    VInt_pushBack(res->usoUR,0.0);
    VInt_pushBack(res->usoURN,0.0);
    VInt_pushBack(res->usoUM,0.0);
    VInt_pushBack(res->usoUCI,0.0);
    VInt_pushBack(res->usoUCP,0.0);

}

int Res_getRound(Results *res)
{
    assert(res);
    return res->nRounds;
}

void Res_setRound(Results *res, int value)
{
    assert(res);
    res->nRounds += value;
}


void Res_setTotalTimeSep(Results *res, double value)
{
    assert(res);
    res->totalTimeSep = value;
}

void Res_setTotalTimeOpt(Results *res, double value)
{
    assert(res);
    res->totalTimeOpt = value;
}

void Res_setTotalTimeRemAdd(Results *res, double value)
{
    assert(res);
    res->totalTimeRemAdd= value;
}

void Res_setTotalTimeCGraph(Results *res, double value)
{
    assert(res);
    res->totalTimeCGraph = value;
}


double Res_getTotalTimeCGraph(Results *res)
{
    assert(res);
    return res->totalTimeCGraph;
}


double Res_getTotalTimeRemAdd(Results *res)
{
    assert(res);
    return res->totalTimeRemAdd;
}

double Res_getTotalTimeSep(Results *res)
{
    assert(res);
    return res->totalTimeSep;
}

double Res_getTotalTimeOpt(Results *res)
{
    assert(res);
    return res->totalTimeOpt;
}
void Res_setUsoUR(Results *res, int nround, int value)
{
    assert(res);
    int newvalue = VInt_get(res->usoUR,nround)+value;
    VInt_set(res->usoUR, nround, newvalue );
}
void Res_setUsoURN(Results *res, int nround, int value)
{
    assert(res);
    int newvalue = VInt_get(res->usoURN,nround)+value;
    VInt_set(res->usoURN, nround, newvalue );
}
void Res_setUsoUM(Results *res, int nround, int value)
{
    assert(res);
    int newvalue = VInt_get(res->usoUM,nround)+value;
    VInt_set(res->usoUM, nround, newvalue );
}
void Res_setUsoUCI(Results *res, int nround, int value)
{
    assert(res);
    int newvalue = VInt_get(res->usoUCI,nround)+value;
    VInt_set(res->usoUCI, nround, newvalue );
}
void Res_setUsoUCP(Results *res, int nround, int value)
{
    assert(res);
    int newvalue = VInt_get(res->usoUCP,nround)+value;
    VInt_set(res->usoUCP, nround, newvalue );
}
void Res_setTotalUsoUR(Results *res, int value)
{
    assert(res);
    res->totalUsoUR +=value;
}
void Res_setTotalUsoURN(Results *res, int value)
{
    assert(res);
    res->totalUsoURN +=value;
}


void Res_setNJumps(Results *res, int value)
{
    assert(res);
    res->nJumps +=value;
}

int Res_getNJumps(Results *res)
{
    assert(res);
    return res->nJumps;
}

void Res_setTotalUsoUM(Results *res, int value)
{
    assert(res);
    res->totalUsoUM +=value;
}
void Res_setTotalUsoUCI(Results *res, int value)
{
    assert(res);
    res->totalUsoUCI +=value;
}
void Res_setTotalUsoUCP(Results *res, int value)
{
    assert(res);
    res->totalUsoUCP += value;
}

int Res_getUsoUR(Results *res, int nround)
{
    assert(res);
    return VInt_get(res->usoUR,nround);
}
int Res_getUsoURN(Results *res, int nround)
{
    assert(res);
    return VInt_get(res->usoURN,nround);
}
int Res_getUsoUM(Results *res, int nround)
{
    assert(res);
    return VInt_get(res->usoUM,nround);
}
int Res_getUsoUCI(Results *res, int nround)
{
    assert(res);
    return VInt_get(res->usoUCI,nround);
}
int Res_getUsoUCP(Results *res, int nround)
{
    assert(res);
    return VInt_get(res->usoUCP,nround);
}
int Res_getTotalUsoUR(Results *res)
{
    assert(res);
    return res->totalUsoUR;
}
int Res_getTotalUsoURN(Results *res)
{
    assert(res);
    return res->totalUsoURN;
}
int Res_getTotalUsoUM(Results *res)
{
    assert(res);
    return res->totalUsoUM;
}
int Res_getTotalUsoUCI(Results *res)
{
    assert(res);
    return res->totalUsoUCI;
}
int Res_getTotalUsoUCP(Results *res)
{
    assert(res);
    return res->totalUsoUCP;
}

int Res_getNMaxElementsConf(Results *res, int cut)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    return VInt_get(res->nMaxElementsConf,cut);
}

int Res_getMaxViolCutA(Results *res, int nround, int cut, int i)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    int idx = VInt_get(res->idxStartCutsMaxViol[cut],nround);
    return VInt_get(res->maxViolCutIdx[cut],idx+i);
}

double Res_getMaxViolCutB(Results *res, int nround, int cut, int i)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    int idx = VInt_get(res->idxStartCutsMaxViol[cut],nround);
    return VDbl_get(res->maxViolCutCoef[cut],idx+i);

}

int Res_getMinViolCutA(Results *res, int nround, int cut, int i)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    int idx = VInt_get(res->idxStartCutsMinViol[cut],nround);
    return VInt_get(res->minViolCutIdx[cut],idx+i);
}

double Res_getMinViolCutB(Results *res, int nround, int cut, int i)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    int idx = VInt_get(res->idxStartCutsMinViol[cut],nround);
    return VDbl_get(res->minViolCutCoef[cut],idx+i);
}

int Res_lastIdxMaxViolCut(Results *res, int nround, int cut)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    int idx = VInt_size(res->maxViolCutIdx[cut]);
    return idx;

}

void Res_setIdxStartCutsMaxViol(Results *res, int cut)
{

    assert(res);
    cut = cut-LP_CUT_TYPES;
    int idx = VInt_size(res->maxViolCutIdx[cut]);
    VInt_pushBack(res->idxStartCutsMaxViol[cut], idx);
}


void Res_setIdxStartCutsMinViol(Results *res, int cut)
{

    assert(res);
    cut = cut-LP_CUT_TYPES;
    int idx = VInt_size(res->minViolCutIdx[cut]);
    VInt_pushBack(res->idxStartCutsMinViol[cut], idx);
}


void Res_setMaxViolCut(Results *res, int nround, int cut, int i, int a, double b)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    // int idx = VInt_get(res->idxStartCutsMaxViol[cut],nround);
    // if(i>idx){
    VInt_pushBack(res->maxViolCutIdx[cut],a);
    VDbl_pushBack(res->maxViolCutCoef[cut],b);
    // }else{
    //   VInt_set(res->maxViolCutIdx[cut],idx+i, a);
    //  VDbl_set(res->maxViolCutCoef[cut],idx+i, b);
    //}
}

void Res_setMinViolCut(Results *res, int nround, int cut, int i, int a, double b)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;

    VInt_pushBack(res->minViolCutIdx[cut],a);
    VDbl_pushBack(res->minViolCutCoef[cut],b);

    //    int idx = VInt_get(res->idxStartCutsMinViol[cut],nround);
    //  VInt_set(res->minViolCutIdx[cut],idx+i, a);
    //  VDbl_set(res->minViolCutCoef[cut],idx+i, b);
}

double Res_getSumViol(Results *res, int nround, int cut)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    return VDbl_get(res->sumViol[cut], nround);
}


int Res_getNElementsCuts(Results *res, int nround, int cut)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    return VInt_get(res->nElementsCuts[cut], nround);
}

int Res_getNElementsMinViol(Results *res, int nround, int cut)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    return VInt_get(res->nElementsMinViol[cut], nround);
}


int Res_getNElementsMaxViol(Results *res, int nround, int cut)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    return VInt_get(res->nElementsMaxViol[cut], nround);
}

int Res_getNMaxElementsCuts(Results *res, int nround, int cut)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    return VInt_get(res->nMaxElementsCuts[cut], nround);
}

int Res_getNCutsTotal(Results *res, int nround, int cut)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    return VInt_get(res->nCutsTotal[cut], nround);
}



int Res_getNMinElementsCuts(Results *res, int nround, int cut)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    return VInt_get(res->nMinElementsCuts[cut], nround);
}


double Res_getMaxViol(Results *res, int nround, int cut)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    return VDbl_get(res->maxViol[cut], nround);
}

double Res_getMinViol(Results *res, int nround, int cut)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    return VDbl_get(res->minViol[cut], nround);
}


void Res_setTCutsTotal(Results *res, int nround,int cut, double value)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    double newvalue = VDbl_get(res->tCutsTotal[cut],nround) + value;
    VDbl_set(res->tCutsTotal[cut], nround, newvalue );
}


void Res_setTW(Results *res, int nMaxTW, double nMaxTWP, int nS, int nE)
{
    assert(res);
    res->nMaxWindowWithCutEnergeticViolated = nMaxTW;
    res->nMaxWindowWithCutEnergeticViolatedPerc = nMaxTWP;
    res->nS = nS;
    res->nE = nE;
}


void Res_setMinTW(Results *res, int nMinTW, double nMinTWP, int nmS, int nmE)
{
    assert(res);
    res->nMinWindowWithCutEnergeticViolated = nMinTW;
    res->nMinWindowWithCutEnergeticViolatedPerc = nMinTWP;
    res->nmS = nmS;
    res->nmE = nmE;
}


void Res_setNElementsMaxViol(Results *res, int nround, int cut, int value)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    VInt_set(res->nElementsMaxViol[cut], nround, value );
}

void Res_setNElementsMinViol(Results *res, int nround, int cut, int value)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    VInt_set(res->nElementsMinViol[cut], nround, value );

}

void Res_setNCutsTotal(Results *res, int nround, int cut, int value)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    int newvalue = VInt_get(res->nCutsTotal[cut],nround) + value;
    VInt_set(res->nCutsTotal[cut], nround, newvalue );
}


void Res_setNElementsCuts(Results *res, int nround, int cut, int value)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    int newvalue = VInt_get(res->nElementsCuts[cut],nround) + value;
    VInt_set(res->nElementsCuts[cut], nround, newvalue );
}


void Res_setNMaxElementsCuts(Results *res, int nround, int cut, int value)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    VInt_set(res->nMaxElementsCuts[cut], nround, value );
}


void Res_setNMinElementsCuts(Results *res, int nround, int cut, int value)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    VInt_set(res->nMinElementsCuts[cut], nround, value );

}


void Res_setMaxViol(Results *res, int nround, int cut, double value)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    VDbl_set(res->maxViol[cut], nround, value );
}


void Res_setMinViol(Results *res, int nround, int cut, double value)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    VDbl_set(res->minViol[cut], nround, value );
}


void Res_setSumViol(Results *res, int nround, int cut, double value)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    double newvalue = VDbl_get(res->sumViol[cut],nround) + value;
    VDbl_set(res->sumViol[cut], nround, newvalue );
}

void Res_setNSumAllVarWithConf(Results *res, int cut, int value)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    VInt_set(res->nSumAllVarWithConf, cut, value );
}

void Res_setNSumAllElementsConf(Results *res, int cut, int value)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    VInt_set(res->nSumAllElementsConf, cut, value );
}

void Res_setNMaxElementsConf(Results *res, int cut, int value)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    VInt_set(res->nMaxElementsConf, cut, value );
}

double Res_getTCutsTotal(Results *res, int nround, int cut)
{
    assert(res);
    cut = cut-LP_CUT_TYPES;
    return VDbl_get(res->tCutsTotal[cut], nround);
}


int Res_getTW(Results *res)
{
    assert(res);
    return  res->nMaxWindowWithCutEnergeticViolated;
}

double Res_getTWP(Results *res)
{
    assert(res);
    return  res->nMaxWindowWithCutEnergeticViolatedPerc;
}

int Res_getS(Results *res)
{
    assert(res);
    return  res->nS;
}


int Res_getE(Results *res)
{
    assert(res);
    return  res->nE;
}


int Res_getMinTW(Results *res)
{
    assert(res);
    return  res->nMinWindowWithCutEnergeticViolated;
}

double Res_getMinTWP(Results *res)
{
    assert(res);
    return  res->nMinWindowWithCutEnergeticViolatedPerc;
}

int Res_getMinS(Results *res)
{
    assert(res);
    return  res->nmS;
}


int Res_getMinE(Results *res)
{
    assert(res);
    return  res->nmE;
}

void Res_free( Results **_res )
{
    Results *res = *_res;

    for(int i = 0; i< LP_CUT_NEW; i++) {
        VDbl_free(&res->tCutsTotal[i]);
        VInt_free(&res->nCutsTotal[i]);
        VInt_free(&res->nElementsCuts[i]);
        VInt_free(&res->nMaxElementsCuts[i]);
        VInt_free(&res->nMinElementsCuts[i] );
        VInt_free(&res->idxStartCutsMaxViol[i]);
        VInt_free(&res->maxViolCutIdx[i]);
        VDbl_free(&res->maxViolCutCoef[i]);
        VInt_free(&res->nElementsMaxViol[i]);
        VInt_free(&res->idxStartCutsMinViol[i]);
        VInt_free(&res->minViolCutIdx[i]);
        VDbl_free(&res->minViolCutCoef[i] );
        VInt_free(&res->nElementsMinViol[i]);
        VDbl_free(&res->maxViol[i]);
        VDbl_free(&res->minViol[i]);
        VDbl_free(&res->sumViol[i]);
    }

    free(res->tCutsTotal);
    free(res->nCutsTotal);
    free(res->nElementsCuts);
    free(res->nMaxElementsCuts);
    free(res->nMinElementsCuts);
    free(res->idxStartCutsMaxViol);
    free(res->maxViolCutIdx);
    free(res->maxViolCutCoef);
    free(res->nElementsMaxViol);
    free(res->idxStartCutsMinViol);
    free(res->minViolCutIdx);
    free(res->minViolCutCoef );
    free(res->nElementsMinViol);
    free(res->maxViol);
    free(res->minViol);
    free(res->sumViol);

    VInt_free(&res->nSumAllVarWithConf);
    VInt_free(&res->nSumAllElementsConf);
    VInt_free(&res->nMaxElementsConf);
    VInt_free(&res->usoUR);
    VInt_free(&res->usoURN);
    VInt_free(&res->usoUM);
    VInt_free(&res->usoUCI);
    VInt_free(&res->usoUCP);

    free( res );
    *_res = NULL;
}


void Res_writeResults(Parameters *par, LinearProgram *mip, const Solution* sol, Results *res, double bo, FILE *fp, char **argv, int argc, double startT, int sumtpd, int sumtms)
{

    assert(par != NULL);
    assert(mip != NULL);
    assert(sol != NULL);

    /*Writing solution*/
    int cg = Par_getCutCG( par );
    int cggpur2 = Par_getCutCGGPUR2( par );
    int cp = Par_getCutPrec( par );
    int cr = Par_getCutRR( par );
    int cc = Par_getCutCLIQUE( par );
    int co = Par_getCutODDHOLES( par );
    int cgcpu = Par_getCutCGCPU( par );
    double slk = Par_getSlack( par );
    double mc = Par_getMaxCut(par);
    int mi = Par_getMaxInstant(par);
    int ju = Par_getJump(par);
   // double gap = lp_get_gap(mip);

    char filelp[256]=" ";
    sprintf( filelp, "%s_cg_%d_cggpur2_%d_cp_%d_cr_%d_cc_%d_co_%d_cgcpu_%d_l_%d_mc_%f_rc_%f_slk_%f_mi_%d_jmp_%d.lp", Par_getName(par), cg,cggpur2,cp,cr,cc,co,cgcpu, Par_getLifting(par),  mc, Par_getMaxReducedCost(par), slk,mi,ju);
    printf("\n%s_cg_%d_cggpur2_%d_cp_%d_cr_%d_cc_%d_co_%d_cgcpu_%d_l_%d_mc_%f_rc_%f_slk_%f_mi_%d_jmp_%d\n",  Par_getName(par), cg,cggpur2,cp,cr,cc,co,cgcpu, Par_getLifting(par),  mc, Par_getMaxReducedCost(par), slk,mi,ju);
    lp_write_lp(mip,filelp);

    int nCutClique = 0, nCutODDHOLES = 0, nCutRR=0, nCutPrec =0,  nCutCG = 0, nCutCGGPU = 0, nCutCGGPUR2 = 0, nCutCGCPU = 0; //nCutJS =0,
    double tCutClique = 0, tCutODDHOLES = 0, tCutRR = 0, tCutPrec = 0, tCutCG = 0, tCutCGGPU = 0,tCutCGGPUR2 = 0,tCutCGCPU = 0; // tCutJS = 0,
    int maxCutClique = 0, maxCutODDHOLES = 0, maxCutRR = 0, maxCutPrec =0, maxCutCGCPU =0,  maxCutCG =0, maxCutCGGPU =0, maxCutCGGPUR2 =0, maxCut =0, minCut = INT_MAX, minCutClique = INT_MAX, minCutODDHOLES = INT_MAX, minCutRR = INT_MAX, minCutCG = INT_MAX, minCutCGGPU = INT_MAX, minCutCGGPUR2 = INT_MAX,  minCutPrec = INT_MAX, minCutCGCPU = INT_MAX; //maxCutJS =0,  minCutJS = INT_MAX,
    //  int maxElementsConf = 0;
    int contUR = 0, contURN = 0, contUM = 0, contUCI = 0, contUCP=0, nJumps=0;
    if(VERBOSE) {
        for(int r = 0; r < Res_getRound(res); r++) {
            nCutCGCPU += Res_getNCutsTotal(res,r,LPC_CGCPU);
            nCutClique += Res_getNCutsTotal(res,r,LPC_CLIQUE);
            nCutODDHOLES += Res_getNCutsTotal(res,r,LPC_ODDHOLES);
            nCutRR += Res_getNCutsTotal(res,r,LPC_RR);
            nCutPrec += Res_getNCutsTotal(res,r,LPC_PREC);
            nCutCG += Res_getNCutsTotal(res,r,LPC_CG);
//            nCutJS += Res_getNCutsTotal(res,r,LPC_JS);
            nCutCGGPU += Res_getNCutsTotal(res,r,LPC_CGGPU);
            nCutCGGPUR2 += Res_getNCutsTotal(res,r,LPC_CGGPUR2);
            tCutCGCPU += Res_getTCutsTotal(res,r,LPC_CGCPU);
            tCutClique += Res_getTCutsTotal(res,r,LPC_CLIQUE);
            tCutODDHOLES += Res_getTCutsTotal(res,r,LPC_ODDHOLES);
            tCutCG += Res_getTCutsTotal(res,r,LPC_CG);
            tCutCGGPU += Res_getTCutsTotal(res,r,LPC_CGGPU);
            tCutCGGPUR2 += Res_getTCutsTotal(res,r,LPC_CGGPUR2);
            tCutRR += Res_getTCutsTotal(res,r,LPC_RR);
            tCutPrec += Res_getTCutsTotal(res,r,LPC_PREC);
//            tCutJS += Res_getTCutsTotal(res,r,LPC_JS);
            maxCut = Res_getNMaxElementsCuts(res,r,LPC_CLIQUE);
            if(maxCutClique < maxCut)
                maxCutClique =maxCut;
            maxCut = Res_getNMaxElementsCuts(res,r,LPC_ODDHOLES);
            if(maxCutODDHOLES < maxCut)
                maxCutODDHOLES =maxCut;
            maxCut = Res_getNMaxElementsCuts(res,r,LPC_RR);
            if(maxCutRR < maxCut)
                maxCutRR =maxCut;
            maxCut = Res_getNMaxElementsCuts(res,r,LPC_PREC);
            if(maxCutPrec < maxCut)
                maxCutPrec =maxCut;
/*           maxCut = Res_getNMaxElementsCuts(res,r,LPC_JS);
            if(maxCutJS < maxCut)
                maxCutJS =maxCut;*/
            maxCut = Res_getNMaxElementsCuts(res,r,LPC_CG);
            if(maxCutCG < maxCut)
                maxCutCG =maxCut;
            maxCut = Res_getNMaxElementsCuts(res,r,LPC_CGGPU);
            if(maxCutCGGPU < maxCut)
                maxCutCGGPU =maxCut;
            maxCut = Res_getNMaxElementsCuts(res,r,LPC_CGGPUR2);
            if(maxCutCGGPUR2 < maxCut)
                maxCutCGGPUR2 =maxCut;
            maxCut = Res_getNMaxElementsCuts(res,r,LPC_CGCPU);
            if(maxCutCGCPU < maxCut)
                maxCutCGCPU =maxCut;
            minCut = Res_getNMinElementsCuts(res,r,LPC_CGCPU);
            if(minCutCGCPU > minCut)
                minCutCGCPU =minCut;
            minCut = Res_getNMinElementsCuts(res,r,LPC_CLIQUE);
            if(minCutClique > minCut)
                minCutClique =minCut;
            minCut = Res_getNMinElementsCuts(res,r,LPC_ODDHOLES);
            if(minCutODDHOLES > minCut)
                minCutODDHOLES =minCut;
            minCut = Res_getNMinElementsCuts(res,r,LPC_RR);
            if(minCutRR > minCut)
                minCutRR =minCut;
            minCut = Res_getNMinElementsCuts(res,r,LPC_PREC);
            if(minCutPrec > minCut)
                minCutPrec =minCut;
/*            minCut = Res_getNMinElementsCuts(res,r,LPC_JS);
            if(minCutJS > minCut)
                minCutJS =minCut; */
            minCut = Res_getNMinElementsCuts(res,r,LPC_CG);
            if(minCutCG > minCut)
                minCutCG =minCut;
            minCut = Res_getNMinElementsCuts(res,r,LPC_CGGPU);
            if(minCutCGGPU > minCut)
                minCutCGGPU =minCut;
            minCut = Res_getNMinElementsCuts(res,r,LPC_CGGPUR2);
            if(minCutCGGPUR2 > minCut)
                minCutCGGPUR2 =minCut;

            //  if(maxElementsConf < Res_getNMaxElementsConf(res,LPC_CLIQUE))
            //    maxElementsConf = Res_getNMaxElementsConf(res,LPC_CLIQUE);
        }
        contUR=Res_getTotalUsoUR(res), contURN=Res_getTotalUsoURN(res), contUM=Res_getTotalUsoUM(res), contUCI=Res_getTotalUsoUCI(res), contUCP = Res_getTotalUsoUCP(res), nJumps = Res_getNJumps(res);
    }

    double timesep = Res_getTotalTimeSep(res); //tCutClique+tCutPrec+tCutRR+tCutCG+tCutCGGPU+tCutCGGPUR2;
    double bc = Inst_getSumTPD(par->inst) == 0 ? 0 : (double) bo/Inst_getSumTPD(par->inst);
    double timeopt = Res_getTotalTimeOpt(res);// ((omp_get_wtime()-startT)-timesep);
    double roundedbo =  ceil(bo);
    double bcrounded = Inst_getSumTPD(par->inst) == 0 ? 0 : (double) roundedbo/Inst_getSumTPD(par->inst);
    double timeremadd =  Res_getTotalTimeRemAdd(res);
    double timecgraph = Res_getTotalTimeCGraph(res);
    // fprintf( fp, " %s ; %ld ; %f ; %d ; ",  argv[3], Sol_getCost(MipC_getSol(mipC)),  (omp_get_wtime()-startT), sumtpd);
    //   printf("\n roundUR %d, contUR %d, roundURN %d, contURN %d, roundUM %d, contUM %d, roundUCI %d, contUCI %d, roundUCP %d, contUCP %d \n", Res_getUsoUR(mipC->res,nround), Res_getTotalUsoUR(mipC->res), Res_getUsoURN(mipC->res,nround), Res_getTotalUsoURN(mipC->res), Res_getUsoUM(mipC->res,nround), Res_getTotalUsoUM(mipC->res),Res_getUsoUCI(mipC->res,nround), Res_getTotalUsoUCI(mipC->res), Res_getUsoUCP(mipC->res,nround), Res_getTotalUsoUCP(mipC->res));


    fprintf( fp, " %s ; %f ; %d ; %d ; %ld ; %ld ; ",  argv[3], (omp_get_wtime()-startT), sumtpd, sumtms, Sol_getTPD(sol),Sol_getTMS(sol));

    if(VERBOSE) {
        for(int n = 0; n< argc ; n++) {
            /*if (strcmp(argv[n],"-mrcpsp") == 0) {
                //fprintf( fp, " %s ; ", argv[n]);
                n++;
                fprintf( fp, "%d ; ", atoi(argv[n]));
                continue;
            }
            if (strcmp(argv[n],"-rcpsp") == 0) {
                //fprintf( fp, " %s ; ", argv[n]);
                n++;
                fprintf( fp, "%d ; ", atoi(argv[n]));
                continue;
            }
            if (strcmp(argv[n],"-mmrcmpsp") == 0) {
                //fprintf( fp, " %s ; ", argv[n]);
                n++;
                fprintf( fp, "%d ; ", atoi(argv[n]));
                continue;
            }
            if (strcmp(argv[n],"-continuous") == 0) {
                    // fprintf( fp, " %s ", argv[n]);
                    n++;
                    fprintf( fp, "%d ; ", atoi(argv[n]));
                    continue;
            }*/
            if (strcmp(argv[n],"-cutRR") == 0) {
                // fprintf( fp, " %s ", argv[n]);
                n++;
                // fprintf( fp, "%d ; ", atoi(argv[n]));
                fprintf( fp, "%d ; %f ; %d ; %d ; ",  nCutRR, tCutRR, minCutRR, maxCutRR);
                continue;
            }
            if (strcmp(argv[n],"-cutPREC") == 0) {
                //fprintf( fp, " %s ", argv[n]);
                n++;
                //fprintf( fp, "%d ; ", atoi(argv[n]));
                fprintf( fp, "%d ; %f ; %d ; %d ;", nCutPrec,tCutPrec, minCutPrec, maxCutPrec );
                continue;
            }
         /*    if (strcmp(argv[n],"-cutJS") == 0) {
                //fprintf( fp, " %s ", argv[n]);
                n++;
                //fprintf( fp, "%d ; ", atoi(argv[n]));
                fprintf( fp, "%d ; %f ; %d ; %d ;", nCutJS,tCutJS, minCutJS, maxCutJS );
                continue;
            }*/
            if (strcmp(argv[n],"-cutODDHOLES") == 0) {
                // fprintf( fp, " %s ", argv[n]);
                n++;
                // fprintf( fp, "%d ; ", atoi(argv[n]));
                fprintf( fp, "%d ; %f ; %d ; %d ; ", nCutODDHOLES, tCutODDHOLES, minCutODDHOLES, maxCutODDHOLES);//, maxElementsConf);
                continue;
            }
             if (strcmp(argv[n],"-cutCGCPU") == 0) {
                // fprintf( fp, " %s ", argv[n]);
                n++;
                // fprintf( fp, "%d ; ", atoi(argv[n]));
                fprintf( fp, "%d ; %f ; %d ; %d ; ", nCutCGCPU, tCutCGCPU, minCutCGCPU, maxCutCGCPU );//, maxElementsConf);
                continue;
            }
            if (strcmp(argv[n],"-cutCLIQUE") == 0) {
                // fprintf( fp, " %s ", argv[n]);
                n++;
                // fprintf( fp, "%d ; ", atoi(argv[n]));
                fprintf( fp, "%d ; %f ; %d ; %d ; ", nCutClique, tCutClique, minCutClique, maxCutClique);//, maxElementsConf);
                continue;
            }
            if (strcmp(argv[n],"-cutCG") == 0) {
                // fprintf( fp, " %s ", argv[n]);
                n++;
                //fprintf( fp, "%d ; ", atoi(argv[n]));
                fprintf( fp, "%d ; %f ; %d ; %d ;  ", nCutCG, tCutCG, minCutCG, maxCutCG);
                fprintf( fp, "%d ; %d ; %d ; %d ; %d ; %d ; ", contUR, contURN, contUM, contUCI, contUCP, nJumps);
                continue;
            }
            if (strcmp(argv[n],"-cutCGGPU") == 0) {
                // fprintf( fp, " %s ", argv[n]);
                n++;
                //fprintf( fp, "%d ; ", atoi(argv[n]));
                fprintf( fp, "%d ; %f ; %d ; %d ;  ", nCutCGGPU, tCutCGGPU, minCutCGGPU, maxCutCGGPU);
                //                fprintf( fp, "%d ; %d ; %d ; %d ; %d ; %d; ", contUR, contURN, contUM, contUCI, contUCP, nJumps);
                continue;
            }
            if (strcmp(argv[n],"-cutCGGPUR2") == 0) {
                // fprintf( fp, " %s ", argv[n]);
                n++;
                //fprintf( fp, "%d ; ", atoi(argv[n]));
                fprintf( fp, "%d ; %f ; %d ; %d ;  ", nCutCGGPUR2, tCutCGGPUR2, minCutCGGPUR2, maxCutCGGPUR2);
                //                fprintf( fp, "%d ; %d ; %d ; %d ; %d ; %d; ", contUR, contURN, contUM, contUCI, contUCP, nJumps);
                continue;
            }
            /*  if (strcmp(argv[n],"-maxinstant") == 0) {
                  n++;
                  fprintf( fp, "%d ; ", atoi(argv[n]));
                  continue;
              }
              if (strcmp(argv[n],"-jump") == 0) {
                  n++;
                  fprintf( fp, "%d ; ", atoi(argv[n]));
                  continue;
              }
              if (strcmp(argv[n],"-cutDefaultCBC") == 0) {
                  //fprintf( fp, " %s ", argv[n]);
                  n++;
                  fprintf( fp, "%d ; ", atoi(argv[n]));
                  continue;
              }
              if (strcmp(argv[n],"-gomory") == 0) {
                  //fprintf( fp, " %s ", argv[n]);
                  n++;
                  fprintf( fp, "%d ; ", atoi(argv[n]));
                  continue;
              }
              if (strcmp(argv[n],"-reduce") == 0) {
                  //fprintf( fp, " %s ", argv[n]);
                  n++;
                  fprintf( fp, "%d ; ", atoi(argv[n]));
                  continue;
              }
              if (strcmp(argv[n],"-mir") == 0) {
                  //fprintf( fp, " %s ", argv[n]);
                  n++;
                  fprintf( fp, "%d ; ", atoi(argv[n]));
                  continue;
              }
              if (strcmp(argv[n],"-twomir") == 0) {
                  //fprintf( fp, " %s ", argv[n]);
                  n++;
                  fprintf( fp, "%d ; ", atoi(argv[n]));
                  continue;
              }
              if (strcmp(argv[n],"-landp") == 0) {
                  //fprintf( fp, " %s ", argv[n]);
                  n++;
                  fprintf( fp, "%d ; ", atoi(argv[n]));
                  continue;
              }
              if (strcmp(argv[n],"-zerohalf") == 0) {
                  //fprintf( fp, " %s ", argv[n]);
                  n++;
                  fprintf( fp, "%d ; ", atoi(argv[n]));
                  continue;
              }
              if (strcmp(argv[n],"-knapsack") == 0) {
                  //fprintf( fp, " %s ", argv[n]);
                  n++;
                  fprintf( fp, "%d ; ", atoi(argv[n]));
                  continue;
              }

              if (strcmp(argv[n],"-flow") == 0) {
                  //fprintf( fp, " %s ", argv[n]);
                  n++;
                  fprintf( fp, "%d ; ", atoi(argv[n]));
                  continue;
              }
              if (strcmp(argv[n],"-maxNode") == 0) {
                  //fprintf( fp, " %s ", argv[n]);
                  n++;
                  fprintf( fp, "%d ; ", atoi(argv[n]));
                  continue;
              }
              if (strcmp(argv[n],"-maxReducedCost") == 0) {
                  // fprintf( fp, " %s ", argv[n]);
                  n++;
                  fprintf( fp, "%f ; ", atof(argv[n]));
                  continue;
              }
              if (strcmp(argv[n],"-slack") == 0) {
                  // fprintf( fp, " %s ", argv[n]);
                  n++;
                  fprintf( fp, "%f ; ", atof(argv[n]));
                  continue;
              }
              if (strcmp(argv[n],"-lifting") == 0) {
                  // fprintf( fp, " %s ", argv[n]);
                  n++;
                  fprintf( fp, "%d ; ", atoi(argv[n]));
                  continue;
              }
              if (strcmp(argv[n],"-maxcuts") == 0) {
                  // fprintf( fp, " %s ", argv[n]);
                  n++;
                  double value = atoi(argv[n]);
                  if(value == 1)
                      mipC->maxcuts =  Inst_nJobs(mipC->inst);
                  else if(value == 2)
                      mipC->maxcuts =  lp_cols(mipC->mip);
                  else
                      mipC->maxcuts =  value;
                  fprintf( fp, "%f ; ", value);
                  continue;
              }*/
        }
        fprintf( fp, "%f ; %f ; %f ; %f ; %d ; %g ; %g ; %g ; %g ;\n", timesep,  timeopt, timeremadd, timecgraph, Res_getRound(res), bo, bc, roundedbo, bcrounded);
    } else
        fprintf( fp, "%g ; %g ; %g ; %g ;\n", bo, bc, roundedbo, bcrounded);
}


