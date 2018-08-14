/*
 * solutionGpu.h
 *
 *  Created on: 31/03/2017
 *      Author: danilo
*/

#ifndef SOLUTION_GPU_H_
#define SOLUTION_GPU_H_


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "cut_gpu.h"
#include "gpulib/types.h"
#include "lp.h"
#include <time.h>
#include <sys/time.h>

EXTERN_C_BEGIN



//typedef int TSCoefficients;
//typedef int TSRightSide;
//typedef int TSViolation;
typedef int TSMult;
typedef int TSConst;
typedef int TSPAux;
typedef short int TSGroup;
typedef unsigned char TSVetSol;

typedef struct {
    TSMult *SMult;
    TSConst *SConst;
    TSPAux *SPAux;
} solutionGpu;

/*Zero Half Struct*/
typedef struct{
    TSGroup sizeGroup;
    TSVetSol *VetSol;
}solutionZeroHalf;

typedef struct{
    long    tv_sec;         /* seconds */
    long    tv_usec;        /* and microseconds */
}timeval_t;

solutionZeroHalf* allocationStructSolutionZeroHalf(int sizeGroup, int nRuns);

solutionGpu* allocationStructSolution2(Cut_gpu *c, int numberMaxConst, int nRuns);

solutionGpu* allocationStructSolution1(Cut_gpu *c, int nRuns);

solutionGpu* createGPUsolution2(solutionGpu* h_solution, Cut_gpu* h_cut,int numberMaxConst, int nRuns);

solutionGpu* createGPUsolution1(solutionGpu* h_solution, Cut_gpu* h_cut, int nRuns);

solutionZeroHalf* createGPUsolutionZeroHalf(solutionZeroHalf* h_solution_zero, int sizeGroup, int nRuns);

Cut_gpu* createCutsOfPhaseOne(Cut_gpu *h_cut, solutionGpu *h_solution, int nCuts, int precision, int nRuns);

Cut_gpu* createCutsOfPhaseTwo(Cut_gpu *h_cut, Cut_gpu_aux *cut_aux, solutionGpu *h_solution, int numberMaxConst, int nCuts, int precision, int nRuns);


Cut_gpu_aux* reallocCut(Cut_gpu *h_cut,Cut_gpu_aux *h_cut_aux);

Cut_gpu_aux* reallocCutR2(Cut_gpu *h_cut,Cut_gpu_aux *h_cut_aux, int *cont);

void bubble_sort(int *vetor,int *pos, int n);

int* returnOrdConstrainsNR(Cut_gpu *cut);

float* returnFolga(Cut_gpu *cut);

//void calcSetConstraint (int *setConstraint, int numberMaxConst,int numberConstrains, int *resR1, int *resNR1, int sizeR1, int sizeNR1, int *Similar, float *Folga, int nRuns);
void calcSetConstraint (int *setConstraint, int *pos_R1, int numberMaxConst,int numberConstrains, int *resR1, int *resNR1, int sizeR1, int sizeNR1, int *Similar, float *Folga, int nRuns, int szR );

void returnCutCG(Cut_gpu *ccg, CutCG* ccgr2, int precision);

void shuffle_Set(int *vec, int nSetConstrains, int n);

int returnK(long double fa_zero);

long double closestIntegerDistance(long double v);
long double ldRound( long double v);


int verifyDominanceCG(int *v1, int rhs1, int *v2, int rhs2, int sz);

Cut_gpu* createCutsStrongPhaseOne(Cut_gpu *h_cut, solutionGpu *h_solution, int nCuts, int precision, int nRuns, int nC_initial, int cont_ini);


Cut_gpu* createCutsStrongPhaseTwo(Cut_gpu *h_cut, solutionGpu *h_solution, int numberMaxConst, int nCuts, int precision, int nRuns, int nThreads,int nBlocks);

void DecimalToBinaryCpu(int number,unsigned char v[], int sizegroup);

Cut_gpu* createCutsStrongZeroHalf(Cut_gpu *h_cut, int *h_solution, int sizeGroup, int nBlocks, int nThreads, int precision, int nCuts, int cont_ini, int nC_initial);

Cut_gpu* complementCutPhase1(Cut_gpu *h_cut, int k);

Cut_gpu* complementCutPhase1_RR(Cut_gpu *h_cut, Cut_gpu_aux *h_cut_aux, int nRR_cggpu );

int fat(int n);

int combinacao_Phase1(int n, int p);

EXTERN_C_END

#endif
