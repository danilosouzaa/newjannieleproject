/*
 * gSolutionGpu.cuh
 *
 *  Created on: 31/03/2017
 *      Author: danilo
 */

#ifndef GSOLUTION_GPU_CUH_
#define GSOLUTION_GPU_CUH_

#include "gpulib/gpu.cuh"
#include <curand.h>
#include <curand_kernel.h>

extern "C" {

#include "cut_gpu.h"
#include "solutionGpu.h"
//#include "configGpu.h"

}

#include "lp.h"


__global__ void runGPUR1(Cut_gpu *d_cut, solutionGpu *d_solution, unsigned int *seed, curandState_t* states, int nThreads, int precision);

__global__ void runGPUR1_aleatory(Cut_gpu *d_cut, solutionGpu *d_solution, unsigned int *seed, curandState_t* states, int nThreads, int precision, int maxDenominator);

__global__ void runGPUR2(Cut_gpu *d_cut, solutionGpu *d_solution, unsigned int *seed, curandState_t* states, int numberMaxConst, int *setConstraint,int nThreads, int precision, int maxDenominator);

__device__ void DecimalToBinary(int number,unsigned char v[], int sizegroup);

__global__ void runGPUZerohHalf(Cut_gpu *d_cut, int *d_solution, int nThread,int sizeGroup, int nBlock, int precision);


#endif /* GSOLUTION_CUH_ */
