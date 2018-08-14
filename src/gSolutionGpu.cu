/*
 * gSolution.cu
 *
 *  Created on: 31/03/2017
 *      Author: danilo
 */
#include "gSolutionGpu.cuh"


solutionZeroHalf* createGPUsolutionZeroHalf(solutionZeroHalf* h_solution_zero, int sizeGroup, int nRuns)
{
    size_t size_solution_zero =  sizeof(solutionZeroHalf) +
                                 sizeof(TSVetSol)*(sizeGroup*nRuns);
    solutionZeroHalf* d_sol_zero;
    gpuMalloc((void**)&d_sol_zero,size_solution_zero);
    gpuMemset(d_sol_zero,0,size_solution_zero);
    h_solution_zero->VetSol = (TSVetSol*)(d_sol_zero+1);
    gpuMemcpy(d_sol_zero, h_solution_zero, size_solution_zero, cudaMemcpyHostToDevice);
    return d_sol_zero;
}


solutionGpu* createGPUsolution1(solutionGpu* h_solution, Cut_gpu* h_cut, int nRuns)
{
    size_t size_solution =  sizeof(solutionGpu) +
                            sizeof(TSMult)*(nRuns) +
                            sizeof(TSConst)*(nRuns) +
                            sizeof(TSPAux)*(nRuns);
    solutionGpu *d_sol;
    gpuMalloc((void**)&d_sol, size_solution);
    gpuMemset(d_sol,0,size_solution);
    h_solution->SMult = (TSMult*)(d_sol+1);
    h_solution->SConst= (TSConst*)(h_solution->SMult + (nRuns));
    h_solution->SPAux = (TSPAux*)(h_solution->SConst + (nRuns));
    gpuMemcpy(d_sol, h_solution, size_solution, cudaMemcpyHostToDevice);
    return d_sol;
}


solutionGpu* createGPUsolution2(solutionGpu* h_solution, Cut_gpu* h_cut,int numberMaxConst, int nRuns)
{

    size_t size_solution =  sizeof(solutionGpu) +
                            sizeof(TSMult)*(nRuns*4) +
                            sizeof(TSConst)*(numberMaxConst*nRuns) +
                            sizeof(TSPAux)*(nRuns);

    solutionGpu *d_sol;
    gpuMalloc((void**)&d_sol, size_solution);
    gpuMemset(d_sol,0,size_solution);
    h_solution->SMult = (TSMult*)(d_sol+1);
    h_solution->SConst= (TSConst*)(h_solution->SMult + (nRuns*4));
    h_solution->SPAux = (TSPAux*)(h_solution->SConst + (numberMaxConst*nRuns));
    gpuMemcpy(d_sol, h_solution, size_solution, cudaMemcpyHostToDevice);
    return d_sol;
}


__global__ void runGPUR1(Cut_gpu *d_cut, solutionGpu *d_solution, unsigned int *seed, curandState_t* states, int nThreads, int precision)
{
    int term = threadIdx.x + blockIdx.x*nThreads;
    __shared__ int *constraints;
    __shared__ int pos;
    curand_init(seed[term],term,0,&states[term]);
    int violation = 0,i,j;
    if(threadIdx.x == 0)
    {
        pos = 0;
        constraints = (int*)malloc(sizeof(int)*d_cut->numberConstrains);
        for(i=0; i<d_cut->numberConstrains; i++)
        {
            if((d_cut->typeConstraints[i] == RES_RR)  )
            {
                constraints[pos] = i;
                pos++;
            }
        }
    }
    __syncthreads();
    int res = constraints[threadIdx.x%pos];
    int sz_u = d_cut->ElementsConstraints[res+1] - d_cut->ElementsConstraints[res] + 1;
    //int *Coef = (int*)malloc(sizeof(int)*(d_cut->numberVariables));
    int *u = (int*)malloc(sizeof(int)* sz_u );
    int n1=-1, d1=-1,el, rhs, aux,value_tes ;
    int nBest=-1, dBest=-1, violation_best=0;
    aux = 0;
    for(j = d_cut->ElementsConstraints[ res ] ; j < d_cut->ElementsConstraints[ res +1 ]; j++)
    {
        u[aux] = d_cut->Coefficients[j];
        aux++;
    }
    u[sz_u - 1] = d_cut->rightSide[res];
//    __syncthreads();
    value_tes = sz_u;
    for( i = 0; i < sz_u; i++ )
    {
        for( j = 0; j < sz_u; j++ ) // sempre 1 elemento à frente
        {
            // se o (x > (x+1)) então o x passa pra frente (ordem crescente)
//            if((u[i]%u[j]==0)&&(u[i]!=0)&&(u[j]!=0)&&(i!=j))
//            {
//                u[j]=0;
//                value_tes--;
//
//            }
//            if((u[j]%u[i]==0)&&(u[i]!=0)&&(u[j]!=0)&&(i!=j))
//            {
//                u[i]=0;
//                value_tes--;
//
//            }
            if((u[i]!=0)&&(u[j]!=0)&&(i!=j))
            {
                if(u[i]%u[j]==0)
                {
                    u[j]=0;
                    value_tes--;
                }
                if(u[j]%u[i]==0)
                {
                    u[i]=0;
                    value_tes--;
                }
            }
            if (( u[i] < u[j] )&&(j>i))
            {
                aux = u[i];
                u[i] = u[j];
                u[j] = aux;
            }
        }
    }
    sz_u = value_tes;
    //for(j = d_cut->ElementsConstraints[ res ] ; j < d_cut->ElementsConstraints[ res +1 ]; j++)
    for(j = 0 ; j < sz_u; j++)
    {
        //d1 = d_cut->Coefficients[j];
        d1 = u[j];
        n1 = 1;
        while(n1<d1)
        {
            rhs = 0;
            violation = 0;
            value_tes = 0;
            rhs = d_cut->rightSide[ res ]* n1;
            if(rhs%d1!=0)
            {
                for(i = d_cut->ElementsConstraints[ res ]; i<d_cut->ElementsConstraints[ res + 1 ]; i++)
                {
                    el = d_cut->Elements[i];
                    aux = d_cut->Coefficients[i] * n1;
                    if( ((aux>0&&d1<0)||(aux<0&&d1>0))&&(aux%d1!=0))
                    {
                        aux = (aux/d1) -1;
                    }
                    else
                    {
                        aux = aux/d1;
                    }
                    //aux = aux< 0 ? (aux/d1) - 1 : aux/d1;
                    value_tes += aux*d_cut->xAsterisc[el];
                }

                if( ((rhs>0&&d1<0)||(rhs<0&&d1>0))&&(rhs%d1!=0))
                {
                    rhs = (rhs/d1) -1;
                }
                else
                {
                    rhs = rhs/d1;
                }

                if(value_tes>rhs*precision)
                {
                    violation = value_tes - (rhs*precision);
                    if(violation>violation_best)
                    {
                        violation_best = violation;
                        nBest=n1;
                        dBest=d1;
                    }
                }
            }
            n1++;
        }
    }
    free(u);
    if(violation_best!=0)
    {
        d_solution->SConst[term] = res;
        d_solution->SMult[term] = nBest;
        d_solution->SPAux[term] = dBest;
    }
    else
    {
        d_solution->SConst[term] = -1;
        d_solution->SMult[term] = -1;
        d_solution->SPAux[term] = -1;
    }

    //free(Coef);
    if(threadIdx.x == 0)
    {
        free(constraints);
    }
}


__global__ void runGPUR1_aleatory(Cut_gpu *d_cut, solutionGpu *d_solution, unsigned int *seed, curandState_t* states, int nThreads, int precision,int maxDenominator)
{

    int term = threadIdx.x + blockIdx.x*nThreads;
    __shared__ int *constraints;
    __shared__ int pos;
    curand_init(seed[term],term,0,&states[term]);

    int violation = 0, cont = 0,i;
    if(threadIdx.x == 0)
    {
        pos = 0;
        constraints = (int*)malloc(sizeof(int)*d_cut->numberConstrains);
        for(i=0; i<d_cut->numberConstrains; i++)
        {
            if(d_cut->typeConstraints[i] == RES_RR)   //&& (d_cut->typeConstraints[i] == LPC_ENERGETIC))

            {
                constraints[pos] = i;
                pos++;
            }
        }

    }
    __syncthreads();

    int res = constraints[threadIdx.x%pos];

    cont = 0;
    int n1=-1, d1=-1,el, rhs, aux,value_tes;
    int nBest=-1, dBest=-1, violation_best=0;
    while((cont<20)&&(violation_best==0))
    {
        cont++;
        d1 = curand(&states[term])%maxDenominator + 2;
        n1 = 1;
        while(n1<d1)
        {
            rhs = 0;
            violation = 0;
            value_tes = 0;
            //printf("%d/%d\n",n1,d1);
            for(i = d_cut->ElementsConstraints[ res ]; i<d_cut->ElementsConstraints[ res + 1 ]; i++)
            {
                el = d_cut->Elements[i];
                aux = d_cut->Coefficients[i] * n1;
                if( ((aux>0&&d1<0)||(aux<0&&d1>0))&&(aux%d1!=0))
                {
                    aux = (aux/d1) -1;
                }
                else
                {
                    aux = aux/d1;
                }
                value_tes += aux*d_cut->xAsterisc[el];
            }
            rhs = d_cut->rightSide[ res ]* n1;
            if( ((rhs>0&&d1<0)||(rhs<0&&d1>0))&&(rhs%d1!=0))
            {
                rhs = (rhs/d1) -1;
            }
            else
            {
                rhs = rhs/d1;
            }

            if(value_tes>rhs*precision)
            {
                violation = value_tes - (rhs*precision);
                if(violation>violation_best)
                {
                    violation_best = violation;
                    nBest=n1;
                    dBest=d1;
                }
            }
            n1++;
        }
    }

    if(violation_best!=0)
    {
        d_solution->SConst[term] = res;
        d_solution->SMult[term] = nBest;
        d_solution->SPAux[term] = dBest;
    }
    else
    {
        d_solution->SConst[term] = -1;
        d_solution->SMult[term] = -1;
        d_solution->SPAux[term] = -1;
    }

    if(threadIdx.x == 0)
    {
        free(constraints);
    }
}



__global__ void runGPUR2(Cut_gpu *d_cut, solutionGpu *d_solution, unsigned int *seed, curandState_t* states, int numberMaxConst, int *setConstraint,int nThreads, int precision, int maxDenominator)
{
    int term = threadIdx.x + blockIdx.x*nThreads;
    int mult_1, mult_2, rest_a,rest_b, i, j, el, rhs1, rhs2, value_tes, violation = 0, aux, n1_best = -1, n2_best = -1, d1_best = -1, qnt_1 = -1, d2_best=-1;//, cont=0;
    curand_init(seed[term],term,0,&states[term]);
    int Numerator[20];
    int Denominator[20];
    short int *Coef = (short int*)malloc(sizeof(short int)*(d_cut->numberVariables)); //free
    short int *Coef2 = (short int*)malloc(sizeof(short int)*(d_cut->numberVariables)); //free
    for(i=0; i<20; i++)
    {
        Denominator[i]= curand(&states[term])%maxDenominator + 2;
        Numerator[i] = curand(&states[term])%(Denominator[i]-1);

    }
    for(mult_1 = 0; mult_1 < 20; mult_1++)
    {
        memset(Coef,0,sizeof(short int)*d_cut->numberVariables);
        rhs1 = 0;
        for(rest_a = 0; rest_a< numberMaxConst; rest_a++)
        {
            for(i=d_cut->ElementsConstraints[ setConstraint[term*numberMaxConst + rest_a] ]; i<d_cut->ElementsConstraints[ setConstraint[term*numberMaxConst + rest_a] + 1]; i++)
            {

                el = d_cut->Elements[i];
                Coef[el] += d_cut->Coefficients[i] * Numerator[mult_1];
            }
            rhs1 += d_cut->rightSide[ setConstraint[term*numberMaxConst + rest_a] ] * Numerator[mult_1];
            for(mult_2 = 0; mult_2<20; mult_2++)
            {
                memset(Coef2,0,sizeof(short int)*d_cut->numberVariables);
                value_tes = 0;
                rhs2 = 0;
                for(rest_b = rest_a + 1; rest_b < numberMaxConst; rest_b++)
                {
                    for(j=d_cut->ElementsConstraints[ setConstraint[term*numberMaxConst + rest_b] ]; j<d_cut->ElementsConstraints[ setConstraint[term*numberMaxConst + rest_b] + 1]; j++)
                    {
                        el = d_cut->Elements[j];
                        Coef2[el] += d_cut->Coefficients[j] * Numerator[mult_2];
                    }
                    rhs2 += d_cut->rightSide[ setConstraint[term*numberMaxConst + rest_b] ]* Numerator[mult_2];
                }



                for(j=0; j<d_cut->numberVariables; j++)
                {
                    //if((double)Coef[j]/(double)Denominator[mult_1] + (double)Coef2[j]/(double)Denominator[mult_2] < 0 )
                    if( (Coef[j]*Denominator[mult_2] + Coef2[j]*Denominator[mult_1])/(Denominator[mult_1]*Denominator[mult_2]) < 0 )
                    {
                        aux = (Coef[j]*Denominator[mult_2] + Coef2[j]*Denominator[mult_1])/(Denominator[mult_1]*Denominator[mult_2]) - 1;

                    }
                    else
                    {
                        aux = (Coef[j]*Denominator[mult_2] + Coef2[j]*Denominator[mult_1])/(Denominator[mult_1]*Denominator[mult_2]);
                    }
                    value_tes += aux*d_cut->xAsterisc[j];
                }
                if( (rhs1*Denominator[mult_2] + rhs2*Denominator[mult_1])/(Denominator[mult_1]*Denominator[mult_2]) < 0)
                {
                    aux = (rhs1*Denominator[mult_2] + rhs2*Denominator[mult_1])/(Denominator[mult_1]*Denominator[mult_2]) - 1;
                }
                else
                {
                    aux = (rhs1*Denominator[mult_2] + rhs2*Denominator[mult_1])/(Denominator[mult_1]*Denominator[mult_2]);
                }


                if((value_tes>aux*precision)&&(value_tes-(aux*precision)>violation))
                {
                    violation = value_tes-(aux*precision);
                    n1_best = Numerator[mult_1];
                    d1_best = Denominator[mult_1];
                    n2_best = Numerator[mult_2];
                    d2_best = Denominator[mult_2];
                    qnt_1 = rest_a;
                }


            }
        }

    }
    __syncthreads();

    if(violation>0)
    {
        for(i=0; i<numberMaxConst; i++)
        {
            d_solution->SConst[i + threadIdx.x*numberMaxConst + blockIdx.x*numberMaxConst*nThreads] = setConstraint[term*numberMaxConst + i];//CPU ja vai ter
        }

        d_solution->SPAux[threadIdx.x + blockIdx.x*nThreads] = qnt_1;
        d_solution->SMult[0 + threadIdx.x*4 + blockIdx.x*4*nThreads] = n1_best;
        d_solution->SMult[1 + threadIdx.x*4 + blockIdx.x*4*nThreads] = d1_best;
        d_solution->SMult[2 + threadIdx.x*4 + blockIdx.x*4*nThreads] = n2_best;
        d_solution->SMult[3 + threadIdx.x*4 + blockIdx.x*4*nThreads] = d2_best;

    }
    else
    {
        for(i=0; i<numberMaxConst; i++)
        {
            d_solution->SConst[i + threadIdx.x*numberMaxConst + blockIdx.x*numberMaxConst*nThreads] = -1;
        }
        d_solution->SPAux[threadIdx.x + blockIdx.x*nThreads] = 0;
        d_solution->SMult[0 + threadIdx.x*4 + blockIdx.x*4*nThreads] = -1;
        d_solution->SMult[1 + threadIdx.x*4 + blockIdx.x*4*nThreads] = -1;
        d_solution->SMult[2 + threadIdx.x*4 + blockIdx.x*4*nThreads] = -1;
        d_solution->SMult[3 + threadIdx.x*4 + blockIdx.x*4*nThreads] = -1;
    }



    free(Coef);
    free(Coef2);
    __syncthreads();

}

__device__ void DecimalToBinary(int number,unsigned char v[], int sizegroup)
{
    int i=0;
    for(i = 0; i<sizegroup; i++)
    {
        v[i] = 0;
    }
    i = sizegroup - 1;
    int n = number;
    while (n > 0)
    {
        // storing remainder in binary array
        v[i] = n % 2;
        n = n / 2;
        i--;
    }


    //int j;
    //for ( j = 0; j < sizegroup; j++)
    //    printf("%d", v[j]);
    //printf("\n");
    //return v;
}

__global__ void runGPUZerohHalf(Cut_gpu *d_cut, int *d_solution, int nThread,int sizeGroup, int nBlock, int precision)
{
    unsigned char binaryTest[20];//Criar inicialmente com 20 para testes.
    int nInitial, nFinal, ite, i, j;
    __shared__ int nRunsPerThreads;
    int idUnique = threadIdx.x + blockIdx.x*nThread;
    short int *Coef = (short int*)malloc(sizeof(short int)*(d_cut->numberVariables));

    if(threadIdx.x == 0)
    {
        nRunsPerThreads = (powf(2,sizeGroup))/(nBlock*nThread);
        nRunsPerThreads++;
        //printf("t: %d\n",d_cut->numberConstrains);
        //printf("teste: %d \n",binaryTest[1]);
//            printf("Entrou no block: %d %d\n",blockIdx.x, d_solution[threadIdx.x]);
//            printf("nRunsPerThreads: %d\n",nRunsPerThreads);
    }
    __syncthreads();
    nInitial = idUnique*nRunsPerThreads;
    nFinal = ((idUnique+1)*nRunsPerThreads) - 1;
    int best_violation = 0;
    int best_number = -1;
    int rhs,el, aux;

    for(ite = nInitial; ite <= nFinal; ite++)
    {
        DecimalToBinary(ite,binaryTest,sizeGroup);
        memset(Coef,0,sizeof(short int)*d_cut->numberVariables);
        rhs = 0;
        for( i = 0 ; i<sizeGroup; i++)
        {
            if(binaryTest[i]==1)
            {
                for(j = d_cut->ElementsConstraints[ i ]; j < d_cut->ElementsConstraints[ i + 1 ]; j++)
                {
                    el = d_cut->Elements[j];
                    Coef[el] += d_cut->Coefficients[j];
                }
                rhs += d_cut->rightSide[i];
                //printf("ENTRA!!");
            }
        }
        aux = 0;
        for(j = 0; j<d_cut->numberVariables; j++)
        {
            if(Coef[j]<0)
            {
                Coef[j] = (Coef[j]/2) - 1;
            }
            else
            {
                Coef[j] = (Coef[j]/2);
            }
            aux += Coef[j]*d_cut->xAsterisc[j];
        }
        if(rhs<0)
        {
            rhs = (rhs/2)-1;
        }
        else
        {
            rhs = (rhs/2);
        }
        if((aux>rhs*precision)&&(aux-(rhs*precision)>best_violation))
        {
            best_violation = aux-(rhs*precision);
            best_number = ite;
        }
    }
    d_solution[idUnique] = best_number;
    free(Coef);
    //printf("idUnique: %d initial: %d final: %d\n", idUnique, nInitial, nFinal);
}
