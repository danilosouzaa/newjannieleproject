/*
 * solutionGpu.c
 *
 *  Created on: 29/03/2017
 *      Author: danilo
*/

#include "solutionGpu.h"

//TSMult *SMult;
//TSConst *SConst;
//TSPAux *SPAux;

solutionZeroHalf* allocationStructSolutionZeroHalf(int sizeGroup, int nRuns)
{
    size_t size_solution_zeroHalf =  sizeof(solutionZeroHalf) +
                                     sizeof(TSVetSol)*(sizeGroup*nRuns);
    solutionZeroHalf* sol_zero;
    sol_zero = (solutionZeroHalf*)malloc(size_solution_zeroHalf);
    assert(sol_zero!=NULL);
    memset(sol_zero,0,size_solution_zeroHalf);
    sol_zero->VetSol = (TSVetSol*)(sol_zero+1);
    return sol_zero;
}


solutionGpu* allocationStructSolution2(Cut_gpu *c, int numberMaxConst, int nRuns)
{
    size_t size_solution =  sizeof(solutionGpu) +
                            sizeof(TSMult)*(nRuns*4) +
                            sizeof(TSConst)*(numberMaxConst*nRuns) +
                            sizeof(TSPAux)*(nRuns);
    solutionGpu *sol;
    sol = (solutionGpu*)malloc(size_solution);
    assert(sol!=NULL);
    memset(sol,0,size_solution);
    sol->SMult = (TSMult*)(sol+1);
    sol->SConst= (TSConst*)(sol->SMult + (nRuns*4));
    sol->SPAux = (TSPAux*)(sol->SConst + (numberMaxConst*nRuns));
    return sol;
}

solutionGpu* allocationStructSolution1(Cut_gpu *c, int nRuns)  //for phase 1
{
    size_t size_solution =  sizeof(solutionGpu) +
                            sizeof(TSMult)*(nRuns) +
                            sizeof(TSConst)*(nRuns) +
                            sizeof(TSPAux)*(nRuns);
    solutionGpu *sol;
    sol = (solutionGpu*)malloc(size_solution);
    assert(sol!=NULL);
    memset(sol,0,size_solution);
    sol->SMult = (TSMult*)(sol+1);
    sol->SConst= (TSConst*)(sol->SMult + (nRuns));
    sol->SPAux = (TSPAux*)(sol->SConst + (nRuns));
    return sol;
}

int returnK(long double fa_zero)
{
    long double value;
    int arround;
    value = 1/fa_zero;
    arround = value;
    if(fabsl(value -(long double)arround) <= 1e-6)
    {
        arround--;
    }
    return arround;

}
/*
Cut_gpu* createCutsStrongPhaseOne(Cut_gpu *h_cut, solutionGpu *h_solution, int nCuts, int precision, int nRuns)
{
    int *Coef1 = (int*)malloc( sizeof(int)*(h_cut->numberVariables) );
    int *Coef_ref = (int*)malloc( sizeof(int)*(h_cut->numberVariables) );
    long double *Coef1_mult = (long double*)malloc( sizeof(long double)*(h_cut->numberVariables) );
    long double *p_frac = (long double*)malloc( sizeof(long double)*(h_cut->numberVariables) );

    double *violation = (double*)malloc(sizeof(double)*nCuts);
    double *violation_mult = (double*)malloc(sizeof(double)*nCuts);
    int i, k,j,ite, d1, n1,el, constraint, rhs, c_aux = 0, aux=0;
    long double rhs_mult, lhs,lhs_mult;

    for(i = 0; i<nRuns; i++)
    {
        if(h_solution->SConst[i]!=-1)
        {

            memset(Coef1,0,sizeof(int)*h_cut->numberVariables);
            memset(Coef_ref,0,sizeof(int)*h_cut->numberVariables);
            for(j=0; j<h_cut->numberVariables; j++)
            {
                Coef1_mult[j]=0;
                p_frac[j]=0;
            }
            //memset(Coef1_mult,0,sizeof(long double)*h_cut->numberVariables);
            lhs = 0;
            lhs_mult = 0;
            c_aux = 0;
            aux = 0;
            violation[c_aux] = 0;
            violation_mult[c_aux] = 0;
            constraint = h_solution->SConst[i];
            n1 = h_solution->SMult[i];
            d1 = h_solution->SPAux[i];

            rhs = h_cut->rightSide[constraint] * n1;
            rhs_mult =  ((long double)h_cut->rightSide[constraint] * (long double)n1)/(long double)d1;
            if((rhs<0&&d1>0) || (rhs>0&&d1<0))
            {
                rhs = (rhs/d1)-1;
            }
            else
            {
                rhs = rhs/d1;
            }
            rhs_mult = rhs_mult - rhs;

//            printf("rhs %d multiplier %d/%d\n",h_cut->rightSide[constraint], n1,d1);
            for(j = h_cut->ElementsConstraints[constraint]; j< h_cut->ElementsConstraints[constraint + 1]; j++)
            {
                el = h_cut->Elements[j];
                Coef1[el] = h_cut->Coefficients[j] * n1;

                Coef1_mult[el] = ((long double)h_cut->Coefficients[j] * (long double) n1)/(long double)d1;//fazer
//                printf("el %d coef %d\n",el,h_cut->Coefficients[j]);
                if((Coef1[el]<0&&d1>0) || (Coef1[el]>0&&d1<0))
                {

                    Coef1[el] = (Coef1[el]/d1)-1;
                }
                else
                {
                    Coef1[el] = (Coef1[el]/d1);
                }
                p_frac[el] = Coef1_mult[el] - Coef1[el];
                Coef1_mult[el] = Coef1_mult[el] - Coef1[el];// fracionaria
                if( ( (Coef1_mult[el] - rhs_mult)/(1.0-rhs_mult) ) > 0 )
                {
                    Coef1_mult[el] = Coef1[el] + (Coef1_mult[el] - rhs_mult)/(1.0-rhs_mult);
                }
                else
                {
                    Coef1_mult [el] = Coef1[el];
                }
                lhs += (long double)Coef1[el]* (long double)h_cut->xAsterisc[el];
                lhs_mult+= (long double)Coef1_mult[el]*(long double)h_cut->xAsterisc[el];
                //printf("Parte fracionaria %lf \n",Coef1_mult[el]);
                //if(Coef1_mult[el]>0)
                //    printf("%LF %LF \t",Coef1_mult[el], p_frac[el]);



            }
//            printf("\n");
            k = returnK(rhs_mult);
//            printf("%d %LF\n", k, rhs_mult);

            assert( rhs_mult<=(1.0/((long double)k)) - 1e-8 );
            assert( rhs_mult>=(1.0/((long double)k+1)) + 1e-8 );



            violation[c_aux] = lhs - (rhs*precision);
            violation_mult[c_aux] = lhs_mult - (rhs*precision);
            printf("Violation initial: %lf\n",violation[c_aux]);
            printf("Violation Strong: %lf\n", violation_mult[c_aux]);

            int *sizeRows = (int*)malloc((k+1)*sizeof(int));
            memset(sizeRows,0, (k+1)*sizeof(int));

            int **mat = (int**)malloc( (k+1)*sizeof(int*));
            for(ite = 0; ite < (k+1); ite++)
                mat[ite] = (int*)malloc(h_cut->numberVariables*sizeof(int));
            for(ite=0; ite<h_cut->numberVariables; ite++)
            {
                for(j = 0; j<=k; j++)
                {
                    mat[ite][j] = 0;
                }
            }

            for(ite=0; ite<h_cut->numberVariables; ite++)
            {
                if((p_frac[ite]<=rhs_mult)&&(Coef1_mult[ite]!=0))
                {
                    mat[0][aux] = ite;
                    aux++;
                }
            }
            sizeRows[0] = aux;
            aux = 0;
            for(j=1; j<=k; j++)
            {
                for(ite=0; ite<h_cut->numberVariables; ite++)
                {   if(Coef1_mult[ite]==0){
                        continue;
                    }
                    if((p_frac[ite] - 1e-6 >=(rhs_mult + (j-1.0)*(1.0-rhs_mult)/(long double)k))&&(p_frac[ite] + 1e-6 <=(rhs_mult + j*(1.0-rhs_mult)/(long double) k)))
                    {
                        mat[j][aux] = ite;
                        aux++;
                    }
                }
                sizeRows[j] = aux;
                aux = 0;
            }
//            printf("antes\n");
            for(ite=0; ite<sizeRows[0]; ite++)
            {
                el = mat[0][ite];
                if(Coef1_mult[el]<0)
                {

                    aux = Coef1_mult[el] - 1;
                }
                else
                {
                    aux = Coef1_mult[el];
                }
//                printf("Coef %LF aux %d  k %d  j %d el %d\n",Coef1_mult[el], aux,k,0,el);
                Coef_ref[el] = aux*(k+1);

            }
//            printf("depois\n");
            for(j=1; j<=k; j++)
            {
                for(ite = 0; ite<sizeRows[j]; ite++)
                {
                    el = mat[j][ite];

                    if(Coef1_mult[el]<0)
                    {
                        aux = Coef1_mult[el] - 1;
                    }
                    else
                    {
                        aux = Coef1_mult[el];
                    }
//                    printf("Coef %LF aux %d  k %d  j %d el %d\n",Coef1_mult[el], aux,k,j,el);
                    Coef_ref[el] = aux*(k+1) + j;
                }

            }
            rhs = (k+1)*rhs;
            lhs = 0;

            for(ite = 0; ite<h_cut->numberVariables; ite++)
            {
                lhs += (long double)Coef_ref[ite]*(long double)h_cut->xAsterisc[ite];
            }
            violation_mult[c_aux] = lhs - (rhs*precision);
            printf("Violation 2: %lf \n",violation_mult[c_aux]);
            for (ite = 0; ite < k+1; ite++)
            {
                free (mat[ite]);
            }
            free (mat);
            free(sizeRows);
            c_aux++;
        }

    }

    free(p_frac);
    free(Coef1);
    free(Coef1_mult);
    free(Coef_ref);
    free(violation);
    free(violation_mult);
    return h_cut;
}*/

long double closestIntegerDistance(long double v)
{
    long double distLow = v - floorl(v);
    long double distUP = ceill(v) - v;
    if (distLow<distUP)
        return distLow;
    return distUP;
}

long double ldRound( long double v)
{
    return v + 0.5;
}
/*
Cut_gpu* createCutsStrongPhaseTwo(Cut_gpu *h_cut, solutionGpu *h_solution, int numberMaxConst, int nCuts, int precision, int nRuns, int nThreads,int nBlocks)
{
    double value = 0 ;
    double *value_violation = (double*)malloc(sizeof(double)*nCuts);
    double *value_violation_r1 = (double*)malloc(sizeof(double)*nCuts);
    double *violation_mult = (double*)malloc(sizeof(double)*nCuts);
    int *Coef1 = (int*)malloc(sizeof(int)*(h_cut->numberVariables));
    int *Coef2 = (int*)malloc(sizeof(int)*(h_cut->numberVariables));
    int *Coef_ref = (int*)malloc(sizeof(int)*(h_cut->numberVariables));
    long double *Coef_mult1 = (long double*)malloc(sizeof(long double)*h_cut->numberVariables);
    long double *p_frac = (long double*)malloc(sizeof(long double)*h_cut->numberVariables);


    int *Coefs_temp = (int*)malloc(sizeof(int)*(nCuts*h_cut->numberVariables));
    int *Elements_temp = (int*)malloc(sizeof(int)*(nCuts*h_cut->numberVariables));
    int *Pos_el_temp = (int*)malloc(sizeof(int)*(nCuts+1) );
    int *rhs_temp = (int*)malloc(sizeof(int)*(nCuts));
    long double rhs_mult1, rhs_mult2;
    int nAll = 0;
    int temp, k_R2;
    int i,k,j,n_1,d_1,n_2,d_2, qnt_1, el, rest_a, rhs1 = 0, rhs2 = 0,aux = 0, cont_aux=0,ite;
    Pos_el_temp[0] = 0;
    for (i=0; i<nThreads; i++)
    {
        for(k=0; k<nBlocks; k++)
        {
            if(h_solution->SConst[0 + i*numberMaxConst + k*numberMaxConst*nThreads]!=-1)
            {
                for(ite=0; ite<h_cut->numberVariables; ite++)
                {
                    Coef_mult1[ite] = 0;
                    Coefs_temp[ite] = 0;
                    p_frac[ite] = 0;
                    Elements_temp[ite] = 0;

                }
                memset(Coef1,0,sizeof(int)*h_cut->numberVariables);
                memset(Coef2,0,sizeof(int)*h_cut->numberVariables);
                memset(Coef_ref,0,sizeof(int)*h_cut->numberVariables);
                rhs1 = 0;
                rhs2 = 0;
                rhs_mult1 = 0;
                rhs_mult2 = 0;
                qnt_1 = h_solution->SPAux[i + k*nThreads];
                n_1 = h_solution->SMult[0 + i*4 + k*4*nThreads];
                d_1 = h_solution->SMult[1 + i*4 + k*4*nThreads];
                n_2 = h_solution->SMult[2 + i*4 + k*4*nThreads];
                d_2 = h_solution->SMult[3 + i*4 + k*4*nThreads];

                for(rest_a = 0; rest_a <= qnt_1; rest_a++)
                {
                    for(j = h_cut->ElementsConstraints[ h_solution->SConst[rest_a + i*numberMaxConst + k*numberMaxConst*nThreads] ]; j< h_cut->ElementsConstraints[ h_solution->SConst[rest_a + i*numberMaxConst + k*numberMaxConst*nThreads] + 1]; j++)
                    {
                        el = h_cut->Elements[j];
                        Coef1[el] += h_cut->Coefficients[j] * n_1;
                        //value_tes += Coef[el] * d_cut->xAsterisc[el];
                    }
                    rhs1 += h_cut->rightSide[ h_solution->SConst[rest_a+ i*numberMaxConst + k*numberMaxConst*nThreads] ] * n_1;

                }




                for(rest_a = qnt_1 + 1; rest_a < numberMaxConst; rest_a++)
                {
                    for(j=h_cut->ElementsConstraints[ h_solution->SConst[rest_a+ i*numberMaxConst + k*numberMaxConst*nThreads] ]; j<h_cut->ElementsConstraints[  h_solution->SConst[rest_a + i*numberMaxConst + k*numberMaxConst*nThreads] + 1]; j++)
                    {
                        el = h_cut->Elements[j];
                        Coef2[el] += h_cut->Coefficients[j] * n_2;
                        //value_tes += Coef[el] * d_cut->xAsterisc[el];
                    }
                    rhs2 += h_cut->rightSide[ h_solution->SConst[rest_a+ i*numberMaxConst + k*numberMaxConst*nThreads] ] * n_2;
                }

                value_violation[cont_aux] = 0;
                value_violation_r1[cont_aux] = 0;
                violation_mult[cont_aux] = 0;
                for(j=0; j<h_cut->numberVariables; j++)
                {
                    //temp = (long double)Coef1[j]/(long double)d_1 + (long double)Coef2[j]/(long double)d_2;
                    Coef_mult1[j] = (long double)Coef1[j]/(long double)d_1 + (long double)Coef2[j]/(long double)d_2; //Coef Double
//                    if(closestIntegerDistance(Coef_mult1[j]) <= 1e-6 ){
//                        temp = ldRound(Coef_mult1[j]);
//                    }else{
//                        temṕ = floorl(Coef_mult1[j]);
//                    }
                    if( (Coef1[j]*d_2 + Coef2[j]*d_1)/(d_1*d_2)<0){
                        temp = ((Coef1[j]*d_2 + Coef2[j]*d_1)/(d_1*d_2)) - 1;
                    }else{
                        temp = (Coef1[j]*d_2 + Coef2[j]*d_1)/(d_1*d_2);
                    }

                    p_frac[j] = Coef_mult1[j] - temp; //Coef Frac

                    if(temp != 0)
                    {
                        Coefs_temp[aux] = temp;
                        Coef1[j] = Coefs_temp[aux];
                        Elements_temp[aux] = j;
                        value_violation[cont_aux] += (Coefs_temp[aux])*h_cut->xAsterisc[j];
                        aux++;
                    }
                    else
                    {
                        Coef1[j] = 0;
                    }

                }
                Pos_el_temp[cont_aux + 1] = aux;
                //temp = ((double)rhs1/(double)d_1 + (double)rhs2/(double)d_2);

                rhs_mult1 = ((double)rhs1/(double)d_1 + (double)rhs2/(double)d_2);// RHS Double
                if((rhs1*d_2 + rhs2*d_1)/(d_1*d_2)<0){
                    temp = ((rhs1*d_2 + rhs2*d_1)/(d_1*d_2))-1;
                }else{
                    temp = (rhs1*d_2 + rhs2*d_1)/(d_1*d_2);
                }

//                if(closestIntegerDistance(rhs_mult1) <= 1e-6 ){
//                        temp = ldRound(rhs_mult1);
//                }else{
//                        temp = floorl(rhs_mult1);
//                }

                rhs_mult2 = rhs_mult1 - temp; //Frac RHS

                for(ite=0; ite<h_cut->numberVariables; ite++)
                {
                    if((p_frac[ite] - rhs_mult2)/(1.0 - rhs_mult2)>0)
                    {
                        Coef_mult1[ite] = Coef1[ite] + (p_frac[ite] - rhs_mult2)/(1.0 - rhs_mult2);
                    }
                    else
                    {
                        Coef_mult1[ite] = Coef1[ite];
                    }
                    value_violation_r1[cont_aux] += (Coef_mult1[ite])*h_cut->xAsterisc[ite];
                }

                value_violation[cont_aux] -= (rhs_temp[cont_aux])*precision;
                value_violation_r1[cont_aux] -= rhs_temp[cont_aux]*precision;
                printf("Violation Initial: %lf\n",value_violation[cont_aux]);
                printf("Violation reforço 1: %lf\n",value_violation_r1[cont_aux]);


                if(rhs_mult2>0)
                {
                    k_R2 = returnK(rhs_mult2);



                    //if(value_violation[cont_aux]>=precision){
                    //    printf("Value violation: %f %d\n", value_violation[cont_aux],precision);
                    //    for(rest_a = 0; rest_a< numberMaxConst; rest_a++){
                    //        int testando = h_solution->SConst[rest_a+ i*numberMaxConst + k*numberMaxConst*nThreads];
                    //        printf("%d \t", testando);
                    //        printf("Tipo: %d \n",h_cut->typeConstraints[testando]);
                    //    }
                    //    printf("\n");
                    //:}

                    int *sizeRows = (int*)malloc((k_R2+1)*sizeof(int));
                    memset(sizeRows,0, (k_R2+1)*sizeof(int));
                    int **mat = (int**)malloc( (k_R2+1)*sizeof(int*));
                    for(ite = 0; ite < (k_R2+1); ite++)
                        mat[ite] = (int*)malloc(h_cut->numberVariables*sizeof(int));
                    for(ite=0; ite< (h_cut->numberVariables); ite++)
                    {
                        for(j=0; j<=k_R2; j++)
                        {
                            mat[ite][j] = 0;
                        }
                    }
                    temp = 0;
                    for(ite=0; ite<h_cut->numberVariables; ite++)
                    {
                        if((p_frac[ite]<=rhs_mult2)&&(Coef_mult1[ite]!=0))
                        {
                            mat[0][temp] = ite;
                            temp++;
                        }
                    }
                    sizeRows[0] = temp;
                    temp = 0;
                    for(j=1; j<=k_R2; j++)
                    {
                        for(ite=0; ite<h_cut->numberVariables; ite++)
                        {
                            if((p_frac[ite]>(rhs_mult2 + (j-1)*(1-rhs_mult2)/k_R2))&&(p_frac[ite]<=(rhs_mult2 + j*(1-rhs_mult2)/k_R2))&&(Coef_mult1[ite]!=0))
                            {
                                mat[j][temp] = ite;
                                temp++;
                            }
                        }
                        sizeRows[j] = temp;
                        temp = 0;
                    }
//            printf("antes\n");
                    for(ite=0; ite<sizeRows[0]; ite++)
                    {
                        el = mat[0][ite];
                        if(Coef_mult1[el]<0)
                        {

                            temp = Coef_mult1[el] - 1;
                        }
                        else
                        {
                            temp = Coef_mult1[el];
                        }
//                printf("Coef %LF aux %d  k %d  j %d el %d\n",Coef1_mult[el], aux,k,0,el);
                        Coef_ref[el] = temp*(k_R2+1);

                    }
//            printf("depois\n");
                    for(j=1; j<=k_R2; j++)
                    {
                        for(ite = 0; ite<sizeRows[j]; ite++)
                        {
                            el = mat[j][ite];

                            if(Coef_mult1[el]<0)
                            {
                                temp = Coef_mult1[el] - 1;
                            }
                            else
                            {
                                temp = Coef_mult1[el];
                            }
//                    printf("Coef %LF aux %d  k %d  j %d el %d\n",Coef1_mult[el], aux,k,j,el);
                            Coef_ref[el] = temp*(k_R2+1) + j;
                        }

                    }
                    rhs_temp[cont_aux] = (k_R2+1)*rhs_temp[cont_aux];
                    int lhs = 0;

                    for(ite = 0; ite<h_cut->numberVariables; ite++)
                    {
                        lhs += (long double)Coef_ref[ite]*(long double)h_cut->xAsterisc[ite];
                    }
                    violation_mult[cont_aux] = lhs - (rhs_temp[cont_aux]*precision);
                    printf("Violation 2: %lf \n",violation_mult[cont_aux]);
                    for (ite = 0; ite < k_R2+1; ite++)
                    {
                        free (mat[ite]);
                    }
                    free (mat);
                    free(sizeRows);
                }
                cont_aux++;
            }
        }
    }
    free(Coef_ref);
    free(value_violation_r1);
    free(value_violation);
    free(Coef1);
    free(Coef2);
    free(Coef_mult1);
    free(Coefs_temp);
    free(Elements_temp);
    free(Pos_el_temp);
    free(rhs_temp);
    free(p_frac);
    free(violation_mult);
    return h_cut;
}

*/
int verifyDominanceCG(int *v1, int rhs1, int *v2, int rhs2, int sz)
{
    int v1_temp[sz];
    int v2_temp[sz];
    int i, cont = 0;
    for(i=0; i<sz; i++)
    {
        v1_temp[i] = v1[i]*rhs2;
        v2_temp[i] = v2[i]*rhs1;
    }
    rhs2 = rhs2*rhs1;
    rhs1 = rhs2;
    //int dominance = 1;
    for(i=0; i<sz; i++)
    {
        if(v1_temp[i]<v2_temp[i])
        {
            return 0;
        }
        else if(v1_temp[i]>v2_temp[i])
        {
            cont++;
        }
    }
    if(cont>0)
        return 1;
    else
        return 0;
}



Cut_gpu* createCutsStrongPhaseOne(Cut_gpu *h_cut, solutionGpu *h_solution, int nCuts, int precision, int nRuns,int nC_initial, int cont_ini)
{
    int *Coef1 = (int*)malloc( sizeof(int)*(h_cut->numberVariables) );
    int *Coef_ref = (int*)malloc( sizeof(int)*(h_cut->numberVariables) );
    long double *Coef1_mult = (long double*)malloc( sizeof(long double)*(h_cut->numberVariables) );
    long double *p_frac = (long double*)malloc( sizeof(long double)*(h_cut->numberVariables) );
    int *v_aux = (int*)malloc(sizeof(int)*(h_cut->numberVariables + 1));

    Cut_gpu *cuts_generated;

    double *violation = (double*)malloc(sizeof(double)*nCuts);
    double *violation_mult = (double*)malloc(sizeof(double)*nCuts);

    int *Coefs_temp = (int*)malloc(sizeof(int)*(nCuts*h_cut->numberVariables));
    int *Elements_temp = (int*)malloc(sizeof(int)*(nCuts*h_cut->numberVariables));
    int *Pos_el_temp = (int*)malloc(sizeof(int)*(nCuts+1) );
    int *rhs_temp = (int*)malloc(sizeof(int)*(nCuts));
//    int vio_aux=0; //retirar
    int i, k,j,ite, d1, n1,el, constraint, rhs,rhs_ref, c_aux = 0, aux=0, tam;// t_mult;
    long double rhs_mult, lhs;//lhs_mult,f_rhs;
    int cont_aux = 0;

    Pos_el_temp[0] = 0;
//    printf("Initial phase 1!\n");
    for(i = 0; i<nRuns; i++)
    {
        if(h_solution->SConst[i]!=-1)
        {
            memset(Coef1,0,sizeof(int)*h_cut->numberVariables);
            memset(Coef_ref,0,sizeof(int)*h_cut->numberVariables);
            for(j=0; j<h_cut->numberVariables; j++)
            {
                Coef1_mult[j]=0;
                p_frac[j]=0;
            }
            //memset(Coef1_mult,0,sizeof(long double)*h_cut->numberVariables);
            lhs = 0;
//            lhs_mult = 0;
            //c_aux = 0;
            aux = 0;
            tam = 0;
            violation[c_aux] = 0;
            violation_mult[c_aux] = 0;
            constraint = h_solution->SConst[i];
            n1 = h_solution->SMult[i];
            d1 = h_solution->SPAux[i];

            rhs = h_cut->rightSide[constraint] * n1;
            rhs_mult =  ((long double)h_cut->rightSide[constraint] * (long double)n1)/(long double)d1;
            if((rhs<0&&d1>0) || (rhs>0&&d1<0))
            {
                rhs = (rhs/d1)-1;
            }
            else
            {
                rhs = rhs/d1;
            }
            rhs_mult = rhs_mult - rhs;
            //printf("value f(rhs): %LF\n",rhs_mult);
//            if(rhs_mult<0.5)
//            {
//                t_mult = 1;
//                do
//                {
//                    t_mult++;
//                }
//                while(t_mult*rhs_mult<0.5);
//                //printf("t_mult = %d\n",t_mult);
//                aux = (n1*t_mult)/d1;
//                n1 = n1*t_mult - aux*d1;
//                rhs = h_cut->rightSide[constraint] * n1;
//                rhs_mult =  ((long double)h_cut->rightSide[constraint] * (long double)n1)/(long double)d1;
//                if((rhs<0&&d1>0) || (rhs>0&&d1<0))
//                {
//                    rhs = (rhs/d1)-1;
//                }
//                else
//                {
//                    rhs = rhs/d1;
//                }
//                rhs_mult = rhs_mult - rhs;
//            }
            aux = 0;
//            printf("rhs %d multiplier %d/%d\n",h_cut->rightSide[constraint], n1,d1);
            for(j = h_cut->ElementsConstraints[constraint]; j< h_cut->ElementsConstraints[constraint + 1]; j++)
            {
                el = h_cut->Elements[j];
                Coef1[el] = h_cut->Coefficients[j] * n1;
                Coef1_mult[el] = ((long double)h_cut->Coefficients[j] * (long double) n1)/(long double)d1;//fazer
//                printf("el %d coef %d\n",el,h_cut->Coefficients[j]);
                if((Coef1[el]<0&&d1>0) || (Coef1[el]>0&&d1<0))
                {

                    Coef1[el] = (Coef1[el]/d1)-1;
                }
                else
                {
                    Coef1[el] = (Coef1[el]/d1);
                }
                p_frac[el] = Coef1_mult[el] - Coef1[el];
                //Coef1_mult[el] = Coef1_mult[el] - Coef1[el];// fracionaria
                //if( ( (Coef1_mult[el] - rhs_mult)/(1.0-rhs_mult) ) > 0 )
                //{
                //    Coef1_mult[el] = Coef1[el] + (Coef1_mult[el] - rhs_mult)/(1-rhs_mult);
                //}
                //else
                //{
                //Coef1_mult [el] = Coef1[el];
                //}
                lhs += (long double)Coef1[el]* (long double)h_cut->xAsterisc[el];
                //lhs_mult+= (long double)Coef1_mult[el]*(long double)h_cut->xAsterisc[el];
                //printf("Parte fracionaria %lf \n",Coef1_mult[el]);
                //if(Coef1_mult[el]>0)
                //    printf("%LF %LF \t",Coef1_mult[el], p_frac[el]);



            }

            violation[c_aux] = lhs - (rhs*precision);
            //violation_mult[c_aux] = lhs_mult - (rhs*precision);
//            printf("Violation initial: %lf\n",violation[c_aux]);
//            printf("Violation Strong: %lf\n", violation_mult[c_aux]);
            //vio_aux = violation_mult[c_aux];//retirar


//            printf("\n");
            if(rhs_mult - 1e-8  > 0)
            {
                k = returnK(rhs_mult);
//            printf("%d %LF\n", k, rhs_mult);
                //printf("%LF %d\n",rhs_mult, k);
                assert( rhs_mult<=(1.0/((long double)k)) + 1e-8 );
                assert( rhs_mult>=(1.0/((long double)k+1)) - 1e-8 );

//                for(ite = 0; ite<h_cut->numberVariables; ite++ )
//                {
//                    if(Coef1[ite]!=0)
//                    {
//                        printf("%d v_%d \t", Coef1[ite], ite );
//                    }
//                }
//                printf("<= %d\n", rhs);
//                printf("reforço 1\n");
//                for(ite = 0; ite<h_cut->numberVariables; ite++ )
//                {
//                    if(Coef1_mult[ite]!=0)
//                    {
//                        printf("%LF v_%d \t", Coef1_mult[ite], ite );
//                    }
//                }
//                printf("<= %d\n", rhs);

                int *sizeRows = (int*)malloc((k+1)*sizeof(int));
                int **mat = (int**)malloc( (k+1)*sizeof(int*));
                for(ite = 0; ite < (k+1); ite++)
                    mat[ite] = (int*)malloc(h_cut->numberVariables*sizeof(int));


                for(ite=0; ite<h_cut->numberVariables; ite++)
                {
                    if(p_frac[ite] - 1e-6 <=rhs_mult)
                    {
                        mat[0][aux] = ite;
                        aux++;
                    }
                }
                sizeRows[0] = aux;
                aux = 0;
                for(j=1; j<=k; j++)
                {
                    for(ite=0; ite<h_cut->numberVariables; ite++)
                    {
                        if(Coef1_mult[ite]==0.0)
                        {
                            continue;
                        }
                        if((p_frac[ite] - 1e-6 >=(rhs_mult + (j-1.0)*(1.0-rhs_mult)/(long double)k))&&(p_frac[ite] /* + 1e-6*/ <=(rhs_mult + j*(1.0-rhs_mult)/(long double)k)))
                        {
                            mat[j][aux] = ite;
                            aux++;
                        }
                    }
                    sizeRows[j] = aux;
                    aux = 0;
                }
//            printf("antes\n");
                for(ite=0; ite<sizeRows[0]; ite++)
                {
                    el = mat[0][ite];
                    if(Coef1_mult[el]<0)
                    {

                        aux = Coef1_mult[el] - 1;
                    }
                    else
                    {
                        aux = Coef1_mult[el];
                    }
//                printf("Coef %LF aux %d  k %d  j %d el %d\n",Coef1_mult[el], aux,k,0,el);
                    Coef_ref[el] = aux*(k+1);

                }
//            printf("depois\n");
                for(j=1; j<=k; j++)
                {
                    for(ite = 0; ite<sizeRows[j]; ite++)
                    {
                        el = mat[j][ite];

                        if(Coef1_mult[el]<0)
                        {
                            aux = Coef1_mult[el] - 1;
                        }
                        else
                        {
                            aux = Coef1_mult[el];
                        }
//                    printf("Coef %LF aux %d  k %d  j %d el %d\n",Coef1_mult[el], aux,k,j,el);
                        Coef_ref[el] = aux*(k+1) + j;
                    }

                }
                rhs_ref = (k+1)*rhs;
                lhs = 0;

                for(ite = 0; ite<h_cut->numberVariables; ite++)
                {
                    lhs += (long double)Coef_ref[ite]*(long double)h_cut->xAsterisc[ite];
                }
                violation_mult[c_aux] = lhs - (rhs_ref*precision);


//                printf("reforço 2\n");
//                for(ite = 0; ite<h_cut->numberVariables; ite++ )
//                {
//                    if(Coef_ref[ite]!=0)
//                    {
//                        printf("%d v_%d \t", Coef_ref[ite], ite );
//                    }
//                }
//                printf("<= %d\n", rhs_ref);
//                for(ite = 0; ite<h_cut->numberVariables; ite++ )
//                {
//                    printf("v_%d = %d\n", ite, h_cut->xAsterisc[ite]);
//                }
//
                //getchar();

//                printf("Violation 2: %lf \n",violation_mult[c_aux]);

                for (ite = 0; ite < k+1; ite++)
                {
                    free (mat[ite]);
                }
                free (mat);
                free(sizeRows);
                tam = 0;
                if(verifyDominanceCG(Coef_ref,rhs_ref,Coef1,rhs,h_cut->numberVariables)==1)
                {
                    for(j=0; j<h_cut->numberVariables; j++)
                    {
                        if(Coef_ref[j]!=0)
                        {
                            v_aux[tam] = Coef_ref[j];
                            tam++;
                        }
                    }
                    v_aux[tam] = rhs_ref;
                    //tam++;
                    int mdc = CutP_maxDivisorCommonVector(v_aux, tam);
                    for(j=0; j<h_cut->numberVariables; j++)
                    {
                        if(Coef_ref[j]!=0)
                        {
                            Coefs_temp[cont_aux] = Coef_ref[j]/mdc;
                            //printf("%d\t",Coefs_temp[cont_aux]);
                            Elements_temp[cont_aux] = j;
                            cont_aux++;
                        }
                    }
                    Pos_el_temp[c_aux+1] = cont_aux;
                    rhs_temp[c_aux] = rhs_ref/mdc;

                }
                else
                {
                    for(j=0; j<h_cut->numberVariables; j++)
                    {
                        if(Coef1[j]!=0)
                        {
                            v_aux[tam] = Coef1[j];
                            tam++;
                        }
                    }
                    v_aux[tam] = rhs;
                    //tam++;
                    int mdc = CutP_maxDivisorCommonVector(v_aux, tam);
                    for(j=0; j<h_cut->numberVariables; j++)
                    {
                        if(Coef1[j]!=0)
                        {
                            Coefs_temp[cont_aux] = Coef1[j]/mdc;
                            //printf("%d\t",Coefs_temp[cont_aux]);
                            Elements_temp[cont_aux] = j;
                            cont_aux++;
                        }
                    }
                    Pos_el_temp[c_aux+1] = cont_aux;
                    rhs_temp[c_aux] = rhs/mdc;
                }
//                if(vio_aux>violation[c_aux])
//                    exit(0);
            }
            else
            {
                tam = 0;
                for(j=0; j<h_cut->numberVariables; j++)
                {
                    if(Coef1[j]!=0)
                    {
                        v_aux[tam] = Coef1[j];
                        tam++;
                    }
                }
                v_aux[tam] = rhs;
                //tam++;
                int mdc = CutP_maxDivisorCommonVector(v_aux, tam);
                for(j=0; j<h_cut->numberVariables; j++)
                {
                    if(Coef1[j]!=0)
                    {
                        Coefs_temp[cont_aux] = Coef1[j]/mdc;
                        //printf("%d\t",Coefs_temp[cont_aux]);
                        Elements_temp[cont_aux] = j;
                        cont_aux++;
                    }
                }
                Pos_el_temp[c_aux+1] = cont_aux;
                rhs_temp[c_aux] = rhs/mdc;
            }
            c_aux++;
        }
    }

    free(Coef1);
    free(Coef1_mult);
    free(Coef_ref);
    free(violation);
    free(violation_mult);
    free(p_frac);
    free(v_aux);

    cuts_generated = AllocationStructCut(cont_ini + cont_aux, nC_initial + nCuts,h_cut->numberVariables);
    for(i=0; i<cuts_generated->numberVariables; i++)
    {
        cuts_generated->xAsterisc[i] = h_cut->xAsterisc[i];
    }
    cuts_generated->ElementsConstraints[0] = 0;
    for(i=0; i<nC_initial; i++)
    {
        cuts_generated->rightSide[i] = h_cut->rightSide[i];
        cuts_generated->typeConstraints[i] = h_cut->typeConstraints[i];
        cuts_generated->ElementsConstraints[i+1] = h_cut->ElementsConstraints[i+1];
    }
    aux = 1;

    for(i = nC_initial; i<cuts_generated->numberConstrains; i++)
    {
        cuts_generated->rightSide[i] = rhs_temp[i - nC_initial];
        cuts_generated->typeConstraints[i] = LPC_CGGPU;
        cuts_generated->ElementsConstraints[i+1] = Pos_el_temp[aux] + cont_ini;
        aux++;//+ h_cut->cont;
    }
    for(i=0; i<cont_ini; i++)
    {
        cuts_generated->Coefficients[i] = h_cut->Coefficients[i];
        cuts_generated->Elements[i] = h_cut->Elements[i];
    }
    for(i= cont_ini; i<cuts_generated->cont; i++)
    {
        cuts_generated->Coefficients[i] = Coefs_temp[i - cont_ini];
        cuts_generated->Elements[i] = Elements_temp[i - cont_ini];
    }


    int *validated = (int*)malloc(sizeof(int)*cuts_generated->numberConstrains);
    memset(validated,0,sizeof(int)*cuts_generated->numberConstrains);
    aux = 1;
    cont_aux=0;
    int p1,p2, minus_elements = 0;
    k=0;
    for(i=0; i<cuts_generated->numberConstrains - 1; i++)
    {
        for(j=i+1; j<cuts_generated->numberConstrains; j++)
        {
            if(validated[j]==0)
            {
                if((cuts_generated->ElementsConstraints[i+1]-cuts_generated->ElementsConstraints[i])!=(cuts_generated->ElementsConstraints[j+1]-cuts_generated->ElementsConstraints[j])||(cuts_generated->rightSide[i]!=cuts_generated->rightSide[j]))
                {
                    aux=0;
                }
                else
                {
                    p1 = cuts_generated->ElementsConstraints[i];
                    p2 = cuts_generated->ElementsConstraints[j];
                    aux = 1;
                    for(k=0; k<cuts_generated->ElementsConstraints[i+1]-cuts_generated->ElementsConstraints[i]; k++)
                    {

                        if((cuts_generated->Coefficients[p1+k]!=cuts_generated->Coefficients[p2+k])||(cuts_generated->Elements[p1+k]!=cuts_generated->Elements[p2+k]))
                        {
                            aux=0;

                            break;
                        }
                    }
                }
                if((aux==1)&&(j >= h_cut->numberConstrains))
                {
                    validated[j]=1;
                    minus_elements += cuts_generated->ElementsConstraints[j+1]-cuts_generated->ElementsConstraints[j];
                    cont_aux++;
                }
            }
        }
    }

//    printf("Number of repeat Strengthening: %d \n", cont_aux);

    Cut_gpu* new_h_cut;

    new_h_cut = AllocationStructCut(cuts_generated->cont - minus_elements,cuts_generated->numberConstrains-cont_aux,cuts_generated->numberVariables);
    aux = 0;
    cont_aux = 0;
    new_h_cut->ElementsConstraints[0] = 0;
    for(i=0; i<cuts_generated->numberConstrains; i++)
    {
        if(validated[i]==0)
        {
//            if(i>=h_cut->numberConstrains)
//            {
//                printf("violation %f\n",violation[i-h_cut->numberConstrains]);
//            }
            new_h_cut->rightSide[aux] = cuts_generated->rightSide[i];
            new_h_cut->typeConstraints[aux] = cuts_generated->typeConstraints[i];
            for(j = cuts_generated->ElementsConstraints[i]; j < cuts_generated->ElementsConstraints[i+1]; j++)
            {
                new_h_cut->Coefficients[cont_aux] = cuts_generated->Coefficients[j];
                new_h_cut->Elements[cont_aux] = cuts_generated->Elements[j];
                cont_aux++;
            }
            new_h_cut->ElementsConstraints[aux + 1] = cont_aux;

            aux++;
        }

    }
    for(i=0; i<new_h_cut->numberVariables; i++)
    {
        new_h_cut->xAsterisc[i] = cuts_generated->xAsterisc[i];
    }
//
//    printf("End phase 1!\n");
    free(Coefs_temp);
    free(Elements_temp);
    free(Pos_el_temp);
    free(rhs_temp);
    free(validated);
    free(cuts_generated);

    return new_h_cut;
}


void DecimalToBinaryCpu(int number,unsigned char v[], int sizegroup)
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
}

Cut_gpu* createCutsStrongZeroHalf(Cut_gpu *h_cut, int *h_solution, int sizeGroup, int nBlocks, int nThreads, int precision, int nCuts, int cont_ini, int nC_initial)
{
    int *Coef1 = (int*)malloc( sizeof(int)*(h_cut->numberVariables) );
    int *Coef_ref = (int*)malloc( sizeof(int)*(h_cut->numberVariables) );
    long double *Coef1_mult = (long double*)malloc( sizeof(long double)*(h_cut->numberVariables) );
    long double *p_frac = (long double*)malloc( sizeof(long double)*(h_cut->numberVariables) );
    int *v_aux = (int*)malloc(sizeof(int)*(h_cut->numberVariables + 1));


    unsigned char *v_sol_bin = (unsigned char*)malloc(sizeof(unsigned char)*(nBlocks*nThreads));


    Cut_gpu *cuts_generated;

    double *violation = (double*)malloc(sizeof(double)*nCuts);
    double *violation_mult = (double*)malloc(sizeof(double)*nCuts);

    int *Coefs_temp = (int*)malloc(sizeof(int)*(nCuts*h_cut->numberVariables));
    int *Elements_temp = (int*)malloc(sizeof(int)*(nCuts*h_cut->numberVariables));
    int *Pos_el_temp = (int*)malloc(sizeof(int)*(nCuts+1) );
    int *rhs_temp = (int*)malloc(sizeof(int)*(nCuts));
//    int vio_aux=0; //retirar
    int i, k,j,ite, d1, n1,el, constraint, rhs,rhs_ref, c_aux = 0, aux=0, tam;// t_mult;
    long double rhs_mult, lhs;//lhs_mult,f_rhs;
    int cont_aux = 0;

    Pos_el_temp[0] = 0;
//    printf("Initial phase 1!\n");
    int nRuns = nThreads*nBlocks;
    for(i = 0; i<nRuns; i++)
    {
        if(h_solution[i]!=-1)
        {

            memset(Coef1,0,sizeof(int)*h_cut->numberVariables);
            memset(Coef_ref,0,sizeof(int)*h_cut->numberVariables);
            for(j=0; j<sizeGroup; j++)
            {
                v_sol_bin[j] = 0;
            }
            DecimalToBinaryCpu(h_solution[i],v_sol_bin,sizeGroup);

            for(j=0; j<h_cut->numberVariables; j++)
            {
                Coef1_mult[j]=0;
                p_frac[j]=0;
            }
            lhs = 0;
            aux = 0;
            tam = 0;
            violation[c_aux] = 0;
            violation_mult[c_aux] = 0;
            rhs = 0;
            for(j=0; j<sizeGroup; j++)
            {
                rhs += h_cut->rightSide[j]*v_sol_bin[j];
            }
            rhs_mult = (long double)rhs/ (long double) 2;
            rhs = h_cut->rightSide[constraint] * n1;
            if(rhs<0)
            {
                rhs = (rhs/2) - 1;
            }
            else
            {
                rhs = rhs/2;
            }
            rhs_mult = rhs_mult - rhs;//Fractional Part


            aux = 0;
            for(k=0; k<sizeGroup; k++)
            {
                if(v_sol_bin[k]==1)
                {
                    constraint = k;
                    for(j = h_cut->ElementsConstraints[constraint]; j< h_cut->ElementsConstraints[constraint + 1]; j++)
                    {
                        el = h_cut->Elements[j];
                        Coef1[el] += h_cut->Coefficients[j];
                    }
                }
            }
            for(j=0;j<h_cut->numberVariables;j++){
                Coef1_mult[j] = (long double)Coef1[j]/(long double)2;
                if(Coef1[j]<0){
                    Coef1[j] = (Coef1[j]/2) - 1;
                }else{
                    Coef1[j] = Coef1[j]/2;
                }
                p_frac[j] = Coef1_mult[j] - Coef1[j];
                lhs += Coef1[j]*h_cut->xAsterisc[j];
            }

            violation[c_aux] = lhs - (rhs*precision);

//            printf("\n");
            if(rhs_mult - 1e-8  > 0)
            {
                k = returnK(rhs_mult);
//            printf("%d %LF\n", k, rhs_mult);
                //printf("%LF %d\n",rhs_mult, k);
                assert( rhs_mult<=(1.0/((long double)k)) + 1e-8 );
                assert( rhs_mult>=(1.0/((long double)k+1)) - 1e-8 );


                int *sizeRows = (int*)malloc((k+1)*sizeof(int));
                int **mat = (int**)malloc( (k+1)*sizeof(int*));
                for(ite = 0; ite < (k+1); ite++)
                    mat[ite] = (int*)malloc(h_cut->numberVariables*sizeof(int));


                for(ite=0; ite<h_cut->numberVariables; ite++)
                {
                    if(p_frac[ite] - 1e-6 <=rhs_mult)
                    {
                        mat[0][aux] = ite;
                        aux++;
                    }
                }
                sizeRows[0] = aux;
                aux = 0;
                for(j=1; j<=k; j++)
                {
                    for(ite=0; ite<h_cut->numberVariables; ite++)
                    {
                        if(Coef1_mult[ite]==0.0)
                        {
                            continue;
                        }
                        if((p_frac[ite] - 1e-6 >=(rhs_mult + (j-1.0)*(1.0-rhs_mult)/(long double)k))&&(p_frac[ite] /* + 1e-6*/ <=(rhs_mult + j*(1.0-rhs_mult)/(long double)k)))
                        {
                            mat[j][aux] = ite;
                            aux++;
                        }
                    }
                    sizeRows[j] = aux;
                    aux = 0;
                }
//            printf("antes\n");
                for(ite=0; ite<sizeRows[0]; ite++)
                {
                    el = mat[0][ite];
                    if(Coef1_mult[el]<0)
                    {

                        aux = Coef1_mult[el] - 1;
                    }
                    else
                    {
                        aux = Coef1_mult[el];
                    }
//                printf("Coef %LF aux %d  k %d  j %d el %d\n",Coef1_mult[el], aux,k,0,el);
                    Coef_ref[el] = aux*(k+1);

                }
//            printf("depois\n");
                for(j=1; j<=k; j++)
                {
                    for(ite = 0; ite<sizeRows[j]; ite++)
                    {
                        el = mat[j][ite];

                        if(Coef1_mult[el]<0)
                        {
                            aux = Coef1_mult[el] - 1;
                        }
                        else
                        {
                            aux = Coef1_mult[el];
                        }
//                    printf("Coef %LF aux %d  k %d  j %d el %d\n",Coef1_mult[el], aux,k,j,el);
                        Coef_ref[el] = aux*(k+1) + j;
                    }

                }
                rhs_ref = (k+1)*rhs;
                lhs = 0;

                for(ite = 0; ite<h_cut->numberVariables; ite++)
                {
                    lhs += (long double)Coef_ref[ite]*(long double)h_cut->xAsterisc[ite];
                }
                violation_mult[c_aux] = lhs - (rhs_ref*precision);

                for (ite = 0; ite < k+1; ite++)
                {
                    free (mat[ite]);
                }
                free (mat);
                free(sizeRows);
                tam = 0;
                if(verifyDominanceCG(Coef_ref,rhs_ref,Coef1,rhs,h_cut->numberVariables)==1)
                {
                    for(j=0; j<h_cut->numberVariables; j++)
                    {
                        if(Coef_ref[j]!=0)
                        {
                            v_aux[tam] = Coef_ref[j];
                            tam++;
                        }
                    }
                    v_aux[tam] = rhs_ref;
                    //tam++;
                    int mdc = CutP_maxDivisorCommonVector(v_aux, tam);
                    for(j=0; j<h_cut->numberVariables; j++)
                    {
                        if(Coef_ref[j]!=0)
                        {
                            Coefs_temp[cont_aux] = Coef_ref[j]/mdc;
                            //printf("%d\t",Coefs_temp[cont_aux]);
                            Elements_temp[cont_aux] = j;
                            cont_aux++;
                        }
                    }
                    Pos_el_temp[c_aux+1] = cont_aux;
                    rhs_temp[c_aux] = rhs_ref/mdc;

                }
                else
                {
                    for(j=0; j<h_cut->numberVariables; j++)
                    {
                        if(Coef1[j]!=0)
                        {
                            v_aux[tam] = Coef1[j];
                            tam++;
                        }
                    }
                    v_aux[tam] = rhs;
                    //tam++;
                    int mdc = CutP_maxDivisorCommonVector(v_aux, tam);
                    for(j=0; j<h_cut->numberVariables; j++)
                    {
                        if(Coef1[j]!=0)
                        {
                            Coefs_temp[cont_aux] = Coef1[j]/mdc;
                            //printf("%d\t",Coefs_temp[cont_aux]);
                            Elements_temp[cont_aux] = j;
                            cont_aux++;
                        }
                    }
                    Pos_el_temp[c_aux+1] = cont_aux;
                    rhs_temp[c_aux] = rhs/mdc;
                }
//                if(vio_aux>violation[c_aux])
//                    exit(0);
            }
            else
            {
                tam = 0;
                for(j=0; j<h_cut->numberVariables; j++)
                {
                    if(Coef1[j]!=0)
                    {
                        v_aux[tam] = Coef1[j];
                        tam++;
                    }
                }
                v_aux[tam] = rhs;
                //tam++;
                int mdc = CutP_maxDivisorCommonVector(v_aux, tam);
                for(j=0; j<h_cut->numberVariables; j++)
                {
                    if(Coef1[j]!=0)
                    {
                        Coefs_temp[cont_aux] = Coef1[j]/mdc;
                        //printf("%d\t",Coefs_temp[cont_aux]);
                        Elements_temp[cont_aux] = j;
                        cont_aux++;
                    }
                }
                Pos_el_temp[c_aux+1] = cont_aux;
                rhs_temp[c_aux] = rhs/mdc;
            }
            c_aux++;
        }
    }
    free(v_sol_bin);
    free(Coef1);
    free(Coef1_mult);
    free(Coef_ref);
    free(violation);
    free(violation_mult);
    free(p_frac);
    free(v_aux);

    cuts_generated = AllocationStructCut(cont_ini + cont_aux, nC_initial + nCuts,h_cut->numberVariables);
    for(i=0; i<cuts_generated->numberVariables; i++)
    {
        cuts_generated->xAsterisc[i] = h_cut->xAsterisc[i];
    }
    cuts_generated->ElementsConstraints[0] = 0;
    for(i=0; i<nC_initial; i++)
    {
        cuts_generated->rightSide[i] = h_cut->rightSide[i];
        cuts_generated->typeConstraints[i] = h_cut->typeConstraints[i];
        cuts_generated->ElementsConstraints[i+1] = h_cut->ElementsConstraints[i+1];
    }
    aux = 1;

    for(i = nC_initial; i<cuts_generated->numberConstrains; i++)
    {
        cuts_generated->rightSide[i] = rhs_temp[i - nC_initial];
        cuts_generated->typeConstraints[i] = LPC_CGGPU;
        cuts_generated->ElementsConstraints[i+1] = Pos_el_temp[aux] + cont_ini;
        aux++;//+ h_cut->cont;
    }
    for(i=0; i<cont_ini; i++)
    {
        cuts_generated->Coefficients[i] = h_cut->Coefficients[i];
        cuts_generated->Elements[i] = h_cut->Elements[i];
    }
    for(i= cont_ini; i<cuts_generated->cont; i++)
    {
        cuts_generated->Coefficients[i] = Coefs_temp[i - cont_ini];
        cuts_generated->Elements[i] = Elements_temp[i - cont_ini];
    }


    int *validated = (int*)malloc(sizeof(int)*cuts_generated->numberConstrains);
    memset(validated,0,sizeof(int)*cuts_generated->numberConstrains);
    aux = 1;
    cont_aux=0;
    int p1,p2, minus_elements = 0;
    k=0;
    for(i=0; i<cuts_generated->numberConstrains - 1; i++)
    {
        for(j=i+1; j<cuts_generated->numberConstrains; j++)
        {
            if(validated[j]==0)
            {
                if((cuts_generated->ElementsConstraints[i+1]-cuts_generated->ElementsConstraints[i])!=(cuts_generated->ElementsConstraints[j+1]-cuts_generated->ElementsConstraints[j])||(cuts_generated->rightSide[i]!=cuts_generated->rightSide[j]))
                {
                    aux=0;
                }
                else
                {
                    p1 = cuts_generated->ElementsConstraints[i];
                    p2 = cuts_generated->ElementsConstraints[j];
                    aux = 1;
                    for(k=0; k<cuts_generated->ElementsConstraints[i+1]-cuts_generated->ElementsConstraints[i]; k++)
                    {

                        if((cuts_generated->Coefficients[p1+k]!=cuts_generated->Coefficients[p2+k])||(cuts_generated->Elements[p1+k]!=cuts_generated->Elements[p2+k]))
                        {
                            aux=0;

                            break;
                        }
                    }
                }
                if((aux==1)&&(j >= h_cut->numberConstrains))
                {
                    validated[j]=1;
                    minus_elements += cuts_generated->ElementsConstraints[j+1]-cuts_generated->ElementsConstraints[j];
                    cont_aux++;
                }
            }
        }
    }

//    printf("Number of repeat Strengthening: %d \n", cont_aux);

    Cut_gpu* new_h_cut;

    new_h_cut = AllocationStructCut(cuts_generated->cont - minus_elements,cuts_generated->numberConstrains-cont_aux,cuts_generated->numberVariables);
    aux = 0;
    cont_aux = 0;
    new_h_cut->ElementsConstraints[0] = 0;
    for(i=0; i<cuts_generated->numberConstrains; i++)
    {
        if(validated[i]==0)
        {
//            if(i>=h_cut->numberConstrains)
//            {
//                printf("violation %f\n",violation[i-h_cut->numberConstrains]);
//            }
            new_h_cut->rightSide[aux] = cuts_generated->rightSide[i];
            new_h_cut->typeConstraints[aux] = cuts_generated->typeConstraints[i];
            for(j = cuts_generated->ElementsConstraints[i]; j < cuts_generated->ElementsConstraints[i+1]; j++)
            {
                new_h_cut->Coefficients[cont_aux] = cuts_generated->Coefficients[j];
                new_h_cut->Elements[cont_aux] = cuts_generated->Elements[j];
                cont_aux++;
            }
            new_h_cut->ElementsConstraints[aux + 1] = cont_aux;

            aux++;
        }

    }
    for(i=0; i<new_h_cut->numberVariables; i++)
    {
        new_h_cut->xAsterisc[i] = cuts_generated->xAsterisc[i];
    }
//
//    printf("End phase 1!\n");
    free(Coefs_temp);
    free(Elements_temp);
    free(Pos_el_temp);
    free(rhs_temp);
    free(validated);
    free(cuts_generated);

    return new_h_cut;





}

Cut_gpu* createCutsStrongPhaseTwo(Cut_gpu *h_cut, solutionGpu *h_solution, int numberMaxConst, int nCuts, int precision, int nRuns, int nThreads,int nBlocks)
{
//    double value = 0 ;
    double *value_violation = (double*)malloc(sizeof(double)*nCuts);
    double *value_violation_r1 = (double*)malloc(sizeof(double)*nCuts);
    double *violation_mult = (double*)malloc(sizeof(double)*nCuts);
    int *Coef1 = (int*)malloc(sizeof(int)*(h_cut->numberVariables));
    int *Coef2 = (int*)malloc(sizeof(int)*(h_cut->numberVariables));
    int *Coef_ref = (int*)malloc(sizeof(int)*(h_cut->numberVariables));
    long double *Coef_mult1 = (long double*)malloc(sizeof(long double)*h_cut->numberVariables);
    long double *p_frac = (long double*)malloc(sizeof(long double)*h_cut->numberVariables);
    int *Coefs_temp = (int*)malloc(sizeof(int)*(nCuts*h_cut->numberVariables));
    int *Elements_temp = (int*)malloc(sizeof(int)*(nCuts*h_cut->numberVariables));
    int *Pos_el_temp = (int*)malloc(sizeof(int)*(nCuts+1) );
    int *rhs_temp = (int*)malloc(sizeof(int)*(nCuts));
    int *v_aux = (int*)malloc(sizeof(int)*(h_cut->numberVariables+1));
    long double rhs_mult1, rhs_mult2;
//    int nAll = 0;
    int temp, k_R2;
    int i,k,j,n_1,d_1,n_2,d_2, qnt_1, el, rest_a, rhs1 = 0, rhs2 = 0,aux = 0, cont_aux=0,ite,tam=0;
//    int n1_aux,  n2_aux;
    Cut_gpu *cuts_generated;
//     printf("Initial phase 2! %d\n",nCuts);
    Pos_el_temp[0] = 0;
    for (i=0; i<nThreads; i++)
    {
        for(k=0; k<nBlocks; k++)
        {
            if(h_solution->SConst[0 + i*numberMaxConst + k*numberMaxConst*nThreads]!=-1)
            {
                for(ite=0; ite<h_cut->numberVariables; ite++)
                {
                    Coef_mult1[ite] = 0;
                    p_frac[ite] = 0;
                }

                rhs1 = 0;
                rhs2 = 0;
                rhs_mult1 = 0;
                rhs_mult2 = 0;
                qnt_1 = h_solution->SPAux[i + k*nThreads];
                n_1 = h_solution->SMult[0 + i*4 + k*4*nThreads];
                d_1 = h_solution->SMult[1 + i*4 + k*4*nThreads];
                n_2 = h_solution->SMult[2 + i*4 + k*4*nThreads];
                d_2 = h_solution->SMult[3 + i*4 + k*4*nThreads];
                //printf("denominator %d %d %d %d %d\n",n_1, d_1,n_2,d_2,h_solution->SConst[0 + i*numberMaxConst + k*numberMaxConst*nThreads]);
                if((d_1==0)||(d_2==0))
                    continue;

                for(rest_a = 0; rest_a <= qnt_1; rest_a++)
                {
                    rhs1 += h_cut->rightSide[ h_solution->SConst[rest_a+ i*numberMaxConst + k*numberMaxConst*nThreads] ] * n_1;
                }

                for(rest_a = qnt_1 + 1; rest_a < numberMaxConst; rest_a++)
                {
                    rhs2 += h_cut->rightSide[ h_solution->SConst[rest_a+ i*numberMaxConst + k*numberMaxConst*nThreads] ] * n_2;
                }
                rhs_mult1 = ((double)rhs1/(double)d_1 + (double)rhs2/(double)d_2);// RHS Double
                if ((rhs1*d_2 + rhs2*d_1)/(d_1*d_2) < 0)
                {
                    temp = (rhs1*d_2 + rhs2*d_1)/(d_1*d_2) - 1;
                }
                else
                {
                    temp = (rhs1*d_2 + rhs2*d_1)/(d_1*d_2);
                }
                rhs_mult2 = rhs_mult1 - temp;
                //printf("rhs_mult: %LF\n",rhs_mult2);
//                if((rhs_mult2 < 0.5)&&(rhs_mult2>0.01))
//                {
//                    int t_mult = 1;
//                    int aux_temp;
//                    do
//                    {
//                        t_mult++;
//                    }
//                    while(t_mult*rhs_mult2<0.5);
//                    n1_aux = n_1;
//                    n2_aux = n_2;
//                     //printf("f_rhs: %LF, t: %d, n1: %d d1: %d n2: %d d2: %d \n",rhs_mult2,t_mult, n_1,d_1,n_2,d_2);
//                    aux_temp = (n_1*t_mult)/d_1;
//                    n_1 = n_1*t_mult - aux_temp*d_1;
//                    aux_temp = (n_2*t_mult)/d_2;
//                    n_2 = n_2*t_mult - aux_temp*d_2;
//                    //printf("f_rhs: %LF, t: %d, n1: %d d1: %d n2: %d d2: %d \n",rhs_mult2,t_mult, n_1,d_1,n_2,d_2);
//                    //getchar();
//
//                }
                //printf("n1: %d d1: %d n2: %d d2: %d\n",n_1,d_1,n_2,d_2);
                rhs1 = 0;
                rhs2 = 0;
                rhs_mult1 = 0;
                rhs_mult2 = 0;
                memset(Coef1,0,sizeof(int)*h_cut->numberVariables);
                memset(Coef2,0,sizeof(int)*h_cut->numberVariables);
                memset(Coef_ref,0,sizeof(int)*h_cut->numberVariables);
                for(rest_a = 0; rest_a <= qnt_1; rest_a++)
                {
                    for(j = h_cut->ElementsConstraints[ h_solution->SConst[rest_a + i*numberMaxConst + k*numberMaxConst*nThreads] ]; j< h_cut->ElementsConstraints[ h_solution->SConst[rest_a + i*numberMaxConst + k*numberMaxConst*nThreads] + 1]; j++)
                    {
                        el = h_cut->Elements[j];
                        Coef1[el] += h_cut->Coefficients[j] * n_1;
                        //printf("Coef1 %d: %d %d\t",el, h_cut->Coefficients[j], h_cut->xAsterisc[el]);
                        //value_tes += Coef[el] * d_cut->xAsterisc[el];
                    }
                    //printf("rhs1: %d\n",h_cut->rightSide[ h_solution->SConst[rest_a+ i*numberMaxConst + k*numberMaxConst*nThreads] ]);
                    //printf("tipo: %d \n", h_cut->typeConstraints[h_solution->SConst[rest_a+ i*numberMaxConst + k*numberMaxConst*nThreads]]);
                    rhs1 += h_cut->rightSide[ h_solution->SConst[rest_a+ i*numberMaxConst + k*numberMaxConst*nThreads] ] * n_1;

                }


                for(rest_a = qnt_1 + 1; rest_a < numberMaxConst; rest_a++)
                {
                    for(j=h_cut->ElementsConstraints[ h_solution->SConst[rest_a+ i*numberMaxConst + k*numberMaxConst*nThreads] ]; j<h_cut->ElementsConstraints[  h_solution->SConst[rest_a + i*numberMaxConst + k*numberMaxConst*nThreads] + 1]; j++)
                    {
                        el = h_cut->Elements[j];
                        Coef2[el] += h_cut->Coefficients[j] * n_2;
                        //printf("Coef2 %d: %d %d\t",el, h_cut->Coefficients[j], h_cut->xAsterisc[el]);
                        //value_tes += Coef[el] * d_cut->xAsterisc[el];
                    }
                    //printf("rhs2: %d\n",h_cut->rightSide[ h_solution->SConst[rest_a+ i*numberMaxConst + k*numberMaxConst*nThreads] ]);
                    //printf("tipo: %d \n", h_cut->typeConstraints[h_solution->SConst[rest_a+ i*numberMaxConst + k*numberMaxConst*nThreads]]);
                    rhs2 += h_cut->rightSide[ h_solution->SConst[rest_a+ i*numberMaxConst + k*numberMaxConst*nThreads] ] * n_2;
                }

                value_violation[cont_aux] = 0;
                value_violation_r1[cont_aux] = 0;
                violation_mult[cont_aux] = 0;

                for(j=0; j<h_cut->numberVariables; j++)
                {
                    //temp = (double)Coef1[j]/(double)d_1 + (double)Coef2[j]/(double)d_2;
                    Coef_mult1[j] = (long double)Coef1[j]/(long double)d_1 + (long double)Coef2[j]/(long double)d_2; //Coef Double
                    //printf("Coef_mult %LF, %d %d %d %d\n",Coef_mult1[j], Coef1[j] , d_1, Coef2[j],d_2);
                    if( (Coef1[j]*d_2+ Coef2[j]*d_1)/(d_1*d_2)<0 )
                    {
                        temp = ( (Coef1[j]*d_2+ Coef2[j]*d_1)/(d_1*d_2) ) -1;
                    }
                    else
                    {
                        temp = (Coef1[j]*d_2+ Coef2[j]*d_1)/(d_1*d_2);
                    }

                    p_frac[j] = Coef_mult1[j] - temp; //Coef Frac
                    //  printf("Coef_mult1 = %LF Temp %d p_frac %LF\n", Coef_mult1[j],temp,p_frac[j]);
                    if(temp != 0)
                    {
                        Coef1[j] = temp;
                        value_violation[cont_aux] += (temp)*h_cut->xAsterisc[j];

                    }
                    else
                    {
                        Coef1[j] = 0;
                    }


                }
                //temp = ((double)rhs1/(double)d_1 + (double)rhs2/(double)d_2);

                rhs_mult1 = ((double)rhs1/(double)d_1 + (double)rhs2/(double)d_2);// RHS Double
                if ((rhs1*d_2 + rhs2*d_1)/(d_1*d_2) < 0)
                {
                    temp = (rhs1*d_2 + rhs2*d_1)/(d_1*d_2) - 1;
                }
                else
                {
                    temp = (rhs1*d_2 + rhs2*d_1)/(d_1*d_2);
                }
                rhs_temp[cont_aux] = temp;
                value_violation[cont_aux] -= (rhs_temp[cont_aux])*precision;




                rhs_mult2 = rhs_mult1 - temp; //Frac RHS
                //printf(" RHS_post: %d \n", rhs_temp[cont_aux]);
                //for(ite=0; ite<h_cut->numberVariables; ite++)
                //{
                //    if((p_frac[ite] - rhs_mult2)/(1.0 - rhs_mult2)>0)
                //    {
                //        Coef_mult1[ite] = Coef1[ite] + (p_frac[ite] - rhs_mult2)/(1.0 - rhs_mult2);
                //    }
                //    else
                //    {
                //        Coef_mult1[ite] = Coef1[ite];
                //    }
                //    value_violation_r1[cont_aux] += (Coef_mult1[ite])*h_cut->xAsterisc[ite];
                //}


                //value_violation_r1[cont_aux] -= rhs_temp[cont_aux]*precision;
                //printf("Violation Initial: %lf\n",value_violation[cont_aux]);
//                printf("Violation reforço 1: %lf\n",value_violation_r1[cont_aux]);
//
//                for(ite = 0; ite<h_cut->numberVariables; ite++ )
//                {
//                    if(Coef1[ite]!=0)
//                    {
//                        printf("%d v_%d \t", Coef1[ite], ite );
//                    }
//                }
//                printf("<= %d\n", rhs_temp[cont_aux]);
//                printf("reforço 1\n");
//                for(ite = 0; ite<h_cut->numberVariables; ite++ )
//                {
//                    if(Coef_mult1[ite]!=0)
//                    {
//                        printf("%LF v_%d \t", Coef_mult1[ite], ite );
//                    }
//                }
//                printf("<= %d\n", rhs_temp[cont_aux]);
                rhs1 = rhs_temp[cont_aux];
                //printf("rhs1: %d\n",rhs1);
                if(rhs_mult2>0)
                {
                    k_R2 = returnK(rhs_mult2);

                    assert( rhs_mult2<=(1.0/((long double)k_R2)) + 1e-8 );
                    assert( rhs_mult2>=(1.0/((long double)k_R2+1)) - 1e-8 );

                    //if(value_violation[cont_aux]>=precision){
                    //    printf("Value violation: %f %d\n", value_violation[cont_aux],precision);
                    //    for(rest_a = 0; rest_a< numberMaxConst; rest_a++){
                    //        int testando = h_solution->SConst[rest_a+ i*numberMaxConst + k*numberMaxConst*nThreads];
                    //        printf("%d \t", testando);
                    //        printf("Tipo: %d \n",h_cut->typeConstraints[testando]);
                    //    }
                    //    printf("\n");
                    //:}
//                    printf("K: %d\n",k_R2);
                    int *sizeRows = (int*)malloc((k_R2+1)*sizeof(int));
                    memset(sizeRows,0,(k_R2+1)*sizeof(int));
                    int **mat = (int**)malloc( (k_R2+1)*sizeof(int*));
                    for(ite = 0; ite < (k_R2+1); ite++)
                        mat[ite] = (int*)malloc(h_cut->numberVariables*sizeof(int));
                    temp = 0;
                    for(ite=0; ite<h_cut->numberVariables; ite++)
                    {
                        if((p_frac[ite] - 1e-6 <=rhs_mult2)&&(Coef_mult1[ite]!=0))
                        {
                            mat[0][temp] = ite;
                            temp++;
                        }
                    }
                    sizeRows[0] = temp;
                    temp = 0;
                    for(j=1; j<=k_R2; j++)
                    {
                        for(ite=0; ite<h_cut->numberVariables; ite++)
                        {
                            if(Coef_mult1[ite]==0)
                            {
                                continue;
                            }

                            if((p_frac[ite] - 1e-6 > (rhs_mult2 + (j-1.0)*(1.0-rhs_mult2)/(long double)k_R2))&&(p_frac[ite] + 1e-6 <=(rhs_mult2 + (j*(1.0-rhs_mult2)/(long double)k_R2))))
                            {
                                mat[j][temp] = ite;
                                temp++;
                            }
                        }
                        sizeRows[j] = temp;
                        temp = 0;
                    }
//            printf("antes\n");
                    for(ite=0; ite<sizeRows[0]; ite++)
                    {
                        el = mat[0][ite];
                        if(Coef_mult1[el]<0)
                        {

                            temp = Coef_mult1[el] - 1;
                        }
                        else
                        {
                            temp = Coef_mult1[el];
                        }
//                printf("Coef %LF aux %d  k %d  j %d el %d\n",Coef1_mult[el], aux,k,0,el);
                        Coef_ref[el] = temp*(k_R2+1);

                    }
//            printf("depois\n");
                    for(j=1; j<=k_R2; j++)
                    {
                        for(ite = 0; ite<sizeRows[j]; ite++)
                        {
                            el = mat[j][ite];

                            if(Coef_mult1[el]<0)
                            {
                                temp = Coef_mult1[el] - 1;
                            }
                            else
                            {
                                temp = Coef_mult1[el];
                            }
//                    printf("Coef %LF aux %d  k %d  j %d el %d\n",Coef1_mult[el], aux,k,j,el);
                            Coef_ref[el] = temp*(k_R2+1) + j;
                        }

                    }
                    rhs_temp[cont_aux] = (k_R2+1)*rhs_temp[cont_aux];
                    int lhs = 0;

                    for(ite = 0; ite<h_cut->numberVariables; ite++)
                    {
                        lhs += (long double)Coef_ref[ite]*(long double)h_cut->xAsterisc[ite];
                    }
                    violation_mult[cont_aux] = lhs - (rhs_temp[cont_aux]*precision);
//                    printf("Violation 2: %lf \n",violation_mult[cont_aux]);

//                    printf("reforço 2\n");
//                    for(ite = 0; ite<h_cut->numberVariables; ite++ )
//                    {
//                        if(Coef_ref[ite]!=0)
//                        {
//                            printf("%d v_%d \t", Coef_ref[ite], ite );
//                        }
//                    }
//                    printf("<= %d\n", rhs_temp[cont_aux]);
//                    for(ite = 0; ite<h_cut->numberVariables; ite++ )
//                    {
//                        printf("v_%d = %d\n", ite, h_cut->xAsterisc[ite]);
//                    }


                    for (ite = 0; ite < k_R2+1; ite++)
                    {
                        free (mat[ite]);
                    }
                    free (mat);
                    free(sizeRows);
                    tam = 0;
                    if(verifyDominanceCG(Coef_ref,rhs_temp[cont_aux],Coef1,rhs1,h_cut->numberVariables)==1)
                    {
                        for(j=0; j<h_cut->numberVariables; j++)
                        {
                            if(Coef_ref[j]!=0)
                            {
                                v_aux[tam] = Coef_ref[j];
                                tam++;
                            }
                        }
                        v_aux[tam] = rhs_temp[cont_aux];
                        //tam++;

                        int mdc = CutP_maxDivisorCommonVector(v_aux, tam);
                        for(j=0; j<h_cut->numberVariables; j++)
                        {
                            if(Coef_ref[j]!=0)
                            {
                                Coefs_temp[aux] = Coef_ref[j]/mdc;
                                //printf("%d\t",Coefs_temp[cont_aux]);
                                Elements_temp[aux] = j;
                                aux++;
                            }
                        }
                        Pos_el_temp[cont_aux+1] = aux;
                        rhs_temp[cont_aux] = rhs_temp[cont_aux]/mdc;

                    }
                    else
                    {
                        tam = 0;
                        for(j=0; j<h_cut->numberVariables; j++)
                        {
                            if(Coef1[j]!=0)
                            {
                                v_aux[tam] = Coef1[j];
                                tam++;
                            }
                        }
                        v_aux[tam] = rhs1;
                        // tam++;

                        int mdc = CutP_maxDivisorCommonVector(v_aux, tam);
                        for(j=0; j<h_cut->numberVariables; j++)
                        {
                            if(Coef1[j]!=0)
                            {
                                Coefs_temp[aux] = Coef1[j]/mdc;
                                //printf("%d\t",Coefs_temp[cont_aux]);
                                Elements_temp[aux] = j;
                                aux++;
                            }
                        }
                        Pos_el_temp[cont_aux+1] = aux;
                        rhs_temp[cont_aux] = rhs1/mdc;
                    }
//                    if(value_violation_r1[cont_aux]>value_violation[cont_aux])
//                    {
//                        exit(0);
//                    }


                }
                else
                {
                    tam = 0;
                    for(j=0; j<h_cut->numberVariables; j++)
                    {
                        if(Coef1[j]!=0)
                        {
                            v_aux[tam] = Coef1[j];
                            tam++;
                        }
                    }
                    v_aux[tam] = rhs1;
                    //tam++;

                    int mdc = CutP_maxDivisorCommonVector(v_aux, tam);
                    for(j=0; j<h_cut->numberVariables; j++)
                    {
                        if(Coef1[j]!=0)
                        {
                            Coefs_temp[aux] = Coef1[j]/mdc;
                            //printf("%d\t",Coefs_temp[cont_aux]);
                            Elements_temp[aux] = j;
                            aux++;
                        }
                    }
                    Pos_el_temp[cont_aux+1] = aux;
                    rhs_temp[cont_aux] = rhs1/mdc;
                }
                cont_aux++;

            }
        }
    }

    free(Coef1);
    free(Coef2);
    free(Coef_mult1);
    free(Coef_ref);
    free(violation_mult);
    free(value_violation);
    free(value_violation_r1);
    free(v_aux);
    free(p_frac);
    cont_aux = aux;
    cuts_generated = AllocationStructCut(h_cut->cont+cont_aux, h_cut->numberConstrains + nCuts,h_cut->numberVariables);
//    printf("Intermediry teste \n");

    for(i=0; i<cuts_generated->numberVariables; i++)
    {
        cuts_generated->xAsterisc[i] = h_cut->xAsterisc[i];
    }
    cuts_generated->ElementsConstraints[0] = 0;
    for(i=0; i<h_cut->numberConstrains; i++)
    {
        cuts_generated->rightSide[i] = h_cut->rightSide[i];
        cuts_generated->typeConstraints[i] = h_cut->typeConstraints[i];
        cuts_generated->ElementsConstraints[i+1] = h_cut->ElementsConstraints[i+1];
    }
    aux = 1;

    for(i = h_cut->numberConstrains; i<cuts_generated->numberConstrains; i++)
    {

        cuts_generated->rightSide[i] = rhs_temp[i - h_cut->numberConstrains];
        cuts_generated->typeConstraints[i] = LPC_CGGPU;
        cuts_generated->ElementsConstraints[i+1] = Pos_el_temp[aux] + h_cut->cont;
        aux++;//+ h_cut->cont;
    }
    for(i=0; i<h_cut->cont; i++)
    {
        cuts_generated->Coefficients[i] = h_cut->Coefficients[i];
        cuts_generated->Elements[i] = h_cut->Elements[i];
    }
    for(i= h_cut->cont; i<cuts_generated->cont; i++)
    {
        cuts_generated->Coefficients[i] = Coefs_temp[i - h_cut->cont];
        cuts_generated->Elements[i] = Elements_temp[i - h_cut->cont];
    }


    int *validated = (int*)malloc(sizeof(int)*cuts_generated->numberConstrains);
    memset(validated,0,sizeof(int)*cuts_generated->numberConstrains);
    aux = 1;
    cont_aux=0;
    int p1,p2, minus_elements = 0;
    k=0;
    for(i=0; i<cuts_generated->numberConstrains - 1; i++)
    {
        for(j=i+1; j<cuts_generated->numberConstrains; j++)
        {
            if(validated[j]==0)
            {
                if((cuts_generated->ElementsConstraints[i+1]-cuts_generated->ElementsConstraints[i])!=(cuts_generated->ElementsConstraints[j+1]-cuts_generated->ElementsConstraints[j])||(cuts_generated->rightSide[i]!=cuts_generated->rightSide[j]))
                {
                    aux=0;
                }
                else
                {
                    p1 = cuts_generated->ElementsConstraints[i];
                    p2 = cuts_generated->ElementsConstraints[j];
                    aux = 1;
                    for(k=0; k<cuts_generated->ElementsConstraints[i+1]-cuts_generated->ElementsConstraints[i]; k++)
                    {

                        if((cuts_generated->Coefficients[p1+k]!=cuts_generated->Coefficients[p2+k])||(cuts_generated->Elements[p1+k]!=cuts_generated->Elements[p2+k]))
                        {
                            aux=0;

                            break;
                        }
                    }
                }
                if((aux==1)&&(j >= h_cut->numberConstrains))
                {
                    validated[j]=1;
                    minus_elements += cuts_generated->ElementsConstraints[j+1]-cuts_generated->ElementsConstraints[j];
                    cont_aux++;
                }
            }
        }
    }

//    printf("Number of repeat Strengthening: %d \n", cont_aux);

    Cut_gpu* new_h_cut;

    new_h_cut = AllocationStructCut(cuts_generated->cont - minus_elements,cuts_generated->numberConstrains-cont_aux,cuts_generated->numberVariables);
    aux = 0;
    cont_aux = 0;
    new_h_cut->ElementsConstraints[0] = 0;
    for(i=0; i<cuts_generated->numberConstrains; i++)
    {
        if(validated[i]==0)
        {
//            if(i>=h_cut->numberConstrains)
//            {
//                printf("violation %f\n",violation[i-h_cut->numberConstrains]);
//            }
            new_h_cut->rightSide[aux] = cuts_generated->rightSide[i];
            new_h_cut->typeConstraints[aux] = cuts_generated->typeConstraints[i];
            for(j = cuts_generated->ElementsConstraints[i]; j < cuts_generated->ElementsConstraints[i+1]; j++)
            {
                new_h_cut->Coefficients[cont_aux] = cuts_generated->Coefficients[j];
                new_h_cut->Elements[cont_aux] = cuts_generated->Elements[j];
                cont_aux++;
            }
            new_h_cut->ElementsConstraints[aux + 1] = cont_aux;

            aux++;
        }

    }
    for(i=0; i<new_h_cut->numberVariables; i++)
    {
        new_h_cut->xAsterisc[i] = cuts_generated->xAsterisc[i];
    }
//    printf("End phase 2!\n");
    free(cuts_generated);
    free(validated);
    free(Coefs_temp);
    free(Elements_temp);
    free(Pos_el_temp);
    free(rhs_temp);

    return new_h_cut;

}



Cut_gpu* createCutsOfPhaseOne(Cut_gpu *h_cut, solutionGpu *h_solution, int nCuts, int precision, int nRuns)
{
    int i,j,el,aux,constraint, rhs, mdc, tam = 0,lhs = 0, n1,d1 ;

    int *Coef1 = (int*)malloc(sizeof(int)*(h_cut->numberVariables)); //free
    int *v_aux = (int*)malloc(sizeof(int)*(h_cut->numberVariables + 1)); //free
    Cut_gpu *cuts_generated; //free

    int *Coefs_temp = (int*)malloc(sizeof(int)*(nCuts*h_cut->numberVariables)); //free
    int *Elements_temp = (int*)malloc(sizeof(int)*(nCuts*h_cut->numberVariables)); //free
    int *Pos_el_temp = (int*)malloc(sizeof(int)*(nCuts+1) ); //free
    int *rhs_temp = (int*)malloc(sizeof(int)*(nCuts)); //free
    int cont_aux = 0, c_aux = 0;
    memset(v_aux,0,sizeof(int)*(h_cut->numberVariables+1));
    // double violation_media = 0;
    // double violation_max = 0;
    //  double violation_min = 100000;
    // double *violation = (double*)malloc(sizeof(double)*nCuts);

    Pos_el_temp[0] = 0;
    for(i = 0; i<nRuns; i++)
    {
        if(h_solution->SConst[i]!=-1)
        {
            aux = 0;
            rhs = 0;
            tam = 0;
            lhs = 0 ;
            //  violation[c_aux] = 0;
            memset(Coef1,0,sizeof(int)*h_cut->numberVariables);
            constraint = h_solution->SConst[i];
            n1 = h_solution->SMult[i];
            d1 = h_solution->SPAux[i];
            for(j = h_cut->ElementsConstraints[constraint]; j< h_cut->ElementsConstraints[constraint + 1]; j++)
            {
                el = h_cut->Elements[j];
                Coef1[el] = h_cut->Coefficients[j] * n1;
                if((Coef1[el]<0&&d1>0) || (Coef1[el]>0&&d1<0))
                {
                    Coef1[el] = (Coef1[el]/d1)-1;
                }
                else
                {
                    Coef1[el] = (Coef1[el]/d1);
                }
                lhs += Coef1[el]*h_cut->xAsterisc[el];
                if(Coef1[el]!=0)
                {
                    v_aux[tam] = Coef1[el];
                    tam++;
                }
            }
            rhs = h_cut->rightSide[constraint] * n1;
            if((rhs<0&&d1>0) || (rhs>0&&d1<0))
            {
                rhs = (rhs/d1)-1;
            }
            else
            {
                rhs = rhs/d1;
            }

//            violation[c_aux] = lhs - (rhs*precision);
            //violation = violation - (rhs*1000);
            //printf("Violation:%f\n",violation);
            v_aux[tam] = rhs;
            //tam++;
            //for(j=0;j<tam;j++){
            //    printf("%d \t", v_aux[j]);
            //}
            //printf("\n");
            mdc = CutP_maxDivisorCommonVector(v_aux, tam);
            for(j=0; j<h_cut->numberVariables; j++)
            {
                if(Coef1[j]!=0)
                {
                    Coefs_temp[cont_aux] = Coef1[j]/mdc;
                    //printf("%d\t",Coefs_temp[cont_aux]);
                    Elements_temp[cont_aux] = j;
                    cont_aux++;
                }
            }
            // getchar();
            // printf("\n %d\n",cont_aux);
            Pos_el_temp[c_aux+1] = cont_aux;
            rhs_temp[c_aux] = rhs/mdc;
            c_aux++;
            //printf("%f\n",violation);
        }

    }
    cuts_generated = AllocationStructCut(h_cut->cont+cont_aux, h_cut->numberConstrains + nCuts,h_cut->numberVariables);
    for(i=0; i<cuts_generated->numberVariables; i++)
    {
        cuts_generated->xAsterisc[i] = h_cut->xAsterisc[i];
    }
    cuts_generated->ElementsConstraints[0] = 0;
    for(i=0; i<h_cut->numberConstrains; i++)
    {
        cuts_generated->rightSide[i] = h_cut->rightSide[i];
        cuts_generated->typeConstraints[i] = h_cut->typeConstraints[i];
        cuts_generated->ElementsConstraints[i+1] = h_cut->ElementsConstraints[i+1];
    }
    aux = 1;
    for(i = h_cut->numberConstrains; i<cuts_generated->numberConstrains; i++)
    {

        cuts_generated->rightSide[i] = rhs_temp[i - h_cut->numberConstrains];
        cuts_generated->typeConstraints[i] = LPC_CGGPU;
        cuts_generated->ElementsConstraints[i+1] = Pos_el_temp[aux] + h_cut->cont;
        aux++;//+ h_cut->cont;
    }
    for(i=0; i<h_cut->cont; i++)
    {
        cuts_generated->Coefficients[i] = h_cut->Coefficients[i];
        cuts_generated->Elements[i] = h_cut->Elements[i];
    }
    for(i= h_cut->cont; i<cuts_generated->cont; i++)
    {
        cuts_generated->Coefficients[i] = Coefs_temp[i - h_cut->cont];
        cuts_generated->Elements[i] = Elements_temp[i - h_cut->cont];
    }

    int *validated = (int*)malloc(sizeof(int)*cuts_generated->numberConstrains); //free
    memset(validated,0,sizeof(int)*cuts_generated->numberConstrains);
    aux = 1;
    cont_aux=0;
    int k,p1,p2, minus_elements = 0;
    for(i=0; i<cuts_generated->numberConstrains - 1; i++)
    {

        for(j=i+1; j<cuts_generated->numberConstrains; j++)
        {
            if(validated[j]==0)
            {
                if((cuts_generated->ElementsConstraints[i+1]-cuts_generated->ElementsConstraints[i])!=(cuts_generated->ElementsConstraints[j+1]-cuts_generated->ElementsConstraints[j])||(cuts_generated->rightSide[i]!=cuts_generated->rightSide[j]))
                {
                    aux=0;
                }
                else
                {
                    p1 = cuts_generated->ElementsConstraints[i];
                    p2 = cuts_generated->ElementsConstraints[j];
                    aux = 1;
                    for(k=0; k<cuts_generated->ElementsConstraints[i+1]-cuts_generated->ElementsConstraints[i]; k++)
                    {

                        if((cuts_generated->Coefficients[p1+k]!=cuts_generated->Coefficients[p2+k])||(cuts_generated->Elements[p1+k]!=cuts_generated->Elements[p2+k]))
                        {
                            aux=0;

                            break;
                        }
                    }
                }
                if((aux==1)&&(j >= h_cut->numberConstrains))
                {
                    validated[j]=1;
                    minus_elements += cuts_generated->ElementsConstraints[j+1]-cuts_generated->ElementsConstraints[j];
                    cont_aux++;
                }
            }
        }
    }

//    printf("Number of repeat: %d \n", cont_aux);
    Cut_gpu* new_h_cut;
    new_h_cut = AllocationStructCut(cuts_generated->cont - minus_elements,cuts_generated->numberConstrains-cont_aux,cuts_generated->numberVariables);
    aux = 0;
    cont_aux = 0;
    new_h_cut->ElementsConstraints[0] = 0;
    for(i=0; i<cuts_generated->numberConstrains; i++)
    {
        if(validated[i]==0)
        {
            // if(i>=h_cut->numberConstrains)
            // {
            //printf("violation %f\n",violation[i-h_cut->numberConstrains]);
            // }
            new_h_cut->rightSide[aux] = cuts_generated->rightSide[i];
            new_h_cut->typeConstraints[aux] = cuts_generated->typeConstraints[i];
            for(j = cuts_generated->ElementsConstraints[i]; j < cuts_generated->ElementsConstraints[i+1]; j++)
            {
                new_h_cut->Coefficients[cont_aux] = cuts_generated->Coefficients[j];
                new_h_cut->Elements[cont_aux] = cuts_generated->Elements[j];
                cont_aux++;
            }
            new_h_cut->ElementsConstraints[aux + 1] = cont_aux;

            aux++;
        }

    }
    for(i=0; i<new_h_cut->numberVariables; i++)
    {
        new_h_cut->xAsterisc[i] = cuts_generated->xAsterisc[i];
    }
    //printf(" Durante: %d\n", new_h_cut->numberConstrains);
    free(validated);
    free(cuts_generated);
    free(Coefs_temp);
    free(Elements_temp);
    free(Pos_el_temp);
    free(rhs_temp);
    free(Coef1);
    free(v_aux);
//    free(violation);
    return new_h_cut;

}

Cut_gpu* createCutsOfPhaseTwo(Cut_gpu *h_cut, Cut_gpu_aux *cut_aux, solutionGpu *h_solution, int numberMaxConst, int nCuts, int precision, int nRuns)
{
    // double value =0 ;
    //double *value_violation = (double*)malloc(sizeof(double)*nCuts);
    int *Coef1 = (int*)malloc(sizeof(int)*(h_cut->numberVariables)); //free
    int *Coef2 = (int*)malloc(sizeof(int)*(h_cut->numberVariables));//free

    int *Coefs_temp = (int*)malloc(sizeof(int)*(nCuts*h_cut->numberVariables)); //free
    int *Elements_temp = (int*)malloc(sizeof(int)*(nCuts*h_cut->numberVariables)); //free
    int *Pos_el_temp = (int*)malloc(sizeof(int)*(nCuts+1) ); //free
    int *rhs_temp = (int*)malloc(sizeof(int)*(nCuts)); //free
    int temp = 0;
    Cut_gpu *cuts_generated; //free
//    printf("PASSA!");
    //  int nAll = 0;
    int i,k,j,n_1,d_1,n_2,d_2, qnt_1, el, rest_a, rhs1 = 0, rhs2 = 0,aux = 0, cont_aux=0;
    Pos_el_temp[0] = 0;
    for (i=0; i<nRuns; i++)
    {
        if(h_solution->SConst[0 + i*numberMaxConst]!=-1)
        {

            memset(Coef1,0,sizeof(int)*h_cut->numberVariables);
            memset(Coef2,0,sizeof(int)*h_cut->numberVariables);
            rhs1 = 0;
            rhs2 = 0;
            qnt_1 = h_solution->SPAux[i];
            n_1 = h_solution->SMult[0 + i*4 ];
            d_1 = h_solution->SMult[1 + i*4 ];
            n_2 = h_solution->SMult[2 + i*4 ];
            d_2 = h_solution->SMult[3 + i*4 ];

            for(rest_a = 0; rest_a <= qnt_1; rest_a++)
            {
                for(j = h_cut->ElementsConstraints[ h_solution->SConst[rest_a + i*numberMaxConst] ]; j< h_cut->ElementsConstraints[ h_solution->SConst[rest_a + i*numberMaxConst] + 1]; j++)
                {
                    el = h_cut->Elements[j];
                    Coef1[el] += h_cut->Coefficients[j] * n_1;

                    //value_tes += Coef[el] * d_cut->xAsterisc[el];
                }
                rhs1 += h_cut->rightSide[ h_solution->SConst[rest_a+ i*numberMaxConst] ] * n_1;
            }

//            for(j=0; j<h_cut->numberVariables; j++)
//            {
//                Coef1[j] = Coef1[j]<0 ? Coef1[j]/d_1 - 1 : Coef1[j]/d_1;
//            }
//            rhs1 = rhs1<0 ? rhs1/d_1-1 : rhs1/d_1;

            for(rest_a = qnt_1 + 1; rest_a < numberMaxConst; rest_a++)
            {
                for(j=h_cut->ElementsConstraints[ h_solution->SConst[rest_a+ i*numberMaxConst ] ]; j<h_cut->ElementsConstraints[  h_solution->SConst[rest_a + i*numberMaxConst] + 1]; j++)
                {
                    el = h_cut->Elements[j];
                    Coef2[el] += h_cut->Coefficients[j] * n_2;
                    //value_tes += Coef[el] * d_cut->xAsterisc[el];
                }
                rhs2 += h_cut->rightSide[ h_solution->SConst[rest_a+ i*numberMaxConst] ] * n_2;
            }
//            for(j=0; j<h_cut->numberVariables; j++)
//                Coef2[j] = Coef2[j]<0 ? Coef2[j]/d_2 - 1 : Coef2[j]/d_2;
//
//            rhs2 = rhs2<0 ? rhs2/d_2-1 : rhs2/d_2;
            //value_violation[cont_aux] = 0;
            for(j=0; j<h_cut->numberVariables; j++)
            {

                temp = (double)Coef1[j]/(double)d_1 + (double)Coef2[j]/(double)d_2;
                if(temp != 0)
                {
                    if(temp<0)
                    {
                        Coefs_temp[aux] = temp - 1;
                    }
                    else
                    {
                        Coefs_temp[aux] = temp;
                    }

                    Elements_temp[aux] = j;
                    //value_violation[cont_aux] += (Coefs_temp[aux])*h_cut->xAsterisc[j];
                    aux++;
                }
            }
            Pos_el_temp[cont_aux + 1] = aux;
            temp = (double)rhs1/(double)d_1 + (double)rhs2/(double)d_2;
            if(temp<0)
            {
                rhs_temp[cont_aux] = temp-1;
            }
            else
            {
                rhs_temp[cont_aux] = temp;
            }
            //value_violation[cont_aux] -= (rhs1+rhs2)*precision;
            cont_aux++;
        }
    }


    cuts_generated = AllocationStructCut(h_cut->cont + aux, h_cut->numberConstrains + nCuts, h_cut->numberVariables);
    for(i=0; i<cuts_generated->numberVariables; i++)
    {
        cuts_generated->xAsterisc[i] = h_cut->xAsterisc[i];
    }
    cuts_generated->ElementsConstraints[0] = 0;
    for(i=0; i<h_cut->numberConstrains; i++)
    {
        cuts_generated->rightSide[i] = h_cut->rightSide[i];
        cuts_generated->typeConstraints[i] = h_cut->typeConstraints[i];
        cuts_generated->ElementsConstraints[i+1] = h_cut->ElementsConstraints[i+1];
    }
    aux = 1;
    for(i = h_cut->numberConstrains; i<cuts_generated->numberConstrains; i++)
    {

        cuts_generated->rightSide[i] = rhs_temp[i - h_cut->numberConstrains];
        cuts_generated->typeConstraints[i] = LPC_CGGPUR2;
        cuts_generated->ElementsConstraints[i+1] = Pos_el_temp[aux] + h_cut->cont;
        aux++;//+ h_cut->cont;
    }
    for(i=0; i<h_cut->cont; i++)
    {
        cuts_generated->Coefficients[i] = h_cut->Coefficients[i];
        cuts_generated->Elements[i] = h_cut->Elements[i];
    }
    for(i= h_cut->cont; i<cuts_generated->cont; i++)
    {
        cuts_generated->Coefficients[i] = Coefs_temp[i - h_cut->cont];
        cuts_generated->Elements[i] = Elements_temp[i - h_cut->cont];
    }

    int *validated = (int*)malloc(sizeof(int)*cuts_generated->numberConstrains); //free
    memset(validated,0,sizeof(int)*cuts_generated->numberConstrains);
    aux = 1;
    cont_aux=0;
    int p1,p2, minus_elements = 0;
    for(i=0; i<cuts_generated->numberConstrains - 1; i++)
    {
        for(j=i+1; j<cuts_generated->numberConstrains; j++)
        {
            if(validated[j]==0)
            {
                if((cuts_generated->ElementsConstraints[i+1]-cuts_generated->ElementsConstraints[i])!=(cuts_generated->ElementsConstraints[j+1]-cuts_generated->ElementsConstraints[j])||(cuts_generated->rightSide[i]!=cuts_generated->rightSide[j]))
                {
                    aux=0;
                }
                else
                {
                    p1 = cuts_generated->ElementsConstraints[i];
                    p2 = cuts_generated->ElementsConstraints[j];
                    aux = 1;
                    for(k=0; k<cuts_generated->ElementsConstraints[i+1]-cuts_generated->ElementsConstraints[i]; k++)
                    {

                        if((cuts_generated->Coefficients[p1+k]!=cuts_generated->Coefficients[p2+k])||(cuts_generated->Elements[p1+k]!=cuts_generated->Elements[p2+k]))
                        {
                            aux=0;

                            break;
                        }
                    }
                }
                if((aux==1)&&(j >= h_cut->numberConstrains))
                {
                    validated[j]=1;
                    minus_elements += cuts_generated->ElementsConstraints[j+1]-cuts_generated->ElementsConstraints[j];
                    cont_aux++;
                }
            }
        }
    }
//    printf("Number of repeat: %d \n", cont_aux);
    Cut_gpu* new_h_cut;
    new_h_cut = AllocationStructCut(cuts_generated->cont - minus_elements,cuts_generated->numberConstrains-cont_aux,cuts_generated->numberVariables);
    aux = 0;
    cont_aux = 0;
    new_h_cut->ElementsConstraints[0] = 0;
    for(i=0; i<cuts_generated->numberConstrains; i++)
    {
        if(validated[i]==0)
        {
            new_h_cut->rightSide[aux] = cuts_generated->rightSide[i];
            new_h_cut->typeConstraints[aux] = cuts_generated->typeConstraints[i];
            // if(i>=h_cut->numberConstrains)
            // {
            //    printf("Violation: %f\n", value_violation[i-h_cut->numberConstrains]);
            // }
            for(j = cuts_generated->ElementsConstraints[i]; j < cuts_generated->ElementsConstraints[i+1]; j++)
            {
                new_h_cut->Coefficients[cont_aux] = cuts_generated->Coefficients[j];
                new_h_cut->Elements[cont_aux] = cuts_generated->Elements[j];
                cont_aux++;
            }

            new_h_cut->ElementsConstraints[aux + 1] = cont_aux;
            aux++;
        }

    }

    for(i=0; i<new_h_cut->numberVariables; i++)
    {
        new_h_cut->xAsterisc[i] = cuts_generated->xAsterisc[i];
    }

    //free(value_violation);
    free(validated);
    free(cuts_generated);
    free(Coefs_temp);
    free(Elements_temp);
    free(Pos_el_temp);
    free(rhs_temp);
    free(Coef1);
    free(Coef2);
    return new_h_cut;
}


Cut_gpu_aux* reallocCut(Cut_gpu *h_cut,Cut_gpu_aux *h_cut_aux)
{
    int i;
    Cut_gpu_aux *cut_aux;
    cut_aux = AllocationStructCutAux(h_cut->numberConstrains,h_cut->numberVariables,h_cut->cont);
    for(i=0; i<h_cut_aux->numberConstrains; i++)
    {
        strcpy(cut_aux->nameConstraints[i].name,h_cut_aux->nameConstraints[i].name);
        cut_aux->intervalMax[i] = h_cut_aux->intervalMax[i];
        cut_aux->intervalMin[i] = h_cut_aux->intervalMin[i];
    }
    for(i=0; i<h_cut->numberConstrains; i++)
    {
        if(h_cut->typeConstraints[i]==LPC_CGGPU)
        {
            //char buffer[20];
            //int aux = *cont;
            //sprintf(buffer,"CUTCGGPUR1(%d)",aux);
            //strcpy(cut_aux->nameConstraints[aux+h_cut_aux->numberConstrains].name,buffer);
            cut_aux->intervalMax[i+h_cut_aux->numberConstrains] = h_cut_aux->intervalMax[0];
            cut_aux->intervalMin[i+h_cut_aux->numberConstrains] = h_cut_aux->intervalMin[0];
            //aux++;
            //*cont = aux;
        }
    }
    //for(i=0; i<h_cut_aux->numberVariables; i++)
    //{
    //    strcpy(cut_aux->nameElements[i].name,h_cut_aux->nameElements[i].name);
    //}
    free(h_cut_aux);
    return cut_aux;
}

Cut_gpu_aux* reallocCutR2(Cut_gpu *h_cut,Cut_gpu_aux *h_cut_aux, int *cont)
{
    int i;
    Cut_gpu_aux *cut_aux;
    cut_aux = AllocationStructCutAux(h_cut->numberConstrains,h_cut->numberVariables,h_cut->cont);
    for(i=0; i<h_cut_aux->numberConstrains; i++)
    {
        strcpy(cut_aux->nameConstraints[i].name,h_cut_aux->nameConstraints[i].name);
        cut_aux->intervalMax[i] = h_cut_aux->intervalMax[i];
        cut_aux->intervalMin[i] = h_cut_aux->intervalMin[i];
    }
    for(i=0; i<h_cut->numberConstrains; i++)
    {

        if(h_cut->typeConstraints[i]==LPC_CGGPUR2)
        {
            char buffer[20];
            int aux = *cont;
            sprintf(buffer,"CUTCGGPUR2(%d)",aux);
            strcpy(cut_aux->nameConstraints[aux+h_cut_aux->numberConstrains].name,buffer);
            cut_aux->intervalMax[aux+h_cut_aux->numberConstrains] = h_cut_aux->intervalMax[0];
            cut_aux->intervalMin[aux+h_cut_aux->numberConstrains] = h_cut_aux->intervalMin[0];
            aux++;
            *cont = aux;
        }
    }
    for(i=0; i<h_cut_aux->numberVariables; i++)
    {
        strcpy(cut_aux->nameElements[i].name,h_cut_aux->nameElements[i].name);
    }
    free(h_cut_aux);
    return cut_aux;

}


void bubble_sort(int *vetor,int *pos, int n)
{
    int k, j, aux;
    for (k = 1; k < n; k++)
    {
        for (j = 0; j < n - 1; j++)
        {
            if (vetor[j] > vetor[j + 1])
            {
                aux          = vetor[j];
                vetor[j]     = vetor[j + 1];
                vetor[j + 1] = aux;
                aux = pos[j];
                pos[j] = pos[j+1];
                pos[j+1] = aux;
            }
        }
    }
}



int* returnOrdConstrainsNR(Cut_gpu *cut)
{
    int *nCofConst;
    nCofConst = (int*)malloc(sizeof(int)*(cut->numberConstrains*cut->numberConstrains));
    memset(nCofConst,0,sizeof(int)*(cut->numberConstrains*cut->numberConstrains));
    int i,j,e, el, cont;
    int nConst = cut->numberConstrains;
    for(i=0; i<nConst; i++)
    {
        for(j=0; j<nConst; j++)
        {
            cont = 0;
            if(i!=j)
            {
                for(e = cut->ElementsConstraints[i]; e < cut->ElementsConstraints[i+1]; e++)
                {
                    for(el =  cut->ElementsConstraints[j]; el <cut->ElementsConstraints[j+1]; el++)
                    {
                        if(cut->Elements[e]==cut->Elements[el])
                        {
                            cont++;
                            break;
                        }
                    }
                }
                nCofConst[i + j*nConst] = cont;
            }
        }
    }
    return nCofConst;
}

float* returnFolga(Cut_gpu *cut)
{

    float *folga;
    folga = (float*)malloc(sizeof(float)*cut->numberConstrains);

    memset(folga,0,sizeof(float)*(cut->numberConstrains));
    int r, el;
    float cont = 0;
    //int nConst = cut->numberConstrains;
    for( r = 0 ; r < cut->numberConstrains ; r++)
    {
        cont=0;
        for(el=cut->ElementsConstraints[r]; el<cut->ElementsConstraints[r+1]; el++)
            cont += cut->Coefficients[el] * (cut->xAsterisc[cut->Elements[el]]/1000);
        folga[r] = cut->rightSide[r] - cont;
        //printf("Folga da restrição %d: %f \n",r,folga[r]);
    }
    return folga;
}

void calcSetConstraint (int *setConstraint, int *pos_R1, int numberMaxConst,int numberConstrains, int *resR1, int *resNR1, int sizeR1, int sizeNR1, int *Similar, float *Folga, int nRuns, int szR )
{
    int i,j,k,aux,el;
    int pos_R1_aux = *pos_R1;
    int szNR = numberMaxConst - szR;
    int *p = (int*)malloc(sizeof(int)*sizeNR1); //free
    for(i = 0; i < nRuns; i++)
    {
        *pos_R1 = (pos_R1_aux + i)%sizeR1;
        memset(p,0,sizeof(int)*sizeNR1);
        for(j = 0; j<szR; j++)
        {
            setConstraint[i*numberMaxConst + j] = resR1[ (*pos_R1 + j)%sizeR1 ];
            for(k= 0; k<sizeNR1; k++)
            {
                p[k] += Similar[setConstraint[i*numberMaxConst + j] + resNR1[k]*numberConstrains];
                if(j==0)
                    p[k] -= 2*Folga[j];
            }
        }
        for (j = 1; j < sizeNR1; j++)
        {
            for (k = 0; k < sizeNR1 - 1; k++)
            {
                if (p[k] > p[k + 1])
                {
                    aux          = p[k];
                    p[k]     = p[k + 1];
                    p[k + 1] = aux;
                    aux          = resNR1[k];
                    resNR1[k]     = resNR1[k + 1];
                    resNR1[k + 1] = aux;
                }
            }
        }
        el = szR;
        for(j = 0; j<szNR; j++)
        {
            setConstraint[i*numberMaxConst + j + el] = resNR1[ j % sizeNR1];
        }
    }
    free(p);
    //*pos_R1 = pos_R1_aux;
}


void returnCutCG(Cut_gpu *ccg, CutCG* ccgr2, int precision)
{
    int i,j;
    int cont=0, aux;
    for(i = 0; i < ccg->numberConstrains; i++)
    {
        if((ccg->typeConstraints[i]==LPC_CGGPU)||(ccg->typeConstraints[i]==LPC_CGGPUR2))
        {
            cont++;
        }
    }

    VecInt **cutElem;
    VecDbl **cutCoef;
    VecDbl *cutrhs;
    VecInt *cutnelem;
    VecDbl *cutviolation;
    VecInt *cutsense;
    VecInt *cutdominated;
    VecStr *cutname;
    double value_violation = 0;
    //int *Coef1 = (int*)malloc(sizeof(int)*(h_cut->numberVariables));
    //int *Coef2 = (int*)malloc(sizeof(int)*(h_cut->numberVariables));
    int c;
    int nAll = 0;
    ALLOCATE_VECTOR(cutElem,VecInt*,cont);
    ALLOCATE_VECTOR(cutCoef,VecDbl*,cont);

    for ( c=0 ; (c<cont) ; ++c )
    {
        cutElem[c] = VInt_create();
        cutCoef[c] = VDbl_create();
        nAll++;
    }

    cutrhs = VDbl_create();
    cutnelem = VInt_create();
    cutviolation = VDbl_create();
    cutsense = VInt_create();
    cutdominated = VInt_create();
    cutname = VStr_create(256);

    c = 0;
    int el;
    for(i = 0; i<ccg->numberConstrains; i++)
    {
        value_violation = 0;
        if((ccg->typeConstraints[i]==LPC_CGGPU)||(ccg->typeConstraints[i]==LPC_CGGPUR2))
        {
            aux = 0;
            for(j = ccg->ElementsConstraints[i]; j < ccg->ElementsConstraints[i+1]; j++)
            {
                VDbl_pushBack(cutCoef[c], ccg->Coefficients[j]);
                el = ccg->Elements[j];
                value_violation += ccg->Coefficients[j]*ccg->xAsterisc[el];
                VInt_pushBack(cutElem[c], ccg->Elements[j]);
                aux++;
            }
            value_violation -= ccg->rightSide[i]*precision;
            VDbl_pushBack(cutrhs,ccg->rightSide[i]);
            VInt_pushBack(cutsense,0);
            VInt_pushBack(cutdominated,0);
            VInt_pushBack(cutnelem,aux);
            VDbl_pushBack(cutviolation,value_violation);
            c++;
        }
    }
    VInt_free(&ccgr2->cutnelem);
    VDbl_free(&ccgr2->cutrhs);
    VInt_free(&ccgr2->cutsense);
    VDbl_free(&ccgr2->cutviolation);
    VInt_free(&ccgr2->cutdominated);
    VStr_free(&ccgr2->cutname);

    for ( c=0 ; (c<ccgr2->nAllocations) ; ++c )
    {
        VInt_free(&ccgr2->cutElem[c]);
        VDbl_free(&ccgr2->cutCoef[c]);
    }

    //VInt_free(&ccg->cutElem[0]);
    //VDbl_free(&ccg->cutCoef[0]);
    free(ccgr2->cutElem);
    free(ccgr2->cutCoef);

    ccgr2->cutCoef = cutCoef;
    ccgr2->cutElem = cutElem;
    ccgr2->nAllocations = nAll;
    ccgr2->cutnelem = cutnelem;
    ccgr2->cutrhs = cutrhs;
    ccgr2->cutviolation = cutviolation;
    ccgr2->cutsense = cutsense;
    ccgr2->cutdominated = cutdominated;
    ccgr2->cutname = cutname;
}


void  shuffle_Set(int *vec, int nSetConstrains, int n)
{
    struct timeval time;
    gettimeofday(&time, NULL);
    srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
    int i, j, aux ;
    int *num_temp = (int*)malloc(sizeof(int)*nSetConstrains); //free
    int *vec_aux = (int*)malloc(sizeof(int)*nSetConstrains); //free
    aux  =  n/nSetConstrains;
    for(i = 0; i < aux ; i++)
    {

        for(j = 0 ; j<nSetConstrains; j++)
        {

            num_temp[j] = rand()%RAND_MAX;
            vec_aux[j] = vec[i*nSetConstrains + j];
        }
        bubble_sort(num_temp,vec_aux,nSetConstrains);
        for(j = 0 ; j<nSetConstrains; j++)
        {
            vec[i*nSetConstrains + j] = vec_aux[j];
        }
    }
    free(num_temp);
    free(vec_aux);
}

int fat(int n)
{
    if(n <= 1)
    {
        return 1;
    }
    else
    {
        return n*fat(n-1);
    }
}

int combinacao_Phase1(int n, int p)
{
    int res = fat(n) / (fat(p)*fat(n-p));
    return res;
}

Cut_gpu* complementCutPhase1(Cut_gpu *h_cut, int k)
{
    int i,j,w,z, el, cons = 0, aux=0;
    int sz_const=0, qnt_add = 0;
    int *constraints = (int*)malloc(sizeof(int)*h_cut->numberConstrains);
    for(i=0; i<h_cut->numberConstrains; i++)
    {
        if(h_cut->typeConstraints[i]==RES_RR)
        {
            constraints[sz_const] = i;
            sz_const++;
        }
    }
    for (i=2; i<=k; i++)
    {
        qnt_add += ( sz_const  - (i - 1) );
    }
    int *coef_news = (int*)malloc(sizeof(int)*(h_cut->numberVariables*qnt_add));
    int *rhs_news = (int*)malloc(sizeof(int)*(qnt_add));
    memset(coef_news,0,sizeof(int)*(h_cut->numberVariables*qnt_add));
    memset(rhs_news,0,sizeof(int)*(qnt_add));

    for(w=2; w<=k; w++)
    {
        for(i=0; i<sz_const - (w-1); i++)
        {
            for(j = i; j < w; j++)
            {
                cons = constraints[j];
                for(z = h_cut->ElementsConstraints[cons]; z< h_cut->ElementsConstraints[ cons + 1]; z++)
                {
                    el = h_cut->Elements[z];
                    coef_news[el + aux*h_cut->numberVariables] += h_cut->Coefficients[z];
                    //value_tes += Coef[el] * d_cut->xAsterisc[el];
                }
                rhs_news[aux]+= h_cut->rightSide[cons];
            }
            aux++;
        }
    }
    aux = 0;
    for(i = 0; i< h_cut->numberVariables*qnt_add; i++)
    {
        if(coef_news[i]!=0)
        {
            aux++;
        }

    }
    Cut_gpu* new_cuts = AllocationStructCut(h_cut->cont+aux,h_cut->numberConstrains + qnt_add,h_cut->numberVariables);

    for(i=0; i<h_cut->numberConstrains; i++)
    {
        new_cuts->rightSide[i] = h_cut->rightSide[i];
        new_cuts->ElementsConstraints[i] = h_cut->ElementsConstraints[i];
        new_cuts->typeConstraints[i] = h_cut->typeConstraints[i];
    }
    new_cuts->ElementsConstraints[i] = h_cut->ElementsConstraints[i];
    for(i=0; i<h_cut->numberVariables; i++)
    {
        new_cuts->xAsterisc[i] = h_cut->xAsterisc[i];
    }
    for(i=0; i<h_cut->cont; i++)
    {
        new_cuts->Coefficients[i] = h_cut->Coefficients[i];
        new_cuts->Elements[i] = h_cut->Elements[i];
    }
    //New constraints
    for(i = h_cut->numberConstrains; i<new_cuts->numberConstrains; i++)
    {
        new_cuts->rightSide[i] = rhs_news[ i - h_cut->numberConstrains ];
        new_cuts->typeConstraints[i] = RES_RR;
    }
    aux = h_cut->cont;
    w  = h_cut->numberConstrains + 1;
    for(i=0; i<qnt_add; i++)
    {
        for(j= 0; j < h_cut->numberVariables; j++)
        {
            if(coef_news[j + i*h_cut->numberVariables] !=0)
            {
                new_cuts->Coefficients[aux] = coef_news[j + i*h_cut->numberVariables];
                new_cuts->Elements[aux] = j;
                aux++;
            }
        }
        new_cuts->ElementsConstraints[w] = aux;
        w++;
    }
    free(h_cut);
    free(constraints);
    free(coef_news);
    free(rhs_news);
    return new_cuts;
}





Cut_gpu* complementCutPhase1_RR(Cut_gpu *h_cut, Cut_gpu_aux *h_cut_aux, int nRR_cggpu )
{
    int i, j;
    int  fin, cont = 0;// ini,
    int a, b,numberC = 0;
    //printf("numero de recurso: %d \n",nRR_cggpu);
    // ini = h_cut_aux->intervalMin[0];
    fin = h_cut_aux->intervalMax[0];
    for(i = 0; i < h_cut_aux->numberConstrains; i++)
    {
        if(h_cut->typeConstraints[i]==RES_RR)
        {
            numberC++;
        }
    }

    int *res_for_time = (int*)malloc(sizeof(int)*((fin+1)*nRR_cggpu));
    int *sz_r = (int*)malloc(sizeof(int)*(fin+1));

    for(i = 0 ; i < ((fin + 1)*nRR_cggpu); i++)
    {
        res_for_time[i] = -1;
    }

    memset(sz_r,0, sizeof(int)*(fin + 1));

    for(i = 0; i < h_cut_aux->numberConstrains; i++)
    {
        if(h_cut->typeConstraints[i]==RES_RR)
        {
            sscanf(h_cut_aux->nameConstraints[i].name,"resR(%d,%d)",&a,&b);
            int s = sz_r[b];
            res_for_time[s + b*nRR_cggpu] = i;
            sz_r[b]++;
            //printf("%d , %d - %s %d\n", a,b,h_cut_aux->nameConstraints[i].name, i );
        }
        cont++;
    }

    int qnt_X_add = 0;

    for(i = 0; i<(fin+1); i++)
    {
        //printf("Instante de tempo: %d\n",i);
        //    for(j=0;j<sz_r[i];j++){
        //        if(j>0){
        //qnt_X_add += combinacao_Phase1(sz_r[i],4);
        if(sz_r[i]>=3)
        {
            a = combinacao_Phase1(sz_r[i],3);
            qnt_X_add += combinacao_Phase1(sz_r[i],3);
        }
        if(sz_r[i]>=2)
        {
            a = combinacao_Phase1(sz_r[i],2);
            qnt_X_add += combinacao_Phase1(sz_r[i],2);
        }
        //        }
        //printf("%d %d \t", res_for_time[j + i*nRR_cggpu], qnt_X_add);
        //    }
        //printf("\n");
    }

    int *coef_news = (int*)malloc(sizeof(int)*(h_cut->numberVariables*qnt_X_add));
    int *rhs_news = (int*)malloc(sizeof(int)*(qnt_X_add));
    memset(coef_news,0,sizeof(int)*(h_cut->numberVariables*qnt_X_add));
    memset(rhs_news,0,sizeof(int)*(qnt_X_add));

    int z,w,y;
    int c, aux = 0;

    for(i =0 ; i<(fin+1); i++)
    {
        if(sz_r[i]>=3)
        {
            for(j = 0; j< sz_r[i]-2; j++)
            {
                a =  res_for_time[j + i*nRR_cggpu];
                for(w = j+1; w < sz_r[i]-1; w++)
                {
                    b = res_for_time[w + i*nRR_cggpu];
                    for(z = w+1; z<sz_r[i]; z++)
                    {
                        c = res_for_time[z + i*nRR_cggpu];


                        for(y = h_cut->ElementsConstraints[a]; y< h_cut->ElementsConstraints[ a + 1]; y++)
                        {
                            int el = h_cut->Elements[y];
                            coef_news[el + aux*h_cut->numberVariables] += h_cut->Coefficients[y];
                            //value_tes += Coef[el] * d_cut->xAsterisc[el];
                        }
                        for(y = h_cut->ElementsConstraints[b]; y< h_cut->ElementsConstraints[ b + 1]; y++)
                        {
                            int el = h_cut->Elements[y];
                            coef_news[el + aux*h_cut->numberVariables] += h_cut->Coefficients[y];
                            //value_tes += Coef[el] * d_cut->xAsterisc[el];
                        }
                        for(y = h_cut->ElementsConstraints[c]; y< h_cut->ElementsConstraints[ c+ 1]; y++)
                        {
                            int el = h_cut->Elements[y];
                            coef_news[el + aux*h_cut->numberVariables] += h_cut->Coefficients[y];
                            //value_tes += Coef[el] * d_cut->xAsterisc[el];
                        }
                        rhs_news[aux]+= h_cut->rightSide[a] + h_cut->rightSide[b] + h_cut->rightSide[c];


                        //printf("%d %d %d\n", a,b,c);

                        aux++;
                    }
                }
            }
        }
        if(sz_r[i]>=2)
        {
            for(j = 0; j< sz_r[i]-1; j++)
            {
                a =  res_for_time[j + i*nRR_cggpu];
                for(w = j+1; w < sz_r[i]; w++)
                {
                    b = res_for_time[w + i*nRR_cggpu];
                    //printf("a: %d b: %d\n", a,b);
                    for(y = h_cut->ElementsConstraints[a]; y< h_cut->ElementsConstraints[ a + 1]; y++)
                    {
                        int el = h_cut->Elements[y];
                        coef_news[el + aux*h_cut->numberVariables] += h_cut->Coefficients[y];
                        //value_tes += Coef[el] * d_cut->xAsterisc[el];
                    }
                    for(y = h_cut->ElementsConstraints[b]; y< h_cut->ElementsConstraints[ b + 1]; y++)
                    {
                        int el = h_cut->Elements[y];
                        coef_news[el + aux*h_cut->numberVariables] += h_cut->Coefficients[y];
                        //value_tes += Coef[el] * d_cut->xAsterisc[el];
                    }
                    rhs_news[aux]+= h_cut->rightSide[a] + h_cut->rightSide[b];
                    aux++;
                }
            }
        }


    }
    //printf("aux: %d, cont: %d\n", aux, qnt_X_add);
    //printf("Intervalo : %d a %d\n",h_cut_aux->intervalMin[0], h_cut_aux->intervalMax[0]);
    aux = 0;
    for(i = 0; i< h_cut->numberVariables*qnt_X_add; i++)
    {
        if(coef_news[i]!=0)
        {
            aux++;
        }

    }
    Cut_gpu* new_cuts = AllocationStructCut(h_cut->cont+aux,h_cut->numberConstrains + qnt_X_add,h_cut->numberVariables);

    for(i=0; i<h_cut->numberConstrains; i++)
    {
        new_cuts->rightSide[i] = h_cut->rightSide[i];
        new_cuts->ElementsConstraints[i] = h_cut->ElementsConstraints[i];
        new_cuts->typeConstraints[i] = h_cut->typeConstraints[i];
    }
    new_cuts->ElementsConstraints[i] = h_cut->ElementsConstraints[i];
    for(i=0; i<h_cut->numberVariables; i++)
    {
        new_cuts->xAsterisc[i] = h_cut->xAsterisc[i];
    }
    for(i=0; i<h_cut->cont; i++)
    {
        new_cuts->Coefficients[i] = h_cut->Coefficients[i];
        new_cuts->Elements[i] = h_cut->Elements[i];
    }
    //New constraints
    for(i = h_cut->numberConstrains; i<new_cuts->numberConstrains; i++)
    {
        new_cuts->rightSide[i] = rhs_news[ i - h_cut->numberConstrains ];
        new_cuts->typeConstraints[i] = RES_RR;
    }
    aux = h_cut->cont;
    w  = h_cut->numberConstrains + 1;
    for(i=0; i<qnt_X_add; i++)
    {
        for(j= 0; j < h_cut->numberVariables; j++)
        {
            if(coef_news[j + i*h_cut->numberVariables] !=0)
            {
                new_cuts->Coefficients[aux] = coef_news[j + i*h_cut->numberVariables];
                new_cuts->Elements[aux] = j;
                aux++;
            }
        }
        new_cuts->ElementsConstraints[w] = aux;
        w++;
    }
    free(h_cut);

//    for(i = ini; i<= fin ; i++)
//    {
//        for(j = 0; j< nRR_cggpu; j++)
//        {
//
//
//        }
//    }

    free(coef_news);
    free(rhs_news);
    free(sz_r);
    free(res_for_time);
    return new_cuts;

}
