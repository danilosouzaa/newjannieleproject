/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * mip_compact: routines for the creation of compact MIP formulations
 *              for the whole problem or for just one project
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */
#include "cut_energetic.h"

#define VERBOSE  1

#define MIN_VIOL 0.02


CutE *CutE_create(const Instance *inst)
{

    CutE *cen;
    ALLOCATE_INI(cen,CutE);

    cen->cutrhs = VDbl_create();
    cen->cutviolation = VDbl_create();
    cen->cutnelem = VInt_create();
    cen->cutsense = VInt_create();
    cen->cutname = VStr_create(256);
    cen->cutdominated = VInt_create();
    cen->nAllocations = 0;

    cen->nMaxWindowWithCutEnergeticViolated = 0;
    cen->nMaxWindowWithCutEnergeticViolatedPerc = 0.0;
    cen->nS = 0;
    cen->nE = 0;


    cen->nMinWindowWithCutEnergeticViolated = INT_MAX;
    cen->nMinWindowWithCutEnergeticViolatedPerc = 0.0;
    cen->nmS = 0;
    cen->nmE = 0;

    cen->inst = inst;

    return cen;
}


void CutE_free( CutE**_cutE )
{

    CutE *cutE = *_cutE;

    for(int nr = 0 ; nr < cutE->nAllocations ; nr++)
    {
        VInt_free(&cutE->cutElem[nr]);
        VDbl_free(&cutE->cutCoef[nr]);
    }

    free(cutE->cutElem);
    free(cutE->cutCoef);

    VStr_free(&cutE->cutname);
    VDbl_free(&cutE->cutrhs);
    VInt_free(&cutE->cutnelem);
    VInt_free(&cutE->cutsense);
    VInt_free(&cutE->cutdominated);
    VDbl_free(&cutE->cutviolation);

    free( cutE );
    *_cutE = NULL;

}

void CutE_add_cuts_energetic_parallel( CutE *cen,  LinearProgram *lp, const int *origCols,  LinearProgram *origLP, const Instance *inst, double timeLeft, int continuous,  int **maxTJM,  double maxcuts, int nround )
{

    int nCols = lp_cols(lp);


    assert(lp);
    assert(inst);
    assert(origLP);


    const double *xf = lp_x(lp);
    char name[256];
    int i;

    /* maximum number of indexes in a var name */
#define MAX_IDX 16
    char prefix[STR_SIZE];
    int idx[MAX_IDX];
    int maxT = 0;

    const int nJobs = Inst_nJobs(inst);
    const int nModes = Inst_nMaxModes(inst);

    // first and last time job has some allocation
    int **firstTJ;
    int **lastTJ;
    ALLOCATE_VECTOR(firstTJ,int*,nJobs)
    ALLOCATE_VECTOR(lastTJ,int*,nJobs)
    for ( int j=0 ; (j<nJobs) ; ++j )
    {
        ALLOCATE_VECTOR_INI(firstTJ[j],int,nModes)
        ALLOCATE_VECTOR_INI(lastTJ[j],int,nModes)
        for(int m = 0 ; m < nModes ; m++)
            firstTJ[j][m] = INT_MAX;

    }



    // checking maximum time where there is some allocation
    for ( i=0 ; (i<lp_cols(lp)) ; ++i )
    {

        lp_col_name( origLP, origCols[i], name );
        parseName( name, prefix, idx );

        if (prefix[0]=='x')
        {
            //printf("%s %s", name); fflush(stdout);
            const int j = idx[0];
            const int m = idx[1];
            const int t = idx[2];
            maxT = MAX( maxT, t );
            firstTJ[j][m]= MIN(firstTJ[j][m], t);
            lastTJ[j][m] = MAX(lastTJ[j][m], t);
        }
    }

    int nTimes = maxT+1;
    // allocation of x
    double ***x;
    ALLOCATE_VECTOR( x, double **, nJobs );
    ALLOCATE_VECTOR( x[0], double *, nJobs*nModes );
    ALLOCATE_VECTOR( x[0][0], double, nJobs*nModes*nTimes );
    for ( int j=1 ; (j<nJobs) ; ++j )
        x[j] = x[j-1] + nModes;

    for ( int j=0; (j<nJobs) ; ++j )
        for ( int m=0 ; (m<nModes) ; ++m )
            x[j][m] = x[0][0] + j*nModes*nTimes + m*nTimes;

    // filling default values
    {
        double *it, *itEnd;
        it = x[0][0];
        itEnd = it + (nJobs*nModes*nTimes);
        for ( ; (it<itEnd) ; it++ )
            *it = 0.0;
    }


    // for all variables
    // filling x from fractional solution
    for ( i=0 ; (i<lp_cols(lp)) ; ++i )
    {

        lp_col_name( origLP, origCols[i], name );
        if (tolower(name[0]=='x'))
        {
            parseName( name, prefix, idx );
            if (prefix[0]=='x')
                x[idx[0]][idx[1]][idx[2]] = xf[i];
        }
    }

    /* x indexes pre-processed problem*/
    int ***xIdx;
    ALLOCATE_VECTOR( xIdx, int **, nJobs );
    ALLOCATE_VECTOR( xIdx[0], int *, nJobs*nModes );
    ALLOCATE_VECTOR( xIdx[0][0], int, nJobs*nModes*nTimes );
    for ( int j=1 ; (j<nJobs) ; ++j )
        xIdx[j] = xIdx[j-1] + nModes;

    for ( int j=0; (j<nJobs) ; ++j )
        for ( int m=0 ; (m<nModes) ; ++m )
            xIdx[j][m] = xIdx[0][0] + j*nModes*nTimes + m*nTimes;
    // filling default values
    {
        int *it, *itEnd;
        it = xIdx[0][0];
        itEnd = it + (nJobs*nModes*nTimes);
        for ( ; (it<itEnd) ; it++ )
            *it = -1;
    }

    /* filling indexes from variables in this pre-processed problem */
    for ( i=0 ; (i<lp_cols(lp)) ; ++i )
    {

        lp_col_name( origLP, origCols[i], name );
        if (tolower(name[0])=='x')
        {
            parseName( name, prefix, idx );
            if (prefix[0]=='x')
                xIdx[idx[0]][idx[1]][idx[2]] = i;
        }
    }



    int ****unitsES;
    ALLOCATE_VECTOR(unitsES,int***,nJobs);
    for(int j = 0 ; j < nJobs ; j++)
    {
        const Job* job = Inst_job(inst,j);
        ALLOCATE_VECTOR(unitsES[j],int**,nModes);
        for(int m = 0 ; m < Job_nModes(job) ; m++)
        {
            const Mode *mode = Job_mode(job,m);
            ALLOCATE_VECTOR(unitsES[j][m],int*,nTimes);
            for(int s = 0 ; s < nTimes ; s++)
            {
                ALLOCATE_VECTOR_INI(unitsES[j][m][s],int,nTimes);
                for(int e = s+1 ; e < nTimes ; e++)
                {
                    int part1 = e-s+1;
                    int part2 = MAX(0,Mode_duration(mode)-MAX(0,s-firstTJ[j][m]));
                    int part3 = MAX(0,Mode_duration(mode)-MAX(0,lastTJ[j][m]+Mode_duration(mode)-e));
                    unitsES[j][m][s][e] = MIN(part1,MIN(part2,part3));
                }

            }
        }
    }

    ALLOCATE_VECTOR(cen->cutElem,VecInt*,Inst_nResR(inst)*nTimes*nTimes);
    ALLOCATE_VECTOR(cen->cutCoef,VecDbl*,Inst_nResR(inst)*nTimes*nTimes);


    int nCut2 =0;

    for(int r = 0 ; r < Inst_nResR(inst) ; r++)
    {
        int e =0;
        for(int s = 0 ; (s < nTimes) && e<(nTimes-1); s++)
        {
            e=s+10;
            IntDblPair *cutIdxCo;
            ALLOCATE_VECTOR(cutIdxCo, IntDblPair, nCols);

            int size = 0;
            double slack =0.0;
            for ( int t=0 ; t< nTimes; ++t )
            {
                for ( int j= 0 ; j< nJobs ; ++j )
                {
                    const Job *job = Inst_job( inst, j );
                    const int nModes = Job_nModes( job );
                    for ( int m=0 ; (m<nModes) ; ++m )
                    {

                        const Mode *mode = Job_mode( job, m );
                        if(Mode_duration(mode)==0) continue;
                        if (!Mode_isFeasible(inst,mode))
                            continue;
                        if(xIdx[j][m][t] == -1 || fabs(x[j][m][t]) == 1.0 || xf[xIdx[j][m][t]] <= 0.0000 || unitsES[j][m][s][e]==0) continue;
                        int rm = Mode_idxResROnMode(inst,mode,r);
                        if(rm==-1) continue;
                        int cjr = Mode_useResR(mode,rm);
                        if(t<=e && t+Mode_duration(mode) >= s  )
                        {

                            /*  int value=0;
                              if( t>=s) {
                                  if(t+Mode_duration(mode) <= e){
                                     value =  Mode_duration(mode);
                                  }else{
                                      value =  Mode_duration(mode) - (t+Mode_duration(mode)-1-e);
                                  }
                              }else{
                                  if(t+Mode_duration(mode) <= e){
                                      value = Mode_duration(mode) -(s-t);
                                  }else{
                                      value = Mode_duration(mode) - (t+Mode_duration(mode)-1-e) - (s-t);
                                  }

                              }*/


                            int idelem =  xIdx[j][m][t];
                            //  printf(" idelem %d x(%d,%d,%d), unitsResourceES[%d][%d][%d][%d] %d value %d duration %d, cjr %d", idelem, j,m,t,j,m,s,e, unitsResourceES[j][m][s][e], value, Mode_duration(mode), cjr );
                            cutIdxCo[size].a = idelem;
                            cutIdxCo[size].b =  (unitsES[j][m][s][e]*cjr);
                            slack += xf[cutIdxCo[size].a] * cutIdxCo[size].b;
                            size++;
                        }
                    }
                }
            }

            slack -= ((e-s+1)*Inst_capResR(inst,r));


            if(slack > 0.0002)
            {
//                  printf("\n r%d s %d e %d :  (e-s)*Inst_capResR(inst,r) %d\n", r, s, e, (e-s+1)*Inst_capResR(inst,r));// getchar();

                /* adding cuts*/

                cen->cutElem[nCut2] = VInt_create();
                cen->cutCoef[nCut2] = VDbl_create();
                cen->nAllocations++;

                char nameCut[STR_SIZE];
                sprintf( nameCut, "cutE(%d,%d,%d)#%d#%d", r, s, e, nround, nCut2 );
            //    printf( "energ cutE(%d,%d,%d)#%d#%d", r, s, e, nround, nCut2 );
                VStr_pushBack(cen->cutname,nameCut);

                for(int sel  = 0 ; sel <size; sel++)
                {
                    int elem = cutIdxCo[sel].a;
                    VInt_pushBack(cen->cutElem[nCut2],elem);
                    VDbl_pushBack(cen->cutCoef[nCut2],cutIdxCo[sel].b);
           //         printf( " %f * %f: %d \n", cutIdxCo[sel].b,xf[elem],elem);
                }
                int *cutIdx = VInt_getPtr(cen->cutElem[nCut2]);
                double *cutCoef = VDbl_getPtr(cen->cutCoef[nCut2]);
                VDbl_pushBack(cen->cutrhs,((e-s+1)*Inst_capResR(inst,r)));
                VInt_pushBack(cen->cutsense,0);
                VInt_pushBack(cen->cutdominated,0);
                VInt_pushBack(cen->cutnelem,size);
                VDbl_pushBack(cen->cutviolation,slack);
                CutP_quick_sort_vec(cutIdx, cutCoef, size);
  //              printf( " rhs%d, sense %d, nelem %d, viol %f \n",((e-s+1)*Inst_capResR(inst,r)), 0,size,slack);
//                getchar();
                nCut2++;
            }

            free(cutIdxCo);
        }
    }


    for(int j = 0 ; j < nJobs ; j++)
    {
        const Job* job = Inst_job(inst,j);
        for(int m = 0 ; m < Job_nModes(job) ; m++)
        {
            for(int s = 0 ; s < nTimes ; s++)
            {
                free(unitsES[j][m][s]);
            }
            free(unitsES[j][m]);
        }
        free(unitsES[j]);
        free(firstTJ[j]);
        free(lastTJ[j]);
    }
    free(unitsES);
    free(firstTJ);
    free(lastTJ);

    free(xIdx[0][0]);
    free(xIdx[0]);
    free(xIdx);
    free(x[0][0]);
    free(x[0]);
    free(x);

    /*  for(int r = 0 ; r < Inst_nResR(inst) ; r++){
          for(int s = 0 ; (s < nTimes); s++){

              for(int e = s ; e<(nTimes) ; e++){
              //for(int e = s ; e < (s+10) && e<(nTimes-1) ; e++){
                  //e=s+10;
                  int units =0;
                  for(int j = 0 ; j < nJobs ; j++){
                      const Job* job = Inst_job(inst,j);
                      for(int m = 0 ; m < Job_nModes(job) ; m++){
                          const Mode *mode = Job_mode(job,m);
                          int rm = Mode_idxResROnMode(inst,mode,r);
                          if(rm==-1) continue;
                          int cjr = Mode_useResR(mode,rm);
                          units += (Inst_unitsResourceES(inst,j,m,s,e)*cjr);
                      }
                  }
                //  printf("s %d e %d : units %d (e-s)*Inst_capResR(inst,r) %d\n", s, e, units, (e-s)*Inst_capResR(inst,r));// getchar();
                   if(units>(e-s)*Inst_capResR(inst,r)){
                      if(units >0){
                          printf("s %d e %d : units %d (e-s)*Inst_capResR(inst,r) %d\n", s, e, units, (e-s)*Inst_capResR(inst,r));// getchar();
                      }
                  }

              }
          }

      }*/



}



void CutE_add_cuts_energetic_model_parallel_orig( CutE *cen,  LinearProgram *lp, const int *origCols,  LinearProgram *origLP, const Instance *inst, double timeLeft, int continuous,  int **maxTJM,  double maxcuts, int nround )
{

    int nCols = lp_cols(lp);


    assert(lp);
    assert(inst);
    assert(origLP);

    const double *xf = lp_x(lp);
    char name[256];
    int i;

    /* maximum number of indexes in a var name */
#define MAX_IDX 16
    char prefix[STR_SIZE];
    int idx[MAX_IDX];
    int auxidZ[MAX_IDX];
    int maxT = 0;

    const int nJobs = Inst_nJobs(inst);
    const int nModes = Inst_nMaxModes(inst);

    // first and last time job has some allocation for each mode
    int **firstTJ;
    int **lastTJ;
    ALLOCATE_VECTOR(firstTJ,int*,nJobs)
    ALLOCATE_VECTOR(lastTJ,int*,nJobs)
    for ( int j=0 ; (j<nJobs) ; ++j )
    {
        ALLOCATE_VECTOR_INI(firstTJ[j],int,nModes)
        ALLOCATE_VECTOR_INI(lastTJ[j],int,nModes)
        for(int m = 0 ; m < nModes ; m++)
            firstTJ[j][m] = INT_MAX;

    }



    // checking maximum time where there is some allocation
    for ( i=0 ; (i<lp_cols(lp)) ; ++i )
    {

        lp_col_name( origLP, origCols[i], name );
        parseName( name, prefix, idx );

        if (prefix[0]=='x')
        {
            //printf("%s %s", name); fflush(stdout);
            const int j = idx[0];
            const int m = idx[1];
            const int t = idx[2];
            maxT = MAX( maxT, t );
            firstTJ[j][m]= MIN(firstTJ[j][m], t);
            lastTJ[j][m] = MAX(lastTJ[j][m], t);
        }
    }

    int nTimes = maxT+1;
    // allocation of x
    double ***x;
    ALLOCATE_VECTOR( x, double **, nJobs );
    ALLOCATE_VECTOR( x[0], double *, nJobs*nModes );
    ALLOCATE_VECTOR( x[0][0], double, nJobs*nModes*nTimes );
    for ( int j=1 ; (j<nJobs) ; ++j )
        x[j] = x[j-1] + nModes;

    for ( int j=0; (j<nJobs) ; ++j )
        for ( int m=0 ; (m<nModes) ; ++m )
            x[j][m] = x[0][0] + j*nModes*nTimes + m*nTimes;

    // filling default values
    {
        double *it, *itEnd;
        it = x[0][0];
        itEnd = it + (nJobs*nModes*nTimes);
        for ( ; (it<itEnd) ; it++ )
            *it = 0.0;
    }


    // for all variables
    // filling x from fractional solution
    for ( i=0 ; (i<lp_cols(lp)) ; ++i )
    {

        lp_col_name( origLP, origCols[i], name );
        if (tolower(name[0]=='x'))
        {
            parseName( name, prefix, idx );
            if (prefix[0]=='x')
                x[idx[0]][idx[1]][idx[2]] = xf[i];
        }
    }

    /* x indexes pre-processed problem*/
    int ***xIdx;
    ALLOCATE_VECTOR( xIdx, int **, nJobs );
    ALLOCATE_VECTOR( xIdx[0], int *, nJobs*nModes );
    ALLOCATE_VECTOR( xIdx[0][0], int, nJobs*nModes*nTimes );
    for ( int j=1 ; (j<nJobs) ; ++j )
        xIdx[j] = xIdx[j-1] + nModes;

    for ( int j=0; (j<nJobs) ; ++j )
        for ( int m=0 ; (m<nModes) ; ++m )
            xIdx[j][m] = xIdx[0][0] + j*nModes*nTimes + m*nTimes;
    // filling default values
    {
        int *it, *itEnd;
        it = xIdx[0][0];
        itEnd = it + (nJobs*nModes*nTimes);
        for ( ; (it<itEnd) ; it++ )
            *it = -1;
    }

    /* filling indexes from variables in this pre-processed problem */
    for ( i=0 ; (i<lp_cols(lp)) ; ++i )
    {

        lp_col_name( origLP, origCols[i], name );
        if (tolower(name[0])=='x')
        {
            parseName( name, prefix, idx );
            if (prefix[0]=='x')
                xIdx[idx[0]][idx[1]][idx[2]] = i;
        }
    }



    int ****unitsES;
    ALLOCATE_VECTOR(unitsES,int***,nJobs);
    for(int j = 0 ; j < nJobs ; j++)
    {
        const Job* job = Inst_job(inst,j);
        ALLOCATE_VECTOR(unitsES[j],int**,nModes);
        for(int m = 0 ; m < Job_nModes(job) ; m++)
        {
            const Mode *mode = Job_mode(job,m);
            ALLOCATE_VECTOR(unitsES[j][m],int*,nTimes);
            for(int s = 0 ; s < nTimes ; s++)
            {
                ALLOCATE_VECTOR_INI(unitsES[j][m][s],int,nTimes);
                for(int e = s+1 ; e < nTimes ; e++)
                {
                    int part1 = e-s+1;
                    int part2 = MAX(0,Mode_duration(mode)-MAX(0,s-firstTJ[j][m]));
                    int part3 = MAX(0,Mode_duration(mode)-MAX(0,lastTJ[j][m]+Mode_duration(mode)-e));
                    unitsES[j][m][s][e] = MIN(part1,MIN(part2,part3));
                }

            }
        }
    }


    ALLOCATE_VECTOR(cen->cutElem,VecInt*,Inst_nResR(inst)*nTimes*nTimes);
    ALLOCATE_VECTOR(cen->cutCoef,VecDbl*,Inst_nResR(inst)*nTimes*nTimes);

    VecInt *vecIdx = VInt_create();
    VecDbl *vecCo = VDbl_create();
    VecInt *vecNCut = VInt_create();
    VecDbl *vecValue = VDbl_create();
    VecStr *vecNames = VStr_create(256);
    int nCut2 =0;

    LinearProgram *mipCutm;

    for(int r = 0 ; r < Inst_nResR(inst) ; r++)
    {
        for(int s = 0 ; (s < nTimes) ; s++)
        {
            for(int e = s+1 ; e<(nTimes-1); e++)
            {
                printf("R %d S %d E %d \n", r,s,e);
                mipCutm = lp_create();
                lp_set_print_messages( mipCutm, 1);

                VecStr *namesbin = VStr_create(STR_SIZE);
                VecInt *lIdx = VInt_create();
                VecDbl *objbin = VDbl_create();
                int nColB =0;

                IntDblPair *cutIdxCo;
                ALLOCATE_VECTOR(cutIdxCo, IntDblPair, nCols);
                int idx[nCols];
                double coef[nCols];

                int idxlr[nCols][2];
                double coeflr[nCols][2];

                int size = 0;
                double slack =0.0;
                int nrow = 0;
                for ( int t=0 ; t< nTimes; ++t )
                {
                    for ( int j= 0 ; j< nJobs ; ++j )
                    {
                        const Job *job = Inst_job( inst, j );
                        const int nModes = Job_nModes( job );
                        for ( int m=0 ; (m<nModes) ; ++m )
                        {
                            if(xIdx[j][m][t] == -1 || xf[ xIdx[j][m][t]] <= 0.0 || unitsES[j][m][s][e] <= 0) continue;
                           // if(xIdx[j][m][t] == -1 || fabs(x[j][m][t]) == 1.0 || xf[ xIdx[j][m][t]] <= 0.0 || unitsES[j][m][s][e]<=0) continue;

                            const Mode *mode = Job_mode( job, m );
                            if(Mode_duration(mode)==0) continue;
                            if (!Mode_isFeasible(inst,mode)) continue;

                            int rm = Mode_idxResROnMode(inst,mode,r);
                            if(rm==-1) continue;

                            int cjr = Mode_useResR(mode,rm);
                            if(cjr<=0 ) continue;

                            if(t<=e && t+Mode_duration(mode) >= s  )
                            {

                                int value=0;
                                if( t>=s)
                                {
                                    if(t+Mode_duration(mode) <= e)
                                    {
                                        value =  Mode_duration(mode);
                                    }
                                    else
                                    {
                                        value =  Mode_duration(mode) - (t+Mode_duration(mode)-1-e);
                                    }
                                }
                                else
                                {
                                    if(t+Mode_duration(mode) <= e)
                                    {
                                        value = Mode_duration(mode) -(s-t);
                                    }
                                    else
                                    {
                                        value = (e-s+1);
                                    }
                                }
                                if(value<=0) continue;
                                printf("  xIdx[j][m][t] %d, xf %f,  idelem %d x(%d,%d,%d), unitES[%d][%d][%d][%d] %d value %d duration %d, cjr %d\n",  xIdx[j][m][t], xf[ xIdx[j][m][t]], nColB, j,m,t,j,m,s,e, unitsES[j][m][s][e], value, Mode_duration(mode), cjr );

                                char vnameZ[STR_SIZE];
                                sprintf( vnameZ, "l(%d,%d,%d)", j, m,t);
                                VStr_pushBack( namesbin, vnameZ );
                                VDbl_pushBack(objbin, (value*cjr*xf[ xIdx[j][m][t]])*-1);
                                VInt_pushBack( lIdx, nColB );
                                idx[size] = nColB;
                                coef[size] = (value*cjr*xf[ xIdx[j][m][t]]);
                                idxlr[nrow][0] = nColB;
                                coeflr[nrow][0] = 1.0;
                                nColB++;
                                size++;


                                sprintf( vnameZ, "r(%d,%d,%d)", j, m,t);
                                VStr_pushBack( namesbin, vnameZ );
                                VDbl_pushBack(objbin,  (unitsES[j][m][s][e]*cjr)*-1);
                                VInt_pushBack( lIdx, nColB );
                                idx[size] = nColB;
                                coef[size] = (unitsES[j][m][s][e]*cjr);
                                idxlr[nrow][1] = nColB;
                                coeflr[nrow][1] = 1.0;
                                nColB++;
                                size++;
                                nrow++;
                            }
                        }
                    }
                }
                //   printf("\n\n");

                if(size>0 && nColB>0)
                {
                    lp_add_bin_cols( mipCutm, nColB, VDbl_getPtr(objbin), VStr_ptr(namesbin) );

                    VecStr *namesint = VStr_create(STR_SIZE);
                    char vnamee[256];
                    sprintf( vnamee, "a");
                    VStr_pushBack( namesint, vnamee );

                    double *objint;
                    ALLOCATE_VECTOR_INI(objint,double,1);

                    char *integer;
                    ALLOCATE_VECTOR(integer,char,1);
                    integer[0] = True;

                    double *lb;
                    ALLOCATE_VECTOR_INI(lb,double,1);
                    double *up;
                    ALLOCATE_VECTOR_INI(up,double,1);
                    FILL(lb,0,1,1);
                    FILL(up,0,1,1);
                    objint[0] = -((e-s+1)*Inst_capResR(inst,r))*-1;

                    lp_add_cols(mipCutm,1,objint,lb,up,integer,VStr_ptr( (VecStr *)namesint));

                    free(integer);
                    free(lb);
                    free(up);
                    free(objint);
                    VStr_free(&namesint);

                    char namewz[STR_SIZE];
                    sprintf( namewz, "cutE(%d,%d,%d)",r,s,e);
                    //  printf("cons(3).(%d,%d)#(%d,%d)",t,r,idj[iy],jmw);
                    lp_add_row( mipCutm, size,idx,coef, namewz, 'G', ((e-s+1)*Inst_capResR(inst,r))+EPS);

                    for(int nr = 0 ; nr < nrow ; nr++)
                    {
                        char name[STR_SIZE];
                        sprintf( name, "oneside(%d)", nr);
                        lp_add_row( mipCutm,2,idxlr[nr],coeflr[nr],name,'E', 1.0);
                    }

                  //  printf("writing model"); fflush(stdout);
                  //  lp_write_lp(mipCutm,"energ.lp");
                    //lp_set_concurrentMIP(mipCutm,1);
                    //lp_set_method(mipCutm,4);
                    //lp_set_seed(mipCutm,100000);
                    int status = lp_optimize(mipCutm);
                   // getchar();

                    if(status == LP_OPTIMAL)
                    {
                      //  lp_write_sol(mipCutm,"energ.sol");
                        slack = lp_obj_value(mipCutm)*-1;
                      //  getchar();
                        if(slack > 0.0002)
                        {

                            int c = 0;

                            const double *xf2 = lp_x(mipCutm);
                            int nc = lp_cols(mipCutm);
                            for ( i=0 ; (i<nc) ; ++i )
                            {
                                char nameC[STR_SIZE];
                                lp_col_name( mipCutm, i, nameC );
                                if (xf2[i]>EPS)
                                {
                                    parseName( nameC, prefix, auxidZ );
                                    int j = auxidZ[0];
                                    int m = auxidZ[1];
                                    int t = auxidZ[2];

                                    if (tolower(nameC[0]=='l') || tolower(nameC[0]=='r'))
                                    {
                                        cutIdxCo[c].a = xIdx[j][m][t];//cZidx[c];
                                        cutIdxCo[c].b = coef[i];
                           //             printf(" %f*%d (%f): %s + ", cutIdxCo[c].b,cutIdxCo[c].a, xf2[i] ,nameC);
                                        c++;
                                    }
                                }
                            }



                            /* adding cuts*/

                            cen->cutElem[nCut2] = VInt_create();
                            cen->cutCoef[nCut2] = VDbl_create();
                            cen->nAllocations++;

                            char nameCut[STR_SIZE];
                            sprintf( nameCut, "cutE(%d,%d,%d)#%d#%d", r, s, e, nround, nCut2 );
                         //   printf( "\n cutE(%d,%d,%d)#%d#%d", r, s, e, nround, nCut2 );
                            VStr_pushBack(cen->cutname,nameCut);

                            for(int sel  = 0 ; sel <c; sel++)
                            {
                                int elem = cutIdxCo[sel].a;
                                VInt_pushBack(cen->cutElem[nCut2],elem);
                                VDbl_pushBack(cen->cutCoef[nCut2],cutIdxCo[sel].b/xf[ cutIdxCo[sel].a]);
                       //         printf( " %f * %f: %d \n", cutIdxCo[sel].b,xf[elem],elem);
                            }
                            int *cutIdx = VInt_getPtr(cen->cutElem[nCut2]);
                            double *cutCoef = VDbl_getPtr(cen->cutCoef[nCut2]);
                            VDbl_pushBack(cen->cutrhs,((e-s+1)*Inst_capResR(inst,r)));
                            VInt_pushBack(cen->cutsense,0);
                            VInt_pushBack(cen->cutdominated,0);
                            VInt_pushBack(cen->cutnelem,c);
                            VDbl_pushBack(cen->cutviolation,slack);
                            CutP_quick_sort_vec(cutIdx, cutCoef, c);
                            if((e-s+1) > cen->nMaxWindowWithCutEnergeticViolated )
                            {
                                cen->nMaxWindowWithCutEnergeticViolated = (e-s+1);
                                cen->nMaxWindowWithCutEnergeticViolatedPerc = (double) (e-s+1)/(double) nTimes;
                                cen->nS = s;
                                cen->nE = e;
                            }
                            if((e-s+1) < cen->nMinWindowWithCutEnergeticViolated )
                            {
                                cen->nMinWindowWithCutEnergeticViolated = (e-s+1);
                                cen->nMinWindowWithCutEnergeticViolatedPerc = (double) (e-s+1)/(double) nTimes;
                                cen->nmS = s;
                                cen->nmE = e;
                            }
  //                          printf( " rhs %d, sense %d, nelem %d, viol %f\n",((e-s+1)*Inst_capResR(inst,r)), 0,c,slack);
//                            getchar();
                            nCut2++;
                        }
                    }
                    free(cutIdxCo);
                    lp_free(&mipCutm);
                    VStr_free(&namesbin);
                    VInt_free(&lIdx);
                    VDbl_free(&objbin);
                }
            }
        }
    }

    VInt_free(&vecIdx);
    VInt_free(&vecNCut);
    VDbl_free(&vecCo);
    VDbl_free(&vecValue);
    VStr_free(&vecNames);

    for(int j = 0 ; j < nJobs ; j++)
    {
        const Job* job = Inst_job(inst,j);
        for(int m = 0 ; m < Job_nModes(job) ; m++)
        {
            for(int s = 0 ; s < nTimes ; s++)
            {
                free(unitsES[j][m][s]);
            }
            free(unitsES[j][m]);
        }
        free(unitsES[j]);
        free(firstTJ[j]);
        free(lastTJ[j]);
    }
    free(unitsES);
    free(firstTJ);
    free(lastTJ);

    free(xIdx[0][0]);
    free(xIdx[0]);
    free(xIdx);
    free(x[0][0]);
    free(x[0]);
    free(x);
}


void CutE_add_cuts_energetic_model_parallel( CutE *cen,  LinearProgram *lp, const int *origCols,  LinearProgram *origLP, const Instance *inst, double timeLeft, int continuous,  int **maxTJM,  double maxcuts, int nround )
{

//    int nCut =  lp_rows(lp);
    int nCols = lp_cols(lp);

    assert(lp);
    assert(inst);
    assert(origLP);

    const double *xf = lp_x(lp);
    char name[256];
    int i;

    /* maximum number of indexes in a var name */
#define MAX_IDX 16
    char prefix[STR_SIZE];
    int idx[MAX_IDX];
    int auxidZ[MAX_IDX];
    int maxT = 0;

    const int nJobs = Inst_nJobs(inst);
    const int nModes = Inst_nMaxModes(inst);

    // first and last time job has some allocation
    int *firstTJ;
    int *lastTJ;
    ALLOCATE_VECTOR(firstTJ,int,nJobs)
    ALLOCATE_VECTOR(lastTJ,int,nJobs)
    int **firstTJM;
    int **lastTJM;
    ALLOCATE_VECTOR(firstTJM,int*,nJobs)
    ALLOCATE_VECTOR(lastTJM,int*,nJobs)
    for ( int j=0 ; (j<nJobs) ; ++j )
    {
        firstTJ[j]= INT_MAX;
        lastTJ[j]= 0;
        ALLOCATE_VECTOR_INI(firstTJM[j],int,nModes)
        ALLOCATE_VECTOR_INI(lastTJM[j],int,nModes)
        for(int m = 0 ; m < nModes ; m++)
            firstTJM[j][m] = INT_MAX;
    }


    // checking maximum time where there is some allocation
    for ( i=0 ; (i<lp_cols(lp)) ; ++i )
    {
        lp_col_name( origLP, origCols[i], name );
        parseName( name, prefix, idx );

        if (prefix[0]=='x')
        {
            const int j = idx[0];
            const int m = idx[1];
            const int t = idx[2];
            maxT = MAX( maxT, t );
            firstTJM[j][m]= MIN(firstTJM[j][m], t);
            lastTJM[j][m] = MAX(lastTJM[j][m], t);
            firstTJ[j]= MIN(firstTJ[j], t);
            lastTJ[j] = MAX(lastTJ[j], t);
        }
    }

    int nTimes = maxT+1;
    // allocation of x
    double ***x;
    ALLOCATE_VECTOR( x, double **, nJobs );
    ALLOCATE_VECTOR( x[0], double *, nJobs*nModes );
    ALLOCATE_VECTOR( x[0][0], double, nJobs*nModes*nTimes );
    for ( int j=1 ; (j<nJobs) ; ++j )
        x[j] = x[j-1] + nModes;

    for ( int j=0; (j<nJobs) ; ++j )
        for ( int m=0 ; (m<nModes) ; ++m )
            x[j][m] = x[0][0] + j*nModes*nTimes + m*nTimes;

    // filling default values
    {
        double *it, *itEnd;
        it = x[0][0];
        itEnd = it + (nJobs*nModes*nTimes);
        for ( ; (it<itEnd) ; it++ )
            *it = 0.0;
    }


    // for all variables
    // filling x from fractional solution
    for ( i=0 ; (i<lp_cols(lp)) ; ++i )
    {
        lp_col_name( origLP, origCols[i], name );
        if (tolower(name[0]=='x'))
        {
            parseName( name, prefix, idx );
            if (prefix[0]=='x')
                x[idx[0]][idx[1]][idx[2]] = xf[i];
        }
    }

    /* x indexes pre-processed problem*/
    int ***xIdx;
    ALLOCATE_VECTOR( xIdx, int **, nJobs );
    ALLOCATE_VECTOR( xIdx[0], int *, nJobs*nModes );
    ALLOCATE_VECTOR( xIdx[0][0], int, nJobs*nModes*nTimes );
    for ( int j=1 ; (j<nJobs) ; ++j )
        xIdx[j] = xIdx[j-1] + nModes;

    for ( int j=0; (j<nJobs) ; ++j )
        for ( int m=0 ; (m<nModes) ; ++m )
            xIdx[j][m] = xIdx[0][0] + j*nModes*nTimes + m*nTimes;
    // filling default values
    {
        int *it, *itEnd;
        it = xIdx[0][0];
        itEnd = it + (nJobs*nModes*nTimes);
        for ( ; (it<itEnd) ; it++ )
            *it = -1;
    }

    /* filling indexes from variables in this pre-processed problem */
    for ( i=0 ; (i<lp_cols(lp)) ; ++i )
    {

        lp_col_name( origLP, origCols[i], name );
        if (tolower(name[0])=='x')
        {
            parseName( name, prefix, idx );
            if (prefix[0]=='x')
                xIdx[idx[0]][idx[1]][idx[2]] = i;
        }
    }



    int ****unitsES;
    ALLOCATE_VECTOR(unitsES,int***,nJobs);
    for(int j = 0 ; j < nJobs ; j++)
    {
        const Job* job = Inst_job(inst,j);
        ALLOCATE_VECTOR(unitsES[j],int**,nModes);
        for(int m = 0 ; m < Job_nModes(job) ; m++)
        {
            const Mode *mode = Job_mode(job,m);
            ALLOCATE_VECTOR(unitsES[j][m],int*,nTimes);
            for(int s = 0 ; s < nTimes ; s++)
            {
                ALLOCATE_VECTOR_INI(unitsES[j][m][s],int,nTimes);
                for(int e = s+1 ; e < nTimes ; e++)
                {
                    int part1 = e-s+1;
                    int part2 = MAX(0,Mode_duration(mode)-MAX(0,s-firstTJM[j][m]));
                    int part3 = MAX(0,Mode_duration(mode)-MAX(0,lastTJM[j][m]+Mode_duration(mode)-e));
                    unitsES[j][m][s][e] = MIN(part1,MIN(part2,part3));
                   // if(j==1&&m==0)
                    //    printf("part1 %d, part2 %d, part3 %d: unitsES[j%d][m%d][s%d][e%d] %d \n",part1,part2,part3, j,m,s,e,unitsES[j][m][s][e]);
                }
            }
        }
    }



    int ****unitsESJ;
    ALLOCATE_VECTOR(unitsESJ,int***,nJobs);
    int ****minMode;
    ALLOCATE_VECTOR(minMode,int***,nJobs);

    for(int j = 0 ; j < nJobs ; j++)
    {
        const Job* job = Inst_job(inst,j);
        ALLOCATE_VECTOR(unitsESJ[j],int**,nTimes);
        ALLOCATE_VECTOR(minMode[j],int**,nTimes);
        for(int s = 0 ; s < nTimes ; s++)
        {
            ALLOCATE_VECTOR_INI(unitsESJ[j][s],int*,nTimes);
            ALLOCATE_VECTOR(minMode[j][s],int*,nTimes);
            for(int e = s+1 ; e < nTimes ; e++)
            {
                 ALLOCATE_VECTOR_INI(unitsESJ[j][s][e],int,Inst_nResR(inst) );
                 ALLOCATE_VECTOR_INI(minMode[j][s][e],int,Inst_nResR(inst));
                 for(int r = 0 ; r < Inst_nResR(inst) ; r++)
                 {

                    int minValue=INT_MAX;
                    for(int m = 0 ; m < Job_nModes(job) ; m++)
                    {
                        const Mode *mode = Job_mode(job,m);
                        int rm = Mode_idxResROnMode(inst,mode,r);
                        int usage = 0;
                        if(rm != -1) usage =  Mode_useResR(mode,rm);

                        if(minValue > unitsES[j][m][s][e] * usage){
                            unitsESJ[j][s][e][r] = unitsES[j][m][s][e];
                            minValue = unitsES[j][m][s][e] * usage;
                            minMode[j][s][e][r] = usage;
                        }

                    }
                }
            }
        }
    }


    ALLOCATE_VECTOR(cen->cutElem,VecInt*,Inst_nResR(inst)*nTimes*nTimes);
    ALLOCATE_VECTOR(cen->cutCoef,VecDbl*,Inst_nResR(inst)*nTimes*nTimes);

    VecInt *vecIdx = VInt_create();
    VecDbl *vecCo = VDbl_create();
    VecInt *vecNCut = VInt_create();
    VecDbl *vecValue = VDbl_create();
    VecStr *vecNames = VStr_create(256);
    int nCut2 =0;

    LinearProgram *mipCutm;

    for(int r = 0 ; r < Inst_nResR(inst) ; r++)
    {
        for(int s = 0 ; (s < nTimes) ; s++)
        {
            for(int e = s+1 ; e<(nTimes); e++)
            {

                mipCutm = lp_create();
                lp_set_print_messages( mipCutm, 1);

                VecStr *namesbin = VStr_create(STR_SIZE);
                VecInt *lIdx = VInt_create();
                VecDbl *objbin = VDbl_create();
                int nColB =0;

                IntDblPair *cutIdxCo;
                ALLOCATE_VECTOR(cutIdxCo, IntDblPair, nCols);

                int idx2[nCols];
                double coef[nCols];

                int idxlr[nJobs*nModes*nTimes][2];
                double coeflr[nJobs*nModes*nTimes][2];

                int size = 0;
                double slack =0.0;
                int contEq = 0;
                for ( int j=0 ; j< nJobs ; j++ )
                {
                    const Job *job = Inst_job( inst, j );
                    const int nModes = Job_nModes( job );

                    int nColBJ = -1;
                    if(unitsESJ[j][s][e][r] > 0.0 ){
                        char vnameZ[STR_SIZE];
                        sprintf( vnameZ, "r(%d)", j);
                        VStr_pushBack( namesbin, vnameZ );
                        double va =  (double )(unitsESJ[j][s][e][r]*minMode[j][s][e][r]);
                        VDbl_pushBack(objbin,va*-1) ;
                        VInt_pushBack( lIdx, nColB );
                        idx2[size] = nColB;
                        coef[size] = va;
                      //  printf(" r(%d) coef duration %d * usage%d \n",j, unitsESJ[j][s][e][r],minMode[j][s][e][r]);
                        nColBJ = nColB;
                        nColB++;
                        size++;
                    }

                    for ( int m=0 ; (m<nModes) ; ++m )
                    {

                        const Mode *mode = Job_mode( job, m );
                        if(Mode_duration(mode)==0) continue;
                        if(!Mode_isFeasible(inst,mode)) continue;
                        int rm = Mode_idxResROnMode(inst,mode,r);
                        int cjr = 0;
                        if(rm!=-1) cjr = Mode_useResR(mode,rm);;
                        if(cjr<=0) continue;
                        for ( int t=0 ; t< nTimes; ++t )
                        {

                            if(xIdx[j][m][t] == -1 || xf[ xIdx[j][m][t]] <= 0.0 ) continue;
                         //   if(xIdx[j][m][t] == -1 ) continue;

                            if(t<=e && t+Mode_duration(mode) >= s  )
                            {

                              int t2 = t+Mode_duration(mode)-1;
                               int value = MAX(0,MIN(t2,e)-MAX(t,s)+1);

/*                                int value=0;
                                if( t>=s)
                                {
                                    if(t+Mode_duration(mode) <= e)
                                    {
                                        value =  Mode_duration(mode);
                                    }
                                    else
                                    {
                                        value =  Mode_duration(mode) - (t+Mode_duration(mode)-1-e);
                                    }
                                }
                                else
                                {
                                    if(t+Mode_duration(mode) <= e)
                                    {
                                        value = Mode_duration(mode) -(s-t);
                                    }
                                    else
                                    {
                                        value = (e-s+1);
                                    }
                                }*/
                                if(value<=0) continue;

                                char vnameZ[STR_SIZE];
                                sprintf( vnameZ, "l(%d,%d,%d)", j, m,t);
                                VStr_pushBack( namesbin, vnameZ );
                                VDbl_pushBack(objbin, (value*cjr*xf[ xIdx[j][m][t]])*-1);
                                VInt_pushBack( lIdx, nColB );
                                idx2[size] = nColB;
                                coef[size] = (value*cjr*xf[ xIdx[j][m][t]]);
                                idxlr[contEq][0] = nColB;
                                coeflr[contEq][0] = 1.0;
                                idxlr[contEq][1] = nColBJ;
                                coeflr[contEq][1] = 1.0;
                                contEq++;
                                nColB++;
                            //    printf("  xIdx[j][m][t] %d, xf %f,  idelem %d x(%d,%d,%d), value %d duration %d, cjr %d coef[size] %f\n",  xIdx[j][m][t], xf[ xIdx[j][m][t]], nColB, j,m,t,value, Mode_duration(mode), cjr, coef[size]);
                                size++;
                            }
                        }
                    }
                }


                if(size>0 && nColB>0)
                {
                    lp_add_bin_cols( mipCutm, nColB, VDbl_getPtr(objbin), VStr_ptr(namesbin) );

                    VecStr *namesint = VStr_create(STR_SIZE);
                    char vnamee[256];
                    sprintf( vnamee, "a");
                    VStr_pushBack( namesint, vnamee );

                    double *objint;
                    ALLOCATE_VECTOR_INI(objint,double,1);

                    char *integer;
                    ALLOCATE_VECTOR(integer,char,1);
                    integer[0] = True;

                    double *lb;
                    ALLOCATE_VECTOR_INI(lb,double,1);
                    double *up;
                    ALLOCATE_VECTOR_INI(up,double,1);
                    FILL(lb,0,1,1);
                    FILL(up,0,1,1);
                    objint[0] = -((e-s+1)*Inst_capResR(inst,r))*-1;

                    double rhs = ((e-s+1)*Inst_capResR(inst,r));

                    lp_add_cols(mipCutm,1,objint,lb,up,integer,VStr_ptr( (VecStr *)namesint));

                    free(integer);
                    free(lb);
                    free(up);
                    free(objint);
                    VStr_free(&namesint);

                    char namewz[STR_SIZE];
                    sprintf( namewz, "cutE(%d,%d,%d)",r,s,e);
                    lp_add_row( mipCutm, size,idx2,coef, namewz, 'G', ((e-s+1)*Inst_capResR(inst,r))+EPS);

                    for(int e = 0 ; e< contEq ; e++)
                    {
                        char name[STR_SIZE];
                        sprintf( name, "oneside(%d)", e);
                        if(idxlr[e][1]!=-1)
                            lp_add_row( mipCutm,2,idxlr[e],coeflr[e],name,'L', 1);
                    }

                   // printf("writing model"); fflush(stdout);
                   // lp_write_lp(mipCutm,"energ.lp");
//lp_set_concurrentMIP(mipCutm,1);
                //lp_set_method(mipCutm,4);

                   //lp_set_seed(mipCutm,100000);
                    int status = lp_optimize(mipCutm);
                  //  getchar();

                    if(status == LP_OPTIMAL)
                    {
                    //    lp_write_sol(mipCutm,"energ.sol");
                        slack = lp_obj_value(mipCutm)*-1;
                    //    getchar();
                        if(slack > 0.0002)
                        {


                            int c = 0;
                            const double *xf2 = lp_x(mipCutm);
                            int nc = lp_cols(mipCutm);
                            for ( i=0 ; (i<nc) ; ++i )
                            {
                                char nameC[STR_SIZE];
                                lp_col_name( mipCutm, i, nameC );
                                if (xf2[i]>EPS)
                                {
                                   parseName( nameC, prefix, auxidZ );
                                   if (tolower(nameC[0]=='l') ) //|| tolower(nameC[0]=='r'))
                                   {
                                        int j = auxidZ[0];
                                        int m = auxidZ[1];
                                        int t = auxidZ[2];

                                        cutIdxCo[c].a = xIdx[j][m][t];
                                        cutIdxCo[c].b = coef[i];// /xf[xIdx[j][m][t]];
                                      //  printf("l %f*%d (%f): %s + ", cutIdxCo[c].b,cutIdxCo[c].a, xf2[i] ,nameC);
                                        c++;
                                    }

                                    if (tolower(nameC[0]=='r') ) //|| tolower(nameC[0]=='r'))
                                    {
                                        rhs -= coef[i];
                                    }
                                    }
                                }

                            /* adding cuts*/
                            cen->cutElem[nCut2] = VInt_create();
                            cen->cutCoef[nCut2] = VDbl_create();
                            cen->nAllocations++;

                            char nameCut[STR_SIZE];
                            sprintf( nameCut, "cutE(%d,%d,%d)#%d#%d", r, s, e, nround, nCut2 );
                         //   printf( "\n cutE(%d,%d,%d)#%d#%d", r, s, e, nround, nCut2 );
                            VStr_pushBack(cen->cutname,nameCut);

                            for(int sel  = 0 ; sel <c; sel++)
                            {
                                int elem = cutIdxCo[sel].a;
                                VInt_pushBack(cen->cutElem[nCut2],elem);
                                VDbl_pushBack(cen->cutCoef[nCut2],cutIdxCo[sel].b);//xf[ cutIdxCo[sel].a]);
                          //      printf( " %f * %f: %d \n", cutIdxCo[sel].b,xf[elem],elem);
                            }
                            int *cutIdx = VInt_getPtr(cen->cutElem[nCut2]);
                            double *cutCoef = VDbl_getPtr(cen->cutCoef[nCut2]);
                            VDbl_pushBack(cen->cutrhs,rhs);//((e-s+1)*Inst_capResR(inst,r)));
                            VInt_pushBack(cen->cutsense,0);
                            VInt_pushBack(cen->cutdominated,0);
                            VInt_pushBack(cen->cutnelem,c);
                            VDbl_pushBack(cen->cutviolation,slack);
                            CutP_quick_sort_vec(cutIdx, cutCoef, c);
                            if((e-s+1) > cen->nMaxWindowWithCutEnergeticViolated )
                            {
                                cen->nMaxWindowWithCutEnergeticViolated = (e-s+1);
                                cen->nMaxWindowWithCutEnergeticViolatedPerc = (double) (e-s+1)/(double) nTimes;
                                cen->nS = s;
                                cen->nE = e;
                            }
                            if((e-s+1) < cen->nMinWindowWithCutEnergeticViolated )
                            {
                                cen->nMinWindowWithCutEnergeticViolated = (e-s+1);
                                cen->nMinWindowWithCutEnergeticViolatedPerc = (double) (e-s+1)/(double) nTimes;
                                cen->nmS = s;
                                cen->nmE = e;
                            }
                         //   printf( " rhs%d, sense %d, nelem %d, viol %f\n",((e-s+1)*Inst_capResR(inst,r)), 0,c,slack);
                            printf( " rhs%f, sense %d, nelem %d, viol %f\n",rhs, 0,c,slack);
                            //getchar();
                            nCut2++;
                        }
                    }
                    free(cutIdxCo);
                    lp_free(&mipCutm);
                    VStr_free(&namesbin);
                    VInt_free(&lIdx);
                    VDbl_free(&objbin);
                }
            }
        }
    }

    VInt_free(&vecIdx);
    VInt_free(&vecNCut);
    VDbl_free(&vecCo);
    VDbl_free(&vecValue);
    VStr_free(&vecNames);


        for(int ncdA = 0; ncdA < VDbl_size(cen->cutrhs) ; ncdA++) {
            for(int ncdB = ncdA+1; ncdB < VDbl_size(cen->cutrhs) ; ncdB++) {
                if(VInt_get(cen->cutdominated,ncdA) == 1) break;
                if(VInt_get(cen->cutdominated,ncdB) == 1) continue;
                CutP_dominanceBetweenTwo( lp, ncdA, VInt_getPtr(cen->cutElem[ncdA]), VDbl_getPtr(cen->cutCoef[ncdA]), VDbl_get(cen->cutrhs,ncdA), VInt_get(cen->cutsense,ncdA), VInt_get(cen->cutnelem,ncdA), ncdB, VInt_getPtr(cen->cutElem[ncdB]), VDbl_getPtr(cen->cutCoef[ncdB]), VDbl_get(cen->cutrhs,ncdB), VInt_get(cen->cutsense,ncdB), VInt_get(cen->cutnelem,ncdB),  cen->cutdominated);
                //getchar();
            }
        }



    for(int j = 0 ; j < nJobs ; j++)
    {
        const Job* job = Inst_job(inst,j);
        for(int m = 0 ; m < Job_nModes(job) ; m++)
        {
            for(int s = 0 ; s < nTimes ; s++)
            {
                free(unitsES[j][m][s]);
            }
            free(unitsES[j][m]);
        }
        for(int s = 0 ; s < nTimes ; s++)
        {
            for(int e = s+1 ; e < nTimes ; e++)
            {
                free(unitsESJ[j][s][e]);
                free(minMode[j][s][e]);
            }
            free(unitsESJ[j][s]);
            free(minMode[j][s]);
        }
        free(unitsESJ[j]);
        free(minMode[j]);
        free(unitsES[j]);
    }
    free(unitsESJ);
    free(minMode);
    free(unitsES);
    free(firstTJ);
    free(lastTJ);

    free(xIdx[0][0]);
    free(xIdx[0]);
    free(xIdx);
    free(x[0][0]);
    free(x[0]);
    free(x);
}



void CutE_add_cuts_energetic_guloso_parallel( CutE *cen,  LinearProgram *lp, const int *origCols,  LinearProgram *origLP, const Instance *inst, double timeLeft, int continuous,  int **maxTJM,  double maxcuts, int nround)
{

    int nCols = lp_cols(lp);


    assert(lp);
    assert(inst);
    assert(origLP);

    const double *xf = lp_x(lp);
    char name[256];
    int i;

    /* maximum number of indexes in a var name */
#define MAX_IDX 16
    char prefix[STR_SIZE];
    int idx[MAX_IDX];
    int maxT = 0;

    const int nJobs = Inst_nJobs(inst);
    const int nModes = Inst_nMaxModes(inst);

    // first and last time job has some allocation
    int *firstTJ;
    int *lastTJ;
    ALLOCATE_VECTOR(firstTJ,int,nJobs)
    ALLOCATE_VECTOR(lastTJ,int,nJobs)
    int **firstTJM;
    int **lastTJM;
    ALLOCATE_VECTOR(firstTJM,int*,nJobs)
    ALLOCATE_VECTOR(lastTJM,int*,nJobs)
    for ( int j=0 ; (j<nJobs) ; ++j )
    {
        firstTJ[j]= INT_MAX;
        lastTJ[j]= 0;
        ALLOCATE_VECTOR_INI(firstTJM[j],int,nModes)
        ALLOCATE_VECTOR_INI(lastTJM[j],int,nModes)
        for(int m = 0 ; m < nModes ; m++)
            firstTJM[j][m] = INT_MAX;
    }


    // checking maximum time where there is some allocation
    for ( i=0 ; (i<lp_cols(lp)) ; ++i )
    {
        lp_col_name( origLP, origCols[i], name );
        parseName( name, prefix, idx );

        if (prefix[0]=='x')
        {
            //printf("%s %s", name); fflush(stdout);
            const int j = idx[0];
            const int m = idx[1];
            const int t = idx[2];
            maxT = MAX( maxT, t );
            firstTJM[j][m]= MIN(firstTJM[j][m], t);
            lastTJM[j][m] = MAX(lastTJM[j][m], t);
            firstTJ[j]= MIN(firstTJ[j], t);
            lastTJ[j] = MAX(lastTJ[j], t);
        }
    }

    int nTimes = maxT+1;
    // allocation of x
    double ***x;
    ALLOCATE_VECTOR( x, double **, nJobs );
    ALLOCATE_VECTOR( x[0], double *, nJobs*nModes );
    ALLOCATE_VECTOR( x[0][0], double, nJobs*nModes*nTimes );
    for ( int j=1 ; (j<nJobs) ; ++j )
        x[j] = x[j-1] + nModes;

    for ( int j=0; (j<nJobs) ; ++j )
        for ( int m=0 ; (m<nModes) ; ++m )
            x[j][m] = x[0][0] + j*nModes*nTimes + m*nTimes;

    // filling default values
    {
        double *it, *itEnd;
        it = x[0][0];
        itEnd = it + (nJobs*nModes*nTimes);
        for ( ; (it<itEnd) ; it++ )
            *it = 0.0;
    }


    // for all variables
    // filling x from fractional solution
    for ( i=0 ; (i<lp_cols(lp)) ; ++i )
    {
        lp_col_name( origLP, origCols[i], name );
        if (tolower(name[0]=='x'))
        {
            parseName( name, prefix, idx );
            if (prefix[0]=='x')
                x[idx[0]][idx[1]][idx[2]] = xf[i];
        }
    }

    /* x indexes pre-processed problem*/
    int ***xIdx;
    ALLOCATE_VECTOR( xIdx, int **, nJobs );
    ALLOCATE_VECTOR( xIdx[0], int *, nJobs*nModes );
    ALLOCATE_VECTOR( xIdx[0][0], int, nJobs*nModes*nTimes );
    for ( int j=1 ; (j<nJobs) ; ++j )
        xIdx[j] = xIdx[j-1] + nModes;

    for ( int j=0; (j<nJobs) ; ++j )
        for ( int m=0 ; (m<nModes) ; ++m )
            xIdx[j][m] = xIdx[0][0] + j*nModes*nTimes + m*nTimes;
    // filling default values
    {
        int *it, *itEnd;
        it = xIdx[0][0];
        itEnd = it + (nJobs*nModes*nTimes);
        for ( ; (it<itEnd) ; it++ )
            *it = -1;
    }

    /* filling indexes from variables in this pre-processed problem */
    for ( i=0 ; (i<lp_cols(lp)) ; ++i )
    {

        lp_col_name( origLP, origCols[i], name );
        if (tolower(name[0])=='x')
        {
            parseName( name, prefix, idx );
            if (prefix[0]=='x')
                xIdx[idx[0]][idx[1]][idx[2]] = i;
        }
    }

    int ****unitsESM;
    ALLOCATE_VECTOR(unitsESM,int***,nJobs);
    for(int j = 0 ; j < nJobs ; j++)
    {
        const Job* job = Inst_job(inst,j);
        ALLOCATE_VECTOR(unitsESM[j],int**,nModes);
        for(int m = 0 ; m < Job_nModes(job) ; m++)
        {
            const Mode *mode = Job_mode(job,m);
            ALLOCATE_VECTOR(unitsESM[j][m],int*,nTimes);
            for(int s = 0 ; s < nTimes ; s++)
            {
                ALLOCATE_VECTOR_INI(unitsESM[j][m][s],int,nTimes);
                for(int e = s+1 ; e < nTimes ; e++)
                {
                     //ALLOCATE_VECTOR_INI(unitsES[j][m][s][e],int,Inst_nResR(inst) );
                   //  for(int r = 0 ; r < Inst_nResR(inst) ; r++)
                    // {
                     //   int rm = Mode_idxResROnMode(inst,mode,r);
                     //   if(rm != -1){
                            int part1 = e-s;
                            int part2 = MAX(0,Mode_duration(mode)-MAX(0,s-firstTJM[j][m]));
                            int part3 = MAX(0,Mode_duration(mode)-MAX(0,lastTJM[j][m]+Mode_duration(mode)-e));
                            unitsESM[j][m][s][e]= MIN(part1,MIN(part2,part3));
                            //unitsES[j][m][s][e][r] = MIN(part1,MIN(part2,part3));

                       //  }
                    //}
                 }

            }
        }
    }
//rr

    int ****unitsESJ;
    ALLOCATE_VECTOR(unitsESJ,int***,nJobs);
    int ****minMode;
    ALLOCATE_VECTOR(minMode,int***,nJobs);

    for(int j = 0 ; j < nJobs ; j++)
    {
        const Job* job = Inst_job(inst,j);
        ALLOCATE_VECTOR(unitsESJ[j],int**,nTimes);
        ALLOCATE_VECTOR(minMode[j],int**,nTimes);
        for(int s = 0 ; s < nTimes ; s++)
        {
            ALLOCATE_VECTOR(unitsESJ[j][s],int*,nTimes);
            ALLOCATE_VECTOR(minMode[j][s],int*,nTimes);
            for(int e = s+1 ; e < nTimes ; e++)
            {
                 ALLOCATE_VECTOR_INI(unitsESJ[j][s][e],int,Inst_nResR(inst) );
                 ALLOCATE_VECTOR_INI(minMode[j][s][e],int,Inst_nResR(inst));
                 for(int r = 0 ; r < Inst_nResR(inst) ; r++)
                 {

                    int minValue=INT_MAX;
                    for(int m = 0 ; m < Job_nModes(job) ; m++)
                    {
                        const Mode *mode = Job_mode(job,m);
                        int rm = Mode_idxResROnMode(inst,mode,r);
                        int usage = 0;
                        if(rm != -1){
                            usage =  Mode_useResR(mode,rm);
                        }
                        if(minValue > unitsESM[j][m][s][e] * usage){
                            unitsESJ[j][s][e][r] = unitsESM[j][m][s][e];
                            minValue = unitsESM[j][m][s][e] * usage;
                            minMode[j][s][e][r] = usage;
                         //   printf(" r(%d) unitsES[%d][%d][%d][%d][%d] %d, unitsES[j][m][s][e][r] %d * minMode[j][s][e][r] %d coef %d\n", j,j,m,s,e,r, unitsES[j][m][s][e][r] , unitsESJ[j][s][e][r], minMode[j][s][e][r],  (unitsESJ[j][s][e][r]*minMode[j][s][e][r]));
                           //     getchar();
                        }
                    }
                }
            }
        }
    }


// getchar();
    //int s =0, e =5;

    ALLOCATE_VECTOR(cen->cutElem,VecInt*,Inst_nResR(inst)*nTimes*nTimes);
    ALLOCATE_VECTOR(cen->cutCoef,VecDbl*,Inst_nResR(inst)*nTimes*nTimes);

    VecInt *vecIdx = VInt_create();
    VecDbl *vecCo = VDbl_create();
    VecInt *vecNCut = VInt_create();
    VecDbl *vecValue = VDbl_create();
    VecStr *vecNames = VStr_create(256);
    int nCut2 =0;

    for(int r = 0 ; r < Inst_nResR(inst) ; r++)
    {
        for(int s = 0 ; (s < nTimes) ; s++)
        {
            for(int e = s+1 ; e<(nTimes); e++)
            {


                IntDblPair *cutIdxCo;
                ALLOCATE_VECTOR(cutIdxCo, IntDblPair, nCols);
                int c =0;

                int idx2[nCols];
                double coef[nCols];
                double coefr[nCols];
                int leftside[nCols];
                double rhs = ((e-s+1)*Inst_capResR(inst,r));

                int size = 0;
                double slack =0.0;
                        //              printf("s %d, e %d \n", s,e);fflush(stdout);

                for ( int j=0 ; j< nJobs ; j++ )
                {
                    const Job *job = Inst_job( inst, j );
                    const int nModes = Job_nModes( job );
                    coefr[j]=0.0;
                   // printf(" r(%d) coef duration %d * usage %d coef %d\n", j, unitsESJ[j][s][e][r], minMode[j][s][e][r],  (unitsESJ[j][s][e][r]*minMode[j][s][e][r]));
                    if(unitsESJ[j][s][e][r] > 0 && minMode[j][s][e][r] > 0 ){
                        coefr[j] = (unitsESJ[j][s][e][r]*minMode[j][s][e][r]);
                       // printf(" r(%d) coef duration %d * usage %d coef %d\n", j, unitsESJ[j][s][e][r], minMode[j][s][e][r],  (unitsESJ[j][s][e][r]*minMode[j][s][e][r]));
                    }


                    double computSum = 0;
                    int contele = 0;
                    for ( int m=0 ; (m<nModes) ; ++m )
                    {

                        const Mode *mode = Job_mode( job, m );
                       // if(Mode_duration(mode)==0) continue;
                        if(!Mode_isFeasible(inst,mode)) continue;
                        int rm = Mode_idxResROnMode(inst,mode,r);
                        int cjr = 0;
                        if(rm!=-1) cjr = Mode_useResR(mode,rm);;
                        if(cjr<=0) continue;
                        for ( int t=0 ; t< nTimes; t++ )
                        {

                            if(xIdx[j][m][t] == -1 || xf[ xIdx[j][m][t]] <= 0.0 ) continue;
                           // if(xIdx[j][m][t] == -1  ) continue;

                            if(t<=e && t+Mode_duration(mode) >= s  )
                            {

                               int t2 = t+Mode_duration(mode)-1;
                               int value = MAX(0,MIN(t2,e)-MAX(t,s)+1);

/*                                int value=0;
                                if( t>=s)
                                {
                                    if(t+Mode_duration(mode) <= e)
                                    {
                                        value =  Mode_duration(mode);
                                    }
                                    else
                                    {
                                        value =  Mode_duration(mode) - (t+Mode_duration(mode)-1-e);
                                    }
                                }
                                else
                                {
                                    if(t+Mode_duration(mode) <= e)
                                    {
                                        value = Mode_duration(mode) -(s-t);
                                    }
                                    else
                                    {
                                        value = (e-s+1);
                                    }
                                }*/
                                if(value<=0) continue;

                                idx2[size] = xIdx[j][m][t];
                                coef[size] = (value*cjr*xf[ xIdx[j][m][t]]);
                            // computSum+=(value*cjr);//coef[size];
                                computSum+=coef[size];
                                leftside[contele] = size;
                                contele++;
                             //   printf("  xIdx[j][m][t] %d,  x(%d,%d,%d), xf %f,  value %d mode duration %d, cjr %d coef %f \n",  xIdx[j][m][t], j,m,t, xf[ xIdx[j][m][t]],value, Mode_duration(mode), cjr,  coef[size]);
                                size++;
                            }
                        }
                    }


                    if(computSum>=coefr[j]){
                        for(int ie = 0; ie < contele ; ie++){
                          cutIdxCo[c].a = idx2[leftside[ie]];
                          cutIdxCo[c].b = coef[leftside[ie]];
                          slack +=coef[leftside[ie]];
                          c++;
                        }
                    }else{
            //            printf("before rhs %f \n", rhs);
                        rhs-=coefr[j];
                      //  printf("rhs %f coefr[%d] %f\n", rhs, j, coefr[j]); fflush(stdout);
                      //  getchar();
                    }

                }
         //       printf("antes slack %f, rhs %f\n", slack, rhs);
                slack-=rhs;
          //      printf("depois slack %f, rhs %f", slack, rhs);
      //    getchar();
                if(slack > 0.0002 && c>0){

                    cen->cutElem[nCut2] = VInt_create();
                    cen->cutCoef[nCut2] = VDbl_create();
                    cen->nAllocations++;

                    char nameCut[STR_SIZE];
                    sprintf( nameCut, "cutE(%d,%d,%d)#%d#%d", r, s, e, nround, nCut2 );
             //       printf( "\n cutE(%d,%d,%d)#%d#%d", r, s, e, nround, nCut2 );
                    VStr_pushBack(cen->cutname,nameCut);

                            for(int sel  = 0 ; sel <c; sel++)
                            {
                                int elem = cutIdxCo[sel].a;
                                VInt_pushBack(cen->cutElem[nCut2],elem);
                                VDbl_pushBack(cen->cutCoef[nCut2],cutIdxCo[sel].b);//xf[ cutIdxCo[sel].a]);
                       //         printf( " %f * %f: %d \n", cutIdxCo[sel].b,xf[elem],elem);
                            }
                            int *cutIdx = VInt_getPtr(cen->cutElem[nCut2]);
                            double *cutCoef = VDbl_getPtr(cen->cutCoef[nCut2]);
                            VDbl_pushBack(cen->cutrhs,rhs);
                            VInt_pushBack(cen->cutsense,0);
                            VInt_pushBack(cen->cutdominated,0);
                            VInt_pushBack(cen->cutnelem,c);
                            VDbl_pushBack(cen->cutviolation,slack);
                            CutP_quick_sort_vec(cutIdx, cutCoef, c);
                            if((e-s+1) > cen->nMaxWindowWithCutEnergeticViolated )
                            {
                                cen->nMaxWindowWithCutEnergeticViolated = (e-s+1);
                                cen->nMaxWindowWithCutEnergeticViolatedPerc = (double) (e-s+1)/(double) nTimes;
                                cen->nS = s;
                                cen->nE = e;
                            }
                            if((e-s+1) < cen->nMinWindowWithCutEnergeticViolated )
                            {
                                cen->nMinWindowWithCutEnergeticViolated = (e-s+1);
                                cen->nMinWindowWithCutEnergeticViolatedPerc = (double) (e-s+1)/(double) nTimes;
                                cen->nmS = s;
                                cen->nmE = e;
                            }
                      //    printf( " rhs %f, sense %d, nelem %d, viol %f\n",rhs, 0,c,slack);
                            nCut2++;

                   //  getchar();
                }
                    free(cutIdxCo);
            }
        }
    }


     /*   for(int ncdA = 0; ncdA < VDbl_size(cen->cutrhs) ; ncdA++) {
            for(int ncdB = ncdA+1; ncdB < VDbl_size(cen->cutrhs) ; ncdB++) {
                if(VInt_get(cen->cutdominated,ncdA) == 1) break;
                if(VInt_get(cen->cutdominated,ncdB) == 1) continue;
                CutP_dominanceBetweenTwo( lp, ncdA, VInt_getPtr(cen->cutElem[ncdA]), VDbl_getPtr(cen->cutCoef[ncdA]), VDbl_get(cen->cutrhs,ncdA), VInt_get(cen->cutsense,ncdA), VInt_get(cen->cutnelem,ncdA), ncdB, VInt_getPtr(cen->cutElem[ncdB]), VDbl_getPtr(cen->cutCoef[ncdB]), VDbl_get(cen->cutrhs,ncdB), VInt_get(cen->cutsense,ncdB), VInt_get(cen->cutnelem,ncdB),  cen->cutdominated);
                //getchar();
            }
        }
*/



    VInt_free(&vecIdx);
    VInt_free(&vecNCut);
    VDbl_free(&vecCo);
    VDbl_free(&vecValue);
    VStr_free(&vecNames);


    for(int j = 0 ; j < nJobs ; j++)
    {
        const Job* job = Inst_job(inst,j);
        free(firstTJM[j]);
        free(lastTJM[j]);
        for(int m = 0 ; m < Job_nModes(job) ; m++)
        {
            for(int s = 0 ; s < nTimes ; s++)
            {
                free(unitsESM[j][m][s]);
            }
            free(unitsESM[j][m]);
        }
        for(int s = 0 ; s < nTimes ; s++)
        {
            for(int e = s+1 ; e < nTimes ; e++)
            {
                free(unitsESJ[j][s][e]);
                free(minMode[j][s][e]);
            }
            free(unitsESJ[j][s]);
            free(minMode[j][s]);
        }
        free(unitsESJ[j]);
        free(minMode[j]);
        free(unitsESM[j]);
    }
    free(unitsESJ);
    free(minMode);
    free(unitsESM);
    free(firstTJ);
    free(lastTJ);
    free(firstTJM);
    free(lastTJM);

    free(xIdx[0][0]);
    free(xIdx[0]);
    free(xIdx);
    free(x[0][0]);
    free(x[0]);
    free(x);


}



void CutE_add_cuts_energetic_guloso_interval (CutE *cen,  LinearProgram *lp, const int *origCols,  LinearProgram *origLP, const Instance *inst, double timeLeft, int continuous,  int **maxTJM,  double maxcuts, int nround, int s, int e )
{

    int nCols = lp_cols(lp);


    assert(lp);
    assert(inst);
    assert(origLP);

    const double *xf = lp_x(lp);
    char name[256];
    int i;

    /* maximum number of indexes in a var name */
#define MAX_IDX 16
    char prefix[STR_SIZE];
    int idx[MAX_IDX];
    int maxT = 0;

    const int nJobs = Inst_nJobs(inst);
    const int nModes = Inst_nMaxModes(inst);

    // first and last time job has some allocation
    int *firstTJ;
    int *lastTJ;
    ALLOCATE_VECTOR(firstTJ,int,nJobs)
    ALLOCATE_VECTOR(lastTJ,int,nJobs)
    int **firstTJM;
    int **lastTJM;
    ALLOCATE_VECTOR(firstTJM,int*,nJobs)
    ALLOCATE_VECTOR(lastTJM,int*,nJobs)
    for ( int j=0 ; (j<nJobs) ; ++j )
    {
        firstTJ[j]= INT_MAX;
        lastTJ[j]= 0;
        ALLOCATE_VECTOR_INI(firstTJM[j],int,nModes)
        ALLOCATE_VECTOR_INI(lastTJM[j],int,nModes)
        for(int m = 0 ; m < nModes ; m++)
            firstTJM[j][m] = INT_MAX;
    }


    // checking maximum time where there is some allocation
    for ( i=0 ; (i<lp_cols(lp)) ; ++i )
    {
        lp_col_name( origLP, origCols[i], name );
        parseName( name, prefix, idx );

        if (prefix[0]=='x')
        {
            //printf("%s %s", name); fflush(stdout);
            const int j = idx[0];
            const int m = idx[1];
            const int t = idx[2];
            maxT = MAX( maxT, t );
            firstTJM[j][m]= MIN(firstTJM[j][m], t);
            lastTJM[j][m] = MAX(lastTJM[j][m], t);
            firstTJ[j]= MIN(firstTJ[j], t);
            lastTJ[j] = MAX(lastTJ[j], t);
        }
    }

    int nTimes = maxT+1;
    // allocation of x
    double ***x;
    ALLOCATE_VECTOR( x, double **, nJobs );
    ALLOCATE_VECTOR( x[0], double *, nJobs*nModes );
    ALLOCATE_VECTOR( x[0][0], double, nJobs*nModes*nTimes );
    for ( int j=1 ; (j<nJobs) ; ++j )
        x[j] = x[j-1] + nModes;

    for ( int j=0; (j<nJobs) ; ++j )
        for ( int m=0 ; (m<nModes) ; ++m )
            x[j][m] = x[0][0] + j*nModes*nTimes + m*nTimes;

    // filling default values
    {
        double *it, *itEnd;
        it = x[0][0];
        itEnd = it + (nJobs*nModes*nTimes);
        for ( ; (it<itEnd) ; it++ )
            *it = 0.0;
    }


    // for all variables
    // filling x from fractional solution
    for ( i=0 ; (i<lp_cols(lp)) ; ++i )
    {
        lp_col_name( origLP, origCols[i], name );
        if (tolower(name[0]=='x'))
        {
            parseName( name, prefix, idx );
            if (prefix[0]=='x')
                x[idx[0]][idx[1]][idx[2]] = xf[i];
        }
    }

    /* x indexes pre-processed problem*/
    int ***xIdx;
    ALLOCATE_VECTOR( xIdx, int **, nJobs );
    ALLOCATE_VECTOR( xIdx[0], int *, nJobs*nModes );
    ALLOCATE_VECTOR( xIdx[0][0], int, nJobs*nModes*nTimes );
    for ( int j=1 ; (j<nJobs) ; ++j )
        xIdx[j] = xIdx[j-1] + nModes;

    for ( int j=0; (j<nJobs) ; ++j )
        for ( int m=0 ; (m<nModes) ; ++m )
            xIdx[j][m] = xIdx[0][0] + j*nModes*nTimes + m*nTimes;
    // filling default values
    {
        int *it, *itEnd;
        it = xIdx[0][0];
        itEnd = it + (nJobs*nModes*nTimes);
        for ( ; (it<itEnd) ; it++ )
            *it = -1;
    }

    /* filling indexes from variables in this pre-processed problem */
    for ( i=0 ; (i<lp_cols(lp)) ; ++i )
    {

        lp_col_name( origLP, origCols[i], name );
        if (tolower(name[0])=='x')
        {
            parseName( name, prefix, idx );
            if (prefix[0]=='x')
                xIdx[idx[0]][idx[1]][idx[2]] = i;
        }
    }

    int ****unitsESM;
    ALLOCATE_VECTOR(unitsESM,int***,nJobs);
    for(int j = 0 ; j < nJobs ; j++)
    {
        const Job* job = Inst_job(inst,j);
        ALLOCATE_VECTOR(unitsESM[j],int**,nModes);
        for(int m = 0 ; m < Job_nModes(job) ; m++)
        {
            const Mode *mode = Job_mode(job,m);
            ALLOCATE_VECTOR(unitsESM[j][m],int*,nTimes);
            ALLOCATE_VECTOR_INI(unitsESM[j][m][s],int,nTimes);
                            int part1 = e-s;
                            int part2 = MAX(0,Mode_duration(mode)-MAX(0,s-firstTJM[j][m]));
                            int part3 = MAX(0,Mode_duration(mode)-MAX(0,lastTJM[j][m]+Mode_duration(mode)-e));
                            unitsESM[j][m][s][e]= MIN(part1,MIN(part2,part3));


        }
    }
//rr

    int ****unitsESJ;
    ALLOCATE_VECTOR(unitsESJ,int***,nJobs);
    int ****minMode;
    ALLOCATE_VECTOR(minMode,int***,nJobs);

    for(int j = 0 ; j < nJobs ; j++)
    {
        const Job* job = Inst_job(inst,j);
        ALLOCATE_VECTOR(unitsESJ[j],int**,nTimes);
        ALLOCATE_VECTOR(minMode[j],int**,nTimes);

            ALLOCATE_VECTOR(unitsESJ[j][s],int*,nTimes);
            ALLOCATE_VECTOR(minMode[j][s],int*,nTimes);

                 ALLOCATE_VECTOR_INI(unitsESJ[j][s][e],int,Inst_nResR(inst) );
                 ALLOCATE_VECTOR_INI(minMode[j][s][e],int,Inst_nResR(inst));
                 for(int r = 0 ; r < Inst_nResR(inst) ; r++)
                 {

                    int minValue=INT_MAX;
                    for(int m = 0 ; m < Job_nModes(job) ; m++)
                    {
                        const Mode *mode = Job_mode(job,m);
                        int rm = Mode_idxResROnMode(inst,mode,r);
                        int usage = 0;
                        if(rm != -1){
                            usage =  Mode_useResR(mode,rm);
                        }
                        if(minValue > unitsESM[j][m][s][e] * usage){
                            unitsESJ[j][s][e][r] = unitsESM[j][m][s][e];
                            minValue = unitsESM[j][m][s][e] * usage;
                            minMode[j][s][e][r] = usage;
                         //   printf(" r(%d) unitsES[%d][%d][%d][%d][%d] %d, unitsES[j][m][s][e][r] %d * minMode[j][s][e][r] %d coef %d\n", j,j,m,s,e,r, unitsES[j][m][s][e][r] , unitsESJ[j][s][e][r], minMode[j][s][e][r],  (unitsESJ[j][s][e][r]*minMode[j][s][e][r]));
                           //     getchar();
                        }
                    }
                }
    }


// getchar();
    //int s =0, e =5;

    ALLOCATE_VECTOR(cen->cutElem,VecInt*,Inst_nResR(inst)*nTimes*nTimes);
    ALLOCATE_VECTOR(cen->cutCoef,VecDbl*,Inst_nResR(inst)*nTimes*nTimes);

    VecInt *vecIdx = VInt_create();
    VecDbl *vecCo = VDbl_create();
    VecInt *vecNCut = VInt_create();
    VecDbl *vecValue = VDbl_create();
    VecStr *vecNames = VStr_create(256);
    int nCut2 =0;

    for(int r = 0 ; r < Inst_nResR(inst) ; r++)
    {


                IntDblPair *cutIdxCo;
                ALLOCATE_VECTOR(cutIdxCo, IntDblPair, nCols);
                int c =0;

                int idx2[nCols];
                double coef[nCols];
                double coefr[nCols];
                int leftside[nCols];
                double rhs = ((e-s+1)*Inst_capResR(inst,r));

                int size = 0;
                double slack =0.0;
                        //              printf("s %d, e %d \n", s,e);fflush(stdout);

                for ( int j=0 ; j< nJobs ; j++ )
                {
                    const Job *job = Inst_job( inst, j );
                    const int nModes = Job_nModes( job );
                    coefr[j]=0.0;
                   // printf(" r(%d) coef duration %d * usage %d coef %d\n", j, unitsESJ[j][s][e][r], minMode[j][s][e][r],  (unitsESJ[j][s][e][r]*minMode[j][s][e][r]));
                    if(unitsESJ[j][s][e][r] > 0 && minMode[j][s][e][r] > 0 ){
                        coefr[j] = (unitsESJ[j][s][e][r]*minMode[j][s][e][r]);
                       // printf(" r(%d) coef duration %d * usage %d coef %d\n", j, unitsESJ[j][s][e][r], minMode[j][s][e][r],  (unitsESJ[j][s][e][r]*minMode[j][s][e][r]));
                    }


                    double computSum = 0;
                    int contele = 0;
                    for ( int m=0 ; (m<nModes) ; ++m )
                    {

                        const Mode *mode = Job_mode( job, m );
                       // if(Mode_duration(mode)==0) continue;
                        if(!Mode_isFeasible(inst,mode)) continue;
                        int rm = Mode_idxResROnMode(inst,mode,r);
                        int cjr = 0;
                        if(rm!=-1) cjr = Mode_useResR(mode,rm);;
                        if(cjr<=0) continue;
                        for ( int t=0 ; t< nTimes; t++ )
                        {

                            if(xIdx[j][m][t] == -1 || xf[ xIdx[j][m][t]] <= 0.0 ) continue;
                           // if(xIdx[j][m][t] == -1  ) continue;

                            if(t<=e && t+Mode_duration(mode) >= s  )
                            {

                               int t2 = t+Mode_duration(mode)-1;
                               int value = MAX(0,MIN(t2,e)-MAX(t,s)+1);

/*                                int value=0;
                                if( t>=s)
                                {
                                    if(t+Mode_duration(mode) <= e)
                                    {
                                        value =  Mode_duration(mode);
                                    }
                                    else
                                    {
                                        value =  Mode_duration(mode) - (t+Mode_duration(mode)-1-e);
                                    }
                                }
                                else
                                {
                                    if(t+Mode_duration(mode) <= e)
                                    {
                                        value = Mode_duration(mode) -(s-t);
                                    }
                                    else
                                    {
                                        value = (e-s+1);
                                    }
                                }*/
                                if(value<=0) continue;

                                idx2[size] = xIdx[j][m][t];
                                coef[size] = (value*cjr*xf[ xIdx[j][m][t]]);
                            // computSum+=(value*cjr);//coef[size];
                                computSum+=coef[size];
                                leftside[contele] = size;
                                contele++;
                             //   printf("  xIdx[j][m][t] %d,  x(%d,%d,%d), xf %f,  value %d mode duration %d, cjr %d coef %f \n",  xIdx[j][m][t], j,m,t, xf[ xIdx[j][m][t]],value, Mode_duration(mode), cjr,  coef[size]);
                                size++;
                            }
                        }
                    }


                    if(computSum>=coefr[j]){
                        for(int ie = 0; ie < contele ; ie++){
                          cutIdxCo[c].a = idx2[leftside[ie]];
                          cutIdxCo[c].b = coef[leftside[ie]];
                          slack +=coef[leftside[ie]];
                          c++;
                        }
                    }else{
            //            printf("before rhs %f \n", rhs);
                        rhs-=coefr[j];
                      //  printf("rhs %f coefr[%d] %f\n", rhs, j, coefr[j]); fflush(stdout);
                      //  getchar();
                    }

                }
         //       printf("antes slack %f, rhs %f\n", slack, rhs);
                slack-=rhs;
          //      printf("depois slack %f, rhs %f", slack, rhs);
      //    getchar();
                if(c>0){

                    cen->cutElem[nCut2] = VInt_create();
                    cen->cutCoef[nCut2] = VDbl_create();
                    cen->nAllocations++;

                    char nameCut[STR_SIZE];
                    sprintf( nameCut, "cutE(%d,%d,%d)#%d#%d", r, s, e, nround, nCut2 );
                    //printf( "\n cutE(%d,%d,%d)#%d#%d", r, s, e, nround, nCut2 );
                    VStr_pushBack(cen->cutname,nameCut);

                            for(int sel  = 0 ; sel <c; sel++)
                            {
                                int elem = cutIdxCo[sel].a;
                                VInt_pushBack(cen->cutElem[nCut2],elem);
                                VDbl_pushBack(cen->cutCoef[nCut2],cutIdxCo[sel].b);//xf[ cutIdxCo[sel].a]);
                              //  printf( " %f * %f: %d \n", cutIdxCo[sel].b,xf[elem],elem);
                            }
                            int *cutIdx = VInt_getPtr(cen->cutElem[nCut2]);
                            double *cutCoef = VDbl_getPtr(cen->cutCoef[nCut2]);
                            VDbl_pushBack(cen->cutrhs,rhs);
                            VInt_pushBack(cen->cutsense,0);
                            VInt_pushBack(cen->cutdominated,0);
                            VInt_pushBack(cen->cutnelem,c);
                            VDbl_pushBack(cen->cutviolation,slack);
                            CutP_quick_sort_vec(cutIdx, cutCoef, c);
                            if((e-s+1) > cen->nMaxWindowWithCutEnergeticViolated )
                            {
                                cen->nMaxWindowWithCutEnergeticViolated = (e-s+1);
                                cen->nMaxWindowWithCutEnergeticViolatedPerc = (double) (e-s+1)/(double) nTimes;
                                cen->nS = s;
                                cen->nE = e;
                            }
                            if((e-s+1) < cen->nMinWindowWithCutEnergeticViolated )
                            {
                                cen->nMinWindowWithCutEnergeticViolated = (e-s+1);
                                cen->nMinWindowWithCutEnergeticViolatedPerc = (double) (e-s+1)/(double) nTimes;
                                cen->nmS = s;
                                cen->nmE = e;
                            }
                        //  printf( " rhs %f, sense %d, nelem %d, viol %f\n",rhs, 0,c,slack);
                            nCut2++;

                   //  getchar();
                }
                    free(cutIdxCo);

    }


     /*   for(int ncdA = 0; ncdA < VDbl_size(cen->cutrhs) ; ncdA++) {
            for(int ncdB = ncdA+1; ncdB < VDbl_size(cen->cutrhs) ; ncdB++) {
                if(VInt_get(cen->cutdominated,ncdA) == 1) break;
                if(VInt_get(cen->cutdominated,ncdB) == 1) continue;
                CutP_dominanceBetweenTwo( lp, ncdA, VInt_getPtr(cen->cutElem[ncdA]), VDbl_getPtr(cen->cutCoef[ncdA]), VDbl_get(cen->cutrhs,ncdA), VInt_get(cen->cutsense,ncdA), VInt_get(cen->cutnelem,ncdA), ncdB, VInt_getPtr(cen->cutElem[ncdB]), VDbl_getPtr(cen->cutCoef[ncdB]), VDbl_get(cen->cutrhs,ncdB), VInt_get(cen->cutsense,ncdB), VInt_get(cen->cutnelem,ncdB),  cen->cutdominated);
                //getchar();
            }
        }
*/



    VInt_free(&vecIdx);
    VInt_free(&vecNCut);
    VDbl_free(&vecCo);
    VDbl_free(&vecValue);
    VStr_free(&vecNames);


    for(int j = 0 ; j < nJobs ; j++)
    {
        const Job* job = Inst_job(inst,j);
        free(firstTJM[j]);
        free(lastTJM[j]);
        for(int m = 0 ; m < Job_nModes(job) ; m++)
        {

            free(unitsESM[j][m][s]);
            free(unitsESM[j][m]);
        }
            free(unitsESJ[j][s][e]);
            free(minMode[j][s][e]);
            free(unitsESJ[j][s]);
            free(minMode[j][s]);

        free(unitsESJ[j]);
        free(minMode[j]);
        free(unitsESM[j]);
    }
    free(unitsESJ);
    free(minMode);
    free(unitsESM);
    free(firstTJ);
    free(lastTJ);
    free(firstTJM);
    free(lastTJM);

    free(xIdx[0][0]);
    free(xIdx[0]);
    free(xIdx);
    free(x[0][0]);
    free(x[0]);
    free(x);


}
