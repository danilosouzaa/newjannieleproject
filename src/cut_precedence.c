
/* GOAL PSP Solver II
 * Resource-Constrained Project Scheduling Problems (RCPSP)
 *
 * Develop as part of the D.Sc. thesis of Araujo, Janniele A. S., with collaboration
 *                                   of Santos, H.G.
 */

#include "cut_precedence.h"
#include "lp.h"
#define VERBOSE  1

/*create and initialize the structure of cutPR*/
CutPR *CutPR_create(const Instance *inst)
{

    CutPR *cpr;
    ALLOCATE_INI(cpr,CutPR);
    cpr->nAllocations = 0;
    cpr->cutrhs = VDbl_create();
    cpr->cutviolation = VDbl_create();
    cpr->cutnelem = VInt_create();
    cpr->cutsense = VInt_create();
    cpr->cutname = VStr_create(256);
    cpr->cutdominated = VInt_create();

    cpr->inst = inst;

    return cpr;
}


double computeSumPred( const Instance *inst,  int *nz, const double coef,  int j, int ***TJ, int **nCont, int start, int end, int t, int s, const double ***x, const int ***xIdx, IntDblPair *cutIdxCo, int lifting );

double computeSumSucc( const Instance *inst,  int *nz, const double coef, int j, int ***TJ, int **nCont, int start, int end, int t, int s, const double ***x, const int ***xIdx, IntDblPair *cutIdxCo );

double computeSumPred( const Instance *inst,  int *nz, const double coef,  int j, int ***TJ, int **nCont, int start, int end, int t, int s, const double ***x, const int ***xIdx, IntDblPair *cutIdxCo, int lifting )
{

    const Job *job = Inst_job( inst, j );
    const int nM = Job_nModes( job );
    double r = 0.0;
    for ( int m=0 ; (m<nM) ; ++m ) {
        //int endt = MIN(end, t+Inst_getMaxDIJ(inst, j, s)-Inst_getMaxDIJM(inst,j,s, m));
        int endt;
        if(lifting)
            endt = t+Inst_getMaxDIJ(inst, j, s)-Inst_getMaxDIJM(inst,j,s, m);
        else
            endt = t+Inst_getMaxDIJ(inst, j, s);

        for ( int tc=start; (tc<=endt) ; ++tc ) {
            if (xIdx[j][m][tc]!=-1) {
                r += x[j][m][tc];
                cutIdxCo[(*nz)].a = xIdx[j][m][tc];
                cutIdxCo[(*nz)].b = coef;
                //  printf(" +%f (j %d, m %d, t %d) id %d ", j,m,tc, coef, xIdx[j][m][tc]);
                (*nz)++;
            }
        }
    }
    return r;
}

double computeSumSucc( const Instance *inst, int *nz, const double coef,  int j, int ***TJ, int **nCont, int start, int end, int t, int s, const double ***x, const int ***xIdx, IntDblPair *cutIdxCo )
{

    const Job *job = Inst_job( inst, s );
    const int nM = Job_nModes( job );
    double r = 0.0;

    int iniT = MIN(end,t+Inst_getMaxDIJ(inst, j, s));
    for ( int m=0 ; (m<nM) ; ++m ) {
        for ( int tc=start+1 ; (tc<=iniT) ; ++tc ) {
            if (xIdx[s][m][tc]!=-1) {
                r += x[s][m][tc];
                cutIdxCo[(*nz)].a = xIdx[s][m][tc];
                cutIdxCo[(*nz)].b = coef;
                //  printf(" %f (s %d, m %d, t %d) id %d ", s,m,tc, coef, xIdx[s][m][tc]);
                (*nz)++;
            }
        }
    }

    return r;
}

/*procedure to separate precedence cuts*/
void CutPR_add_cuts_precedence_parallel( CutPR *cpr, LinearProgram *lp, const int *origCols, LinearProgram *origLP, const Instance *inst, double timeLeft, int continuous, int lifting, double maxcuts, int nround )
{
    //double iniPC = omp_get_wtime();

    //    double tinit = omp_get_wtime();
    assert(inst);
    int nCut = lp_rows(lp);

    const double *xf = lp_x(lp);
    char name[256];
    int i;

#define MAX_IDX 16
    char prefix[STR_SIZE];
    int idx[MAX_IDX];
    int maxT = 0;

    const int nJobs = Inst_nJobs(inst);
    const int nModes = Inst_nMaxModes(inst);

    int ***TJ;
    ALLOCATE_VECTOR( TJ, int **, nJobs );
    for(int tji = 0 ; tji <nJobs ; tji++) {
        ALLOCATE_VECTOR(TJ[tji], int*, nModes);
        for(int m = 0 ; m <nModes ; m++)
            ALLOCATE_VECTOR_INI(TJ[tji][m], int, 100000);
    }
    int **cont;
    ALLOCATE_VECTOR( cont, int*, nJobs );
    for(int tji = 0 ; tji <nJobs ; tji++)
        ALLOCATE_VECTOR_INI(cont[tji], int, nModes);

    int **TJJ;
    ALLOCATE_VECTOR( TJJ, int*, nJobs );
    for(int tji = 0 ; tji <nJobs ; tji++)
        ALLOCATE_VECTOR_INI(TJJ[tji], int, 100000);

    int *contJ;
    ALLOCATE_VECTOR_INI( contJ, int, nJobs );

    int *firstTJ;
    ALLOCATE_VECTOR_INI( firstTJ, int, nJobs );
    int *lastTJ;
    ALLOCATE_VECTOR_INI( lastTJ, int, nJobs );

    for ( i=0 ; (i<lp_cols(lp)) ; ++i ) {
        lp_col_name( origLP, origCols[i], name );
        parseName( name, prefix, idx );
        if (prefix[0]=='x') {
            int j = idx[0];
            const int t = idx[2];
            const int m = idx[1];
            maxT = MAX( maxT, t );
            TJ[j][m][cont[j][m]++] = t;
            TJJ[j][contJ[j]++] = t;
            firstTJ[j] = MIN(firstTJ[j], t);
            lastTJ[j] = MAX(lastTJ[j], t);

        }
    }

    int nTimes = maxT+1;

    double ***x;
    ALLOCATE_VECTOR( x, double **, nJobs );
    ALLOCATE_VECTOR( x[0], double *, nJobs*nModes );
    ALLOCATE_VECTOR( x[0][0], double, nJobs*nModes*nTimes );
    for ( int j=1 ; (j<nJobs) ; ++j )
        x[j] = x[j-1] + nModes;

    for ( int j=0; (j<nJobs) ; ++j )
        for ( int m=0 ; (m<nModes) ; ++m )
            x[j][m] = x[0][0] + j*nModes*nTimes + m*nTimes;
    {
        double *it, *itEnd;
        it = x[0][0];
        itEnd = it + (nJobs*nModes*nTimes);
        for ( ; (it<itEnd) ; it++ )
            *it = 0.0;
    }

    for ( i=0 ; (i<lp_cols(lp)) ; ++i ) {

        lp_col_name( origLP, origCols[i], name );
        if (tolower(name[0]=='x')) {
            parseName( name, prefix, idx );
            if (prefix[0]=='x')
                x[idx[0]][idx[1]][idx[2]] = xf[i];
        }
    }

    int ***xIdx;
    ALLOCATE_VECTOR( xIdx, int **, nJobs );
    ALLOCATE_VECTOR( xIdx[0], int *, nJobs*nModes );
    ALLOCATE_VECTOR( xIdx[0][0], int, nJobs*nModes*nTimes );
    for ( int j=1 ; (j<nJobs) ; ++j )
        xIdx[j] = xIdx[j-1] + nModes;

    for ( int j=0; (j<nJobs) ; ++j )
        for ( int m=0 ; (m<nModes) ; ++m )
            xIdx[j][m] = xIdx[0][0] + j*nModes*nTimes + m*nTimes;
    {
        int *it, *itEnd;
        it = xIdx[0][0];
        itEnd = it + (nJobs*nModes*nTimes);
        for ( ; (it<itEnd) ; it++ )
            *it = -1;
    }

    for ( i=0 ; (i<lp_cols(lp)) ; ++i ) {
        lp_col_name( origLP, origCols[i], name );
        if (tolower(name[0])=='x') {
            parseName( name, prefix, idx );
            if (prefix[0]=='x')
                xIdx[idx[0]][idx[1]][idx[2]] = i;
        }
    }

    int nCols = lp_cols(lp);

    VecInt *vecIdx = VInt_create();
    VecDbl *vecCo = VDbl_create();
    VecInt *vecNCut = VInt_create();
    VecDbl *vecValue = VDbl_create();
    VecStr *vecNames = VStr_create(256);

    int nc = 0;
    for ( int ij=0 ; (ij<nJobs) ; ij++ ) {

        const Job *job = Inst_job(inst,ij);

        if (Job_minDuration(job)==0)
            continue;

        for ( int iss=0; (iss<Inst_getSizePath(inst,ij)) ; ++iss ) {

            int is = Inst_getValuePosPath(inst,ij,iss);
            if(is==ij) continue;

            const Job *jSucc = Inst_job( inst, is );
            if (Job_minDuration(jSucc)==0)
                continue;

            int tmin = MAX(firstTJ[ij], firstTJ[is]-Inst_getMaxDIJ(inst, ij, is));
            int tmax = MIN(lastTJ[ij], lastTJ[is]-Inst_getMaxDIJ(inst, ij, is)-1);
            int aux = 0;
            //int  t = tmin;
            for ( int t= tmin ; (t<=tmax) ; ++t ) {

                int cutNZ = 0;
                IntDblPair *cutIdxCo;
                ALLOCATE_VECTOR(cutIdxCo,IntDblPair,nCols);

                const double sumPred = computeSumPred( inst, &cutNZ, 1.0, ij, TJ, cont, firstTJ[ij], lastTJ[ij], t, is, (const double ***) x,  (const int ***) xIdx, cutIdxCo, lifting );
                // printf("\n: sumPred %f \n", sumPred);
                if(cutNZ == 0 ) {
                    free(cutIdxCo);
                    continue;
                }
                int contAux = cutNZ;
                const double sumSucc = computeSumSucc( inst, &cutNZ, -1.0, ij, TJ, cont, firstTJ[is],lastTJ[is],  t, is, (const double ***) x,  (const int ***) xIdx, cutIdxCo );
                //  printf("\n: sumSuc %f \n", sumSucc);
                if(cutNZ-contAux == 0 ) {
                    free(cutIdxCo);
                    continue;
                }
                //    printf("\n \n") ; getchar();
                //if (sumPred-sumSucc > 0.0005 ) {
                //if (sumPred-sumSucc < 0.0005 ) {
                if (sumPred < (sumSucc-0.0005)) {

                    char *name;
                    ALLOCATE_VECTOR(name,char,256);
                    sprintf( name, "cutPrec(%d,%d,%d)#%d#%d", ij,is,t,nround, lp_rows(lp)+nCut++);
                    VStr_pushBack(vecNames,name);
                    VInt_pushBack(vecNCut,cutNZ);
                    //  printf("%s\n", VStr_get(vecNames[nc],0));
                    //VInt_pushBack(vecIdx,cutNZ);
                    //VDbl_pushBack(vecCo, (double)cutNZ);
                    VInt_pushBack(vecNCut,VInt_size(vecIdx));
                    for(int idxc = 0 ; idxc < cutNZ ; idxc++) {
                        VInt_pushBack(vecIdx,cutIdxCo[idxc].a);
                        VDbl_pushBack(vecCo, cutIdxCo[idxc].b);
                        //    printf("%d * %f ", VInt_get(vecIdx[nc],idxc), VDbl_get(vecCo[nc],idxc));
                    }
                    // getchar();

                    double valueViol =  (sumPred-sumSucc)*-1 ;
                    // double valueViol =  (sumPred-sumSucc);

                    VDbl_pushBack(vecValue,  valueViol);
                    nc++;
                    free(name);
                    aux = 1;
                }
                free(cutIdxCo);
            }
            if(aux) {
                int max = Inst_getSizePath(inst,ij);
                while(is != ij && iss < max-1 ) {
                    //  printf("max %d ij %d is %d, iss %d\n", max, ij, is, iss);
                    iss++;
                    // printf("iss %d\n", iss);
                    is = Inst_getValuePosPath(inst,ij,iss);
                }
            }
        }
    }


    int ncut2 = 0;
    if(nc>0) {
        int *idValueCut;
        double *valueCute;
        ALLOCATE_VECTOR_INI(idValueCut,int,nc);
        ALLOCATE_VECTOR_INI(valueCute,double,nc);

        //  IntDblPair *valueCutIdxCO;
        // ALLOCATE_VECTOR(valueCutIdxCO,IntDblPair,nc);
        int icP = 0;
        if (maxcuts<1) maxcuts = maxcuts*lp_rows(lp);

        if(VERBOSE==2) printf("maxcuts %f. ", maxcuts);


        ALLOCATE_VECTOR(cpr->cutElem,VecInt*,nc);
        ALLOCATE_VECTOR(cpr->cutCoef,VecDbl*,nc);
        //REALLOCATE_VECTOR(cpr->cutElem,VecInt*,nc);
        //REALLOCATE_VECTOR(cpr->cutCoef,VecDbl*,nc);
        int icc2 = 0;
        for(icc2 = 0 ; icc2 < nc && icc2 < maxcuts; icc2++) {
            idValueCut[icc2] = icc2;
            valueCute[icc2] = VDbl_get(vecValue,icc2);
            //  if(icc>0) {
            cpr->cutElem[icc2] = VInt_create();
            cpr->cutCoef[icc2] = VDbl_create();
            cpr->nAllocations++;
            // }
        }
        //  printf("MipC_quick_sort_vec_by_double: start\n");fflush(stdout);
        CutP_quick_sort_vec_by_double ( idValueCut, valueCute, icc2);
        // qsort(valueCutIdxCO, nc, sizeof(IntDblPair), cmp_int_dbl_pair_b);
        //printf("MipC_quick_sort_vec_by_double: end\n");
        if(VERBOSE==2) printf("\n%d precedence cuts were found. ", nc);
        fflush(stdout);
        int *elem = VInt_getPtr(vecIdx);
        double *coef = VDbl_getPtr(vecCo);

        for(int icc = 0 ; icc < nc && icc < maxcuts; icc++) {

            int ic = idValueCut[icc];
            //  printf("nCut %d ", VInt_get(vecNCut,icP));
            int sizeCut = VInt_get(vecNCut,icP++);
            int posCut = VInt_get(vecNCut,icP++);

            //   printf("\n");
            for(int ec = 0 ; ec < sizeCut ; ec++ ) {
                VInt_pushBack(cpr->cutElem[icc],elem[posCut]);
                VDbl_pushBack(cpr->cutCoef[icc],coef[posCut]);
                //   printf("%f*%d ",VDbl_get(cpr->cutCoef[icc],ec),VInt_get(cpr->cutElem[icc],ec));
                posCut++;
            }
            int *cutIdx = VInt_getPtr(cpr->cutElem[icc]);
            double *cutCoef = VDbl_getPtr(cpr->cutCoef[icc]);
            VStr_pushBack(cpr->cutname, VStr_get(vecNames,ic));
            VDbl_pushBack(cpr->cutrhs, 0.0);
            VInt_pushBack(cpr->cutsense, 2);
            VInt_pushBack(cpr->cutdominated,0);
            VInt_pushBack(cpr->cutnelem, sizeCut);
            VDbl_pushBack(cpr->cutviolation, valueCute[ic]);

            // printf("\nSENSE %d RHS %f VI %f\n", 2,0.0, VDbl_get(cpr->cutviolation,icc));fflush(stdout);
            //   getchar();
            CutP_quick_sort_vec(cutIdx, cutCoef, sizeCut);
            ncut2++;
            nCut++;
        }
        free(valueCute);
        free(idValueCut);
    }
    // cpr->nAddCuts=nCut;

    /*   double timeinitdomresource = omp_get_wtime();
       for(int ncdA = 0; ncdA < ncut2 ; ncdA++){
           for(int ncdB = ncdA+1; ncdB < ncut2 ; ncdB++){
               if(VInt_get(cpr->cutdominated,ncdA) == 1) break;
               if(VInt_get(cpr->cutdominated,ncdB) == 1) continue;
               CutP_dominanceBetweenTwo( lp, ncdA, VInt_getPtr(cpr->cutElem[ncdA]), VDbl_getPtr(cpr->cutCoef[ncdA]), VDbl_get(cpr->cutrhs,ncdA), VInt_get(cpr->cutsense,ncdA), VInt_get(cpr->cutnelem,ncdA), ncdB, VInt_getPtr(cpr->cutElem[ncdB]), VDbl_getPtr(cpr->cutCoef[ncdB]), VDbl_get(cpr->cutrhs,ncdB), VInt_get(cpr->cutsense,ncdB), VInt_get(cpr->cutnelem,ncdB),  cpr->cutdominated);
               //   getchar();
           }
       }
       printf("time compute dominance precedence %f\n",omp_get_wtime()-timeinitdomresource);
    */
    VInt_free(&vecIdx);
    VInt_free(&vecNCut);
    VDbl_free(&vecCo);
    VDbl_free(&vecValue);
    VStr_free(&vecNames);

    free(xIdx[0][0]);
    free(xIdx[0]);
    free(xIdx);
    free(x[0][0]);
    free(x[0]);
    free(x);

    for(int tji = 0 ; tji <nJobs ; tji++) {
        free(cont[tji]);
        for(int m = 0 ; m <nModes ; m++)
            free(TJ[tji][m]);
        free(TJ[tji]);
        free(TJJ[tji]);
    }
    free(TJ);
    free(TJJ);
    free(cont);
    free(contJ);

    free(firstTJ);
    free(lastTJ);
    //double endPC = omp_get_wtime()- iniPC;
    //printf("timePC %f \n", endPC); fflush(stdout);


    //  cpr->timeseparation = (double) omp_get_wtime()-tinit;

}


/*procedure to insert precedence cuts at callback*/
void CutPR_add_cuts_precedence_callback(  GRBmodel *model,   void *cbdata, const Instance *inst, double *timeLeft, int lifting, double maxcuts, int  *nCuts)
{
    double tinit = omp_get_wtime();
    CutPR *cpr = CutPR_create(inst);

    int numvars=0;
    int numrows=0;


    int error = GRBgetintattr(model, GRB_INT_ATTR_NUMVARS, &numvars);
    error = GRBgetintattr(model, GRB_INT_ATTR_NUMCONSTRS, &numrows);
    double *xf = malloc(numvars*sizeof(double));
    char **varname;

    GRBcbget(cbdata,GRB_CB_MIPNODE,GRB_CB_MIPNODE_REL,xf);

    varname = malloc(numvars*sizeof(char*));

    for (int j = 0; j < numvars; ++j ) {
        error = GRBgetstrattrelement(model, GRB_STR_ATTR_VARNAME, j, &varname[j]);
    //    error = GRBgetdblattrelement(model, GRB_DBL_ATTR_X, j, &xf[j]);
       // if (xf[j] > 0.00000) {
      //  printf("%s %f\n", varname[j], xf[j]);
       // }
    }
    //getchar();


    int i;

#define MAX_IDX 16
    char prefix[STR_SIZE];
    int idx[MAX_IDX];
    int maxT = 0;

    int nJobs = Inst_nJobs(inst);
    int nModes = Inst_nMaxModes(inst);

   // printf("nJobs %d", nJobs); getchar();
    int ***TJ;
    ALLOCATE_VECTOR( TJ, int **, nJobs );
    for(int tji = 0 ; tji <nJobs ; tji++) {
        ALLOCATE_VECTOR(TJ[tji], int*, nModes);
        for(int m = 0 ; m <nModes ; m++)
            ALLOCATE_VECTOR_INI(TJ[tji][m], int, 100000);
    }
    int **cont;
    ALLOCATE_VECTOR( cont, int*, nJobs );
    for(int tji = 0 ; tji <nJobs ; tji++)
        ALLOCATE_VECTOR_INI(cont[tji], int, nModes);

    int **TJJ;
    ALLOCATE_VECTOR( TJJ, int*, nJobs );
    for(int tji = 0 ; tji <nJobs ; tji++)
        ALLOCATE_VECTOR_INI(TJJ[tji], int, 100000);

    int *contJ;
    ALLOCATE_VECTOR_INI( contJ, int, nJobs );

    int *firstTJ;
    ALLOCATE_VECTOR_INI( firstTJ, int, nJobs );
    int *lastTJ;
    ALLOCATE_VECTOR_INI( lastTJ, int, nJobs );

    for ( i=0 ; (i<numvars) ; ++i ) {
        parseName( varname[i], prefix, idx );
       // printf("%s %f\n", varname[i], xf[i]);
        if (prefix[0]=='x') {
            int j = idx[0];
            const int t = idx[2];
            const int m = idx[1];
            maxT = MAX( maxT, t );
            TJ[j][m][cont[j][m]++] = t;
            TJJ[j][contJ[j]++] = t;
            firstTJ[j] = MIN(firstTJ[j], t);
            lastTJ[j] = MAX(lastTJ[j], t);

        }
    }

    int nTimes = maxT+1;

    double ***x;
    ALLOCATE_VECTOR( x, double **, nJobs );
    ALLOCATE_VECTOR( x[0], double *, nJobs*nModes );
    ALLOCATE_VECTOR( x[0][0], double, nJobs*nModes*nTimes );
    for ( int j=1 ; (j<nJobs) ; ++j )
        x[j] = x[j-1] + nModes;

    for ( int j=0; (j<nJobs) ; ++j )
        for ( int m=0 ; (m<nModes) ; ++m )
            x[j][m] = x[0][0] + j*nModes*nTimes + m*nTimes;
    {
        double *it, *itEnd;
        it = x[0][0];
        itEnd = it + (nJobs*nModes*nTimes);
        for ( ; (it<itEnd) ; it++ )
            *it = 0.0;
    }



    int ***xIdx;
    ALLOCATE_VECTOR( xIdx, int **, nJobs );
    ALLOCATE_VECTOR( xIdx[0], int *, nJobs*nModes );
    ALLOCATE_VECTOR( xIdx[0][0], int, nJobs*nModes*nTimes );
    for ( int j=1 ; (j<nJobs) ; ++j )
        xIdx[j] = xIdx[j-1] + nModes;

    for ( int j=0; (j<nJobs) ; ++j )
        for ( int m=0 ; (m<nModes) ; ++m )
            xIdx[j][m] = xIdx[0][0] + j*nModes*nTimes + m*nTimes;
    {
        int *it, *itEnd;
        it = xIdx[0][0];
        itEnd = it + (nJobs*nModes*nTimes);
        for ( ; (it<itEnd) ; it++ )
            *it = -1;
    }

    for ( i=0 ; (i<numvars) ; ++i ) {
        parseName( varname[i], prefix, idx );
        if (tolower(varname[i][0])=='x') {
            parseName( varname[i], prefix, idx );
            if (prefix[0]=='x'){
                xIdx[idx[0]][idx[1]][idx[2]] = i;
                x[idx[0]][idx[1]][idx[2]] = xf[i];
            }
        }
    }

    free(varname);

    VecInt *vecIdx = VInt_create();
    VecDbl *vecCo = VDbl_create();
    VecInt *vecNCut = VInt_create();
    VecDbl *vecValue = VDbl_create();

    int nc = 0;
    for ( int ij=0 ; (ij<nJobs) ; ij++ ) {

        const Job *job = Inst_job(inst,ij);

        if (Job_minDuration(job)==0)
            continue;

        for ( int iss=0; (iss<Inst_getSizePath(inst,ij)) ; ++iss ) {

            int is = Inst_getValuePosPath(inst,ij,iss);
            if(is==ij) continue;

            const Job *jSucc = Inst_job( inst, is );
            if (Job_minDuration(jSucc)==0)
                continue;

            int tmin = MAX(firstTJ[ij], firstTJ[is]-Inst_getMaxDIJ(inst, ij, is));
            int tmax = MIN(lastTJ[ij], lastTJ[is]-Inst_getMaxDIJ(inst, ij, is)-1);
            int aux = 0;
            //int  t = tmin;
            for ( int t= tmin ; (t<=tmax) ; ++t ) {

                int cutNZ = 0;
                IntDblPair *cutIdxCo;
                ALLOCATE_VECTOR(cutIdxCo,IntDblPair,numvars);

                double sumPred = computeSumPred( inst, &cutNZ, 1.0, ij, TJ, cont, firstTJ[ij], lastTJ[ij], t, is, (const double ***) x,  (const int ***) xIdx, cutIdxCo, lifting );
             //    printf("\n: sumPred %f \n", sumPred);
                if(cutNZ == 0 ) {
                    free(cutIdxCo);
                    continue;
                }
                int contAux = cutNZ;
                double sumSucc = computeSumSucc( inst, &cutNZ, -1.0, ij, TJ, cont, firstTJ[is],lastTJ[is],  t, is, (const double ***) x,  (const int ***) xIdx, cutIdxCo );
             //     printf("\n: sumSuc %f \n", sumSucc);
                if(cutNZ-contAux == 0 ) {
                    free(cutIdxCo);
                    continue;
                }
                //    printf("\n \n") ; getchar();
                //if (sumPred-sumSucc > 0.0005 ) {
                //if (sumPred-sumSucc < 0.0005 ) {
                if (sumPred < (sumSucc-0.0005)) {

                    VInt_pushBack(vecNCut,cutNZ);
                    VInt_pushBack(vecNCut,VInt_size(vecIdx));
                    for(int idxc = 0 ; idxc < cutNZ ; idxc++) {
                        VInt_pushBack(vecIdx,cutIdxCo[idxc].a);
                        VDbl_pushBack(vecCo, cutIdxCo[idxc].b);
                          //  printf("%d * %f ", VInt_get(vecIdx,idxc), VDbl_get(vecCo,idxc));
                    }
                    // getchar();

                    double valueViol =  (sumPred-sumSucc)*-1 ;
                    // double valueViol =  (sumPred-sumSucc);

                    VDbl_pushBack(vecValue,  valueViol);
                    nc++;
                    aux = 1;
                }
                free(cutIdxCo);
            }
            if(aux) {
                int max = Inst_getSizePath(inst,ij);
                while(is != ij && iss < max-1 ) {
                    //  printf("max %d ij %d is %d, iss %d\n", max, ij, is, iss);
                    iss++;
                    // printf("iss %d\n", iss);
                    is = Inst_getValuePosPath(inst,ij,iss);
                }
            }
        }
    }


    int ncut2 = 0;
    if(nc>0) {
        int *idValueCut;
        double *valueCute;
        ALLOCATE_VECTOR_INI(idValueCut,int,nc);
        ALLOCATE_VECTOR_INI(valueCute,double,nc);

        //  IntDblPair *valueCutIdxCO;
        // ALLOCATE_VECTOR(valueCutIdxCO,IntDblPair,nc);
        int icP = 0;

        if (maxcuts<1) maxcuts = maxcuts*numrows;

        if(VERBOSE==2) printf("maxcuts %f. numrows %d", maxcuts, numrows);
     //   printf("maxcuts %f. numrows %d", maxcuts, numrows); //getchar();

        ALLOCATE_VECTOR(cpr->cutElem,VecInt*,nc);
        ALLOCATE_VECTOR(cpr->cutCoef,VecDbl*,nc);
        //REALLOCATE_VECTOR(cpr->cutElem,VecInt*,nc);
        //REALLOCATE_VECTOR(cpr->cutCoef,VecDbl*,nc);
        int icc2 = 0;
        for(icc2 = 0 ; icc2 < nc && icc2 < maxcuts; icc2++) {
            idValueCut[icc2] = icc2;
            valueCute[icc2] = VDbl_get(vecValue,icc2);
            //  if(icc>0) {
            cpr->cutElem[icc2] = VInt_create();
            cpr->cutCoef[icc2] = VDbl_create();
            cpr->nAllocations++;
            // }
        }
        //  printf("MipC_quick_sort_vec_by_double: start\n");fflush(stdout);
        CutP_quick_sort_vec_by_double ( idValueCut, valueCute, icc2);
        // qsort(valueCutIdxCO, nc, sizeof(IntDblPair), cmp_int_dbl_pair_b);
        //printf("MipC_quick_sort_vec_by_double: end\n");
        if(VERBOSE==2) printf("\n%d precedence cuts were found. ", nc);
        fflush(stdout);
        int *elem = VInt_getPtr(vecIdx);
        double *coef = VDbl_getPtr(vecCo);

      //  printf("nc %d maxcuts %f ", nc, maxcuts);  getchar();
        for(int icc = 0 ; icc < nc && icc < maxcuts; icc++) {

            int ic = idValueCut[icc];

            if( valueCute[ic] <= 0.0002)
            {
                printf("Non Violated cut %f",  valueCute[ic]);
                continue;
            }
       //     printf("nCut %d sizeCut %d", VInt_get(vecNCut,icP));

            int sizeCut = VInt_get(vecNCut,icP++);
            int posCut = VInt_get(vecNCut,icP++);
        //    printf("sizeCut %d", sizeCut);
         //   printf("\n");
            for(int ec = 0 ; ec < sizeCut ; ec++ ) {
                VInt_pushBack(cpr->cutElem[icc],elem[posCut]);
                VDbl_pushBack(cpr->cutCoef[icc],coef[posCut]);
           //        printf("%f*%d ",VDbl_get(cpr->cutCoef[icc],ec),VInt_get(cpr->cutElem[icc],ec));
                posCut++;
            }

            int *cutIdx = VInt_getPtr(cpr->cutElem[icc]);
            double *cutCoef = VDbl_getPtr(cpr->cutCoef[icc]);
            CutP_quick_sort_vec(cutIdx, cutCoef, sizeCut);
//             printf("\nSENSE %c RHS %f VI %f\n", sen, 0.0, valueCute[ic]);fflush(stdout);

            lp_add_cut_grb( cbdata,  sizeCut, cutIdx, cutCoef, GRB_GREATER_EQUAL, 0.0 );
            //   getchar();
            ncut2++;
//            nCut++;
            if(*timeLeft < omp_get_wtime()-tinit +1) break;
        }
        free(elem);
        free(coef);
        free(valueCute);
        free(idValueCut);

    }

//    printf("\nCP added in CALLBACK: %d\n", ncut2); fflush(stdout);//getchar();
    *nCuts += ncut2;
    double _time;
    _time = ( (double) *timeLeft - (omp_get_wtime()-tinit) );
    *timeLeft = _time;

    free(xf);
    if(nc>0){
        free(vecIdx);
        free(vecCo);
    }else{
        VInt_free(&vecIdx);
        VDbl_free(&vecCo);
    }
    VInt_free(&vecNCut);
    VDbl_free(&vecValue);

    free(xIdx[0][0]);
    free(xIdx[0]);
    free(xIdx);
    free(x[0][0]);
    free(x[0]);
    free(x);

    for(int tji = 0 ; tji <nJobs ; tji++) {
        free(cont[tji]);
        for(int m = 0 ; m <nModes ; m++)
           free(TJ[tji][m]);
        free(TJ[tji]);
        free(TJJ[tji]);
    }
    free(TJ);
    free(TJJ);
    free(cont);
    free(contJ);

    free(firstTJ);
    free(lastTJ);

    CutPR_free(&cpr);
    //double endPC = omp_get_wtime()- iniPC;
    //printf("timePC %f \n", endPC); fflush(stdout);


    //  cpr->timeseparation = (double) omp_get_wtime()-tinit;

}

/*free memory*/
void CutPR_free( CutPR **_cutPR )
{

    CutPR *cutPR = *_cutPR;

    for(int nr = 0 ; nr < cutPR->nAllocations; nr++) {
        VInt_free(&cutPR->cutElem[nr]);
        VDbl_free(&cutPR->cutCoef[nr]);
    }

    free(cutPR->cutElem);
    free(cutPR->cutCoef);

    VStr_free(&cutPR->cutname);
    VDbl_free(&cutPR->cutrhs);
    VDbl_free(&cutPR->cutviolation);
    VInt_free(&cutPR->cutnelem);
    VInt_free(&cutPR->cutdominated);
    VInt_free(&cutPR->cutsense);

    free( cutPR );
    *_cutPR = NULL;

}
