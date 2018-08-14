/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * mip_compact: routines for the creation of compact MIP formulations
 *              for the whole problem or for just one project
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */
#include "cut_clique.h"
#include "build_cgraph.h"
#include "oddhs.h"
#ifdef GRB
#include "gurobi_c.h"
#endif
#define VERBOSE  1
#define PERCMAXCLIQUE  0.5

#define MIN_VIOL 0.02


CutCL *CutCL_create(const Instance *inst)
{

    CutCL *ccl;
    ALLOCATE_INI(ccl,CutCL);

    //ALLOCATE_VECTOR(ccl->cutElem,VecInt*,1);
    //ALLOCATE_VECTOR(ccl->cutCoef,VecDbl*,1);
    //ccl->cutElem[0] = VInt_create();
    //ccl->cutCoef[0] = VDbl_create();
    ccl->nAllocations = 0;
    ccl->cutrhs = VDbl_create();
    ccl->cutviolation = VDbl_create();
    ccl->cutnelem = VInt_create();
    ccl->cutsense = VInt_create();
    ccl->cutname = VStr_create(256);
    ccl->cutdominated = VInt_create();
    //  ccl->nAddCuts =0;
    // ccl->timeseparation = 0.0;

    ccl->inst = inst;

    return ccl;
}


void CutCL_free( CutCL **_cutCL )
{

    CutCL *cutCL = *_cutCL;

    //    printf("cutCL->nAllocations %d cutCL->cutsense %d", cutCL->nAllocations, VInt_size(cutCL->cutsense));
    //   fflush(stdout);
    //for(int nc = 0; nc < VInt_size(cutCL->cutsense); nc++) {
    for(int nc = 0; nc < cutCL->nAllocations; nc++) {
        VInt_free(&cutCL->cutElem[nc]);
        VDbl_free(&cutCL->cutCoef[nc]);
    }
    free(cutCL->cutElem);
    free(cutCL->cutCoef);

    VStr_free(&cutCL->cutname);
    VDbl_free(&cutCL->cutrhs);
    VInt_free(&cutCL->cutnelem);
    VInt_free(&cutCL->cutsense);
    VInt_free(&cutCL->cutdominated);
    VDbl_free(&cutCL->cutviolation);

    free( cutCL );
    *_cutCL = NULL;

}

void CutCL_add_cuts_conflicts_clique_parallel(  const CGraph *cgraph, CutCL *ccl, LinearProgram *lp, const Instance *inst, double timeLeft,  double maxcuts,  int nround )
{

    double tinit = omp_get_wtime();

    int ncols = lp_cols(lp);
    int nCut =  lp_rows(lp);

    const double *xf = lp_x(lp);
    double maxRC = lp_get_max_reduced_cost(lp);

    CliqueSeparation *sep = clq_sep_create(cgraph);

    int nRows = lp_rows(lp)*PERCMAXCLIQUE;
    if(nRows<maxcuts) nRows = maxcuts;
    clq_sep_set_maxcliques(sep,  nRows);
    clq_sep_set_max_it_bk(sep, 1000);

    double _time =  ( (double) timeLeft - (omp_get_wtime()-tinit) );
   // printf("before clq_sep_set_extend_method %f. \n", _time);fflush(stdout);
    if(maxRC != -1.0) {
        clq_sep_set_rc(sep,lp_reduced_cost(lp));
    } else
        clq_sep_set_extend_method(sep,0);

    _time = ( (double) timeLeft - (omp_get_wtime()-tinit) );
 //   printf("after clq_sep_set_extend_method %f.\n", _time); fflush(stdout);

    if(_time <= clq_sep_get_max_time_bk(sep))
        clq_sep_set_max_time_bk(sep,_time);
   // clq_sep_params_print(sep);


  //  printf("before clq_sep_separate %f. \n", _time); fflush(stdout);
    clq_sep_separate( sep, xf);
  //  printf("after clq_sep_separate %f. \n", _time);  fflush(stdout);

    const CliqueSet *clqSet = clq_sep_get_cliques( sep );
    int nCliques = clq_get_number_of_cliques(clqSet);

    CliqueEliteSet *clqES;
    int nCliquesElite = 0;
    if(nCliques > 0) {
        if(maxcuts<1) {

            maxcuts =  maxcuts*lp_rows(lp);
            clqES = clqES_create( clqSet, ncols,xf, maxcuts, 0.75, lp_cols(lp) );
        } else
            clqES = clqES_create( clqSet, ncols,xf, maxcuts, 0.75, lp_cols(lp));

        nCliquesElite = clqES_size(clqES);
    }


    int i = 0, j=0;


    if(VERBOSE==2) printf("\n%d clique cuts were found. ", nCliques);            fflush(stdout);
    if(VERBOSE==2) printf("maxcuts %f. ", maxcuts);            fflush(stdout);

    ALLOCATE_VECTOR(ccl->cutElem,VecInt*,maxcuts);
    ALLOCATE_VECTOR(ccl->cutCoef,VecDbl*,maxcuts);

    while(j < nCliquesElite && j < nCliques && j < maxcuts && j < nRows) {

        double vi = 0.0;

        const int* iset = clqES_get_clique_elements(clqES,j);
        assert(iset);
        int size = clqES_get_clique_size(clqES,j);
        if(size==0) {
            fflush(stdout);
            j++;
            continue;
        }
        ccl->cutElem[i] = VInt_create();
        ccl->cutCoef[i] = VDbl_create();
        ccl->nAllocations++;
        for(int sel = 0 ; sel<size ; sel++) {
            VInt_pushBack(ccl->cutElem[i],iset[sel]);
            VDbl_pushBack(ccl->cutCoef[i],1.0);
            vi  += xf[iset[sel]];
        }

        char nameCut[STR_SIZE];
        sprintf( nameCut, "cutClique(%d)#%d#%d", i, nround,lp_rows(lp)+nCut++);

        int *cutIdx = VInt_getPtr(ccl->cutElem[i]);
        double *cutCoef = VDbl_getPtr(ccl->cutCoef[i]);
        VStr_pushBack(ccl->cutname, nameCut);
        VDbl_pushBack(ccl->cutrhs, 1.0);
        VInt_pushBack(ccl->cutsense, 0);
        VInt_pushBack(ccl->cutdominated,0);
        VInt_pushBack(ccl->cutnelem, size);
        VDbl_pushBack(ccl->cutviolation, vi-1);

        CutP_quick_sort_vec(cutIdx,cutCoef, size);

        j++;
        i++;
    }

    /*  double timeinitdomresource = omp_get_wtime();
      int dom = 0;
      for(int ncdA = 0; ncdA < i ; ncdA++){
          for(int ncdB = ncdA+1; ncdB < i ; ncdB++){
              if(VInt_get(ccl->cutdominated,ncdA) == 1) break;
              if(VInt_get(ccl->cutdominated,ncdB) == 1) continue;
              CutP_dominanceBetweenTwo( lp, ncdA, VInt_getPtr(ccl->cutElem[ncdA]), VDbl_getPtr(ccl->cutCoef[ncdA]), VDbl_get(ccl->cutrhs,ncdA), VInt_get(ccl->cutsense,ncdA), VInt_get(ccl->cutnelem,ncdA), ncdB, VInt_getPtr(ccl->cutElem[ncdB]), VDbl_getPtr(ccl->cutCoef[ncdB]), VDbl_get(ccl->cutrhs,ncdB), VInt_get(ccl->cutsense,ncdB), VInt_get(ccl->cutnelem,ncdB),  ccl->cutdominated);
           //   getchar();
          }
      }
      printf("time compute dominance clique %f\n",omp_get_wtime()-timeinitdomresource);
    */

    clq_sep_free(&sep);
    if(nCliques>0)
        clqES_free(&clqES);

#undef MAX_IDX
}

void CutCL_add_cuts_clique_callback(  GRBmodel *model,   void *cbdata, const Instance *inst, double *timeLeft, int lifting, double maxcuts, int  *nCuts)
{
    double tinit = omp_get_wtime();
    //double iniPC = omp_get_wtime();
    double maxRC =0.0;
    CutCL *ccl = CutCL_create(inst);
//    int ncc = *nCuts;

    int numvars=0;
    int numrows=0;


    int error = GRBgetintattr(model, GRB_INT_ATTR_NUMVARS, &numvars);
  //  printf("numbers of variables callback: %d\n", numvars);// getchar();

    error = GRBgetintattr(model, GRB_INT_ATTR_NUMCONSTRS, &numrows);
    double *xf;
    ALLOCATE_VECTOR_INI(xf,double,numvars);
    double *rdc;
    ALLOCATE_VECTOR_INI(rdc,double,numvars);
    char **varname;

    GRBcbget(cbdata,GRB_CB_MIPNODE,GRB_CB_MIPNODE_REL,xf);

    varname = malloc(numvars*sizeof(char*));

    for (int j = 0; j < numvars; ++j ) {
        error = GRBgetstrattrelement(model, GRB_STR_ATTR_VARNAME, j, &varname[j]);
        error = GRBgetdblattrelement(model, GRB_DBL_ATTR_RC, j, &rdc[j]);

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


    for ( i=0 ; (i<numvars) ; ++i ) {
        parseName( varname[i], prefix, idx );
       // printf("%s %f\n", varname[i], xf[i]);
        if (prefix[0]=='x') {
  //          int j = idx[0];
            const int t = idx[2];
//            const int m = idx[1];
            maxT = MAX( maxT, t );
        }
    }

    int nTimes = maxT+1;




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


    int ncols = numvars;

    int **conf;
    ALLOCATE_VECTOR(conf, int*, ncols);

    int *nConf;
    ALLOCATE_VECTOR_INI(nConf, int, ncols);

    int **varInConfs;
    ALLOCATE_VECTOR_INI(varInConfs, int*, ncols);


    int nProj = Inst_nProjects(inst);
    int nResources = Inst_nResR(inst);

    int **vectorJobs;
    ALLOCATE_VECTOR(vectorJobs,int*,nJobs);
    for(int i = 0 ; i < nJobs ; i++)
        ALLOCATE_VECTOR_INI(vectorJobs[i],int,ncols);

    IntTriple **vectorProjects;
    ALLOCATE_VECTOR(vectorProjects,IntTriple*,nProj);
    for(int i = 0 ; i < nProj ; i++)
        ALLOCATE_VECTOR_INI(vectorProjects[i],IntTriple,ncols);

    IntTriple ***vectorTimes;
    ALLOCATE_VECTOR(vectorTimes,IntTriple**,nTimes);
    for(int t = 0 ; t < nTimes ; t++)
    {
        ALLOCATE_VECTOR(vectorTimes[t],IntTriple*,nResources);
        for(int r = 0 ; r < nResources ; r++)
            ALLOCATE_VECTOR(vectorTimes[t][r],IntTriple,ncols);
    }

    int *nVJobs;
    ALLOCATE_VECTOR_INI(nVJobs,int,nJobs);
    int *nVProj;
    ALLOCATE_VECTOR_INI(nVProj,int,nProj);
    int **nVTime;
    ALLOCATE_VECTOR(nVTime,int*,nTimes);
    for(int t = 0 ; t < nTimes; t++)
        ALLOCATE_VECTOR_INI(nVTime[t],int,nResources);

    VecInt** conflicts;
    ALLOCATE_VECTOR(conflicts,VecInt*,ncols);

    double _time;
    for(int l = 0; l < ncols ; l++)
    {
        ALLOCATE_VECTOR_INI(varInConfs[l], int, ncols);
        conflicts[l] = VInt_create();
        if (fabs(xf[l])>1e-5 || maxRC != -1.0)
        {
            if(maxRC != -1.0)
                if(rdc[l] > maxRC) continue;

            if (tolower(varname[l][0])=='x')
            {
                parseName( varname[l], prefix, idx );
                if (prefix[0]=='x')
                {

                    int jlp = idx[0];
                    int mlp =idx[1];
                    int tlp = idx[2];

                    const Job *job = Inst_job(inst,jlp);
                    const Mode *mode = Job_mode(job,mlp);
                    if(Mode_duration(mode)==0)continue;

                    int plp = Job_project(job);

                    if(xIdx[jlp][mlp][tlp] == -1) continue;

                    vectorJobs[jlp][nVJobs[jlp]] = l;

                    for(int r = 0; r < nResources ; r++)
                    {

                        int idxROnInstD = Mode_idxResROnMode(inst,mode,r);
                        if( idxROnInstD == -1) continue;
                        int size = nVTime[tlp][r];
                        vectorTimes[tlp][r][size].idx = l;
                        vectorTimes[tlp][r][size].j = jlp;
                        vectorTimes[tlp][r][size].m = mlp;
                        vectorTimes[tlp][r][size].t = tlp;
                        vectorTimes[tlp][r][size].value = xf[l];

                        nVTime[tlp][r]++;

                    }

                    vectorProjects[plp][nVProj[plp]].idx = l;
                    vectorProjects[plp][nVProj[plp]].j = jlp;
                    vectorProjects[plp][nVProj[plp]].m = mlp;
                    vectorProjects[plp][nVProj[plp]].t = tlp;
                    vectorProjects[plp][nVProj[plp]].value = xf[l];

                    //   printf("(%d, %d, %d) %d | ", jlp,mlp,tlp,vectorJobs[jlp][nVJobs[jlp]]  );
                    nVJobs[jlp]++;
                    nVProj[plp]++;

                }
            }
        }
    }

    free(varname);

    for(int t = 0 ; t < nTimes ; t++)
    {
        for(int r = 0 ; r < nResources ; r++)
        {
            if(nVTime[t][r] == 0) continue;
            for(int elem = 0 ; elem < nVTime[t][r] ; elem++)
            {
                int j = vectorTimes[t][r][elem].j, m = vectorTimes[t][r][elem].m;// t1 = vectorTimes[t][r][elem].t;
                int l = vectorTimes[t][r][elem].idx;
                const Job *job = Inst_job(inst,j);
                const Mode *mode = Job_mode(job,m);
                if(Mode_duration(mode)==0)continue;
                //   printf("\nE(%d,%d,%d)", j, m, t1);
                for(int elem2 = 0 ; elem2 < nVTime[t][r] ; elem2++)
                {
                    if(elem == elem2 ) continue;
                    int j2 = vectorTimes[t][r][elem2].j, m2 = vectorTimes[t][r][elem2].m,  t2 = vectorTimes[t][r][elem2].t;
                    if(j==j2 && m==m2 && t ==t2) continue;
                    const Job *job2 = Inst_job(inst,j2);
                    const Mode *mode2 = Job_mode(job2,m2);
                    if(Mode_duration(mode2)==0)continue;

                    int res = Mode_idxResROnMode(inst,mode,r);
                    int res2 = Mode_idxResROnMode(inst,mode2,r);
                    int sumUse = Mode_useResR(mode,res)+ Mode_useResR(mode2,res2);
                    if( sumUse > Inst_capResR(inst,r))
                    {
                        int l2 = vectorTimes[t][r][elem2].idx;
                        if(varInConfs[l][l2]==1) continue;
                        //         printf("R (%d,%d,%d) %d ", j2, m2, t2, l2);
                        VInt_pushBack(conflicts[l],l2);
                        varInConfs[l][l2] = 1;
                    }
                }
            }
        }
    }

    for(int p = 0 ; p < nProj ; p++)
    {
        for(int elem = 0 ; elem < nVProj[p] ; elem++)
        {
            int l = vectorProjects[p][elem].idx;
            int j = vectorProjects[p][elem].j, m = vectorProjects[p][elem].m, t =  vectorProjects[p][elem].t;
            const Job *job = Inst_job(inst,j);
            const Mode *mode = Job_mode(job,m);
            if(Mode_duration(mode)==0)continue;
            //    printf("\nE(%d,%d,%d)", j, m, t);
            for(int elem2 = 0 ; elem2 < nVProj[p] ; elem2++)
            {
                if(elem == elem2 ) continue;
                int j2 = vectorProjects[p][elem2].j, m2 = vectorProjects[p][elem2].m, t2 = vectorProjects[p][elem2].t;
                if(j==j2 && m==m2 && t ==t2) continue;

                const Job *job2 = Inst_job(inst,j2);
                const Mode *mode2 = Job_mode(job2,m2);
                if(Mode_duration(mode2)==0)continue;

                int durationE = Mode_duration(mode);
                int timeEndPred = t+durationE;

                int over = 0;
                if( Job_hasIndSucc(job,j2) )
                {
                    int value = Inst_getMaxDIJ(inst,j,j2)- Job_minDuration(job);
                    int winTime = value ;
                    if(timeEndPred>t2 || ((t2-timeEndPred) < winTime))
                    {
                        //     printf("timeEndPred = t %d + durationE %d \n", t, durationE );
                        //   printf("\ntimeEndPred %d > t2 %d || ((t2 %d - timeEndPred %d) < winTime %d)", timeEndPred,t2,t2,timeEndPred, winTime);
                        int l2 = vectorProjects[p][elem2].idx;
                        if(varInConfs[l][l2]==1) continue;
                        //    printf("timed D(%d,%d,%d) %d ", j2, m2, t2, l2);
                        VInt_pushBack(conflicts[l],l2);
                        varInConfs[l][l2] = 1;
                        continue;
                    }

                    int part1 = 0, cp1 = 0;
                    const Project *project = Inst_project(inst,p);
                    cp1 = Project_releaseDate(project)+Project_criticalPath(project) -  Inst_getMaxDIMJM(inst, job,Mode_index(mode),j2,Mode_index(mode2));
                    part1 = t - cp1 < 0 ? 0 :
                            t - cp1 ;
                    over =  part1;
                    if(over > Inst_getSumTPD(inst))
                    {
                        int l2 = vectorProjects[p][elem2].idx;
                        if(varInConfs[l][l2]==1) continue;
                        //   printf("tpd(%d,%d,%d) %d ", j2, m2, t2, l2);
                        VInt_pushBack(conflicts[l],l2);
                        varInConfs[l][l2] = 1;
                    }
                }
            }

            for(int p2 = 0 ; p2 < nProj ; p2++)
            {
                if(p==p2) continue;
                for(int elem2 = 0 ; elem2 < nVProj[p2] ; elem2++)
                {
                    int j2 = vectorProjects[p2][elem2].j, m2 = vectorProjects[p2][elem2].m, t2 = vectorProjects[p2][elem2].t;
                    if(j==j2 && m==m2 && t ==t2) continue;

                    const Job *job2 = Inst_job(inst,j2);
                    const Mode *mode2 = Job_mode(job2,m2);
                    if(Mode_duration(mode2)==0)continue;

                    const Project *project = Inst_project(inst,p);
                    int cp1 = Project_releaseDate(project)+Project_criticalPath(project) - Inst_getMaxDIM(inst,job,Mode_index(mode));
                    const Project *project2 = Inst_project(inst,Job_project(job2));
                    int cp2 = Project_releaseDate(project2)+Project_criticalPath(project2) - Inst_getMaxDIM(inst,job2,Mode_index(mode2));

                    int part1 = t - cp1 < 0 ? 0 :
                                t - cp1 ;
                    int part2 = t2 - cp2 < 0 ? 0 :
                                t2 - cp2 ;
                    if(part1 + part2 > Inst_getSumTPD(inst))
                    {
                        // printf("\npart1 %d (t%d-cp1 %d) + part2 %d (t2%d-cp2 %d) > Inst_getSumTPD(inst) %d  releaseDatep1 %d + criticalPathp1 %d  DJMj1 %d releaseDatep2 %d + criticalPathp %d DIM2 %d", part1, t, cp1, part2, t2, cp2, Inst_getSumTPD(inst), Project_releaseDate(project),Project_criticalPath(project) , Inst_getMaxDIM(inst,job,Mode_index(mode)), Project_releaseDate(project2),Project_criticalPath(project2),Inst_getMaxDIM(inst,job2,Mode_index(mode2)));
                        int l2 = vectorProjects[p2][elem2].idx;
                        if(varInConfs[l][l2]==1) continue;
                        //  printf("tpd p(%d,%d,%d) %d ", j2, m2, t2, l2);
                        VInt_pushBack(conflicts[l],l2); //if the sum of ends is bigger then tpd
                        varInConfs[l][l2] = 1;
                    }
                }
            }
        }
    }
    //    int weight[lp_cols(lp)];

    for(int l = 0 ; l < ncols ; l++)
    {
        int sizeConflictsJ = VInt_size(conflicts[l]);
        if(sizeConflictsJ != 0)
        {
            nConf[l] = sizeConflictsJ;
            int *confJ = VInt_getPtr(conflicts[l]);
            ALLOCATE_VECTOR(conf[l], int, sizeConflictsJ);
            for(int o = 0 ; o <  sizeConflictsJ; o++)
                conf[l][o]= confJ[o];
            // if(xf[l]<=0.00001)
            //   weight[l] = 1;
            //else
            //  weight[l] = xf[l]*1000.0;

        }//else
        //weight[l] = 0.0;
    }




    /* for(int c = 0; c < lp_cols(lp); c++) {
         if(nConf[c] > 0) {
             Res_setNSumAllVarWithConf(res,LPC_CLIQUE,1);
             Res_setNSumAllElementsConf(res,LPC_CLIQUE,nConf[c]);
             if(nConf[c] > Res_getNMaxElementsConf(res,LPC_CLIQUE))
                 Res_setNMaxElementsConf(res,LPC_CLIQUE,nConf[c]);
         }
     }
     */

    /*  printf("Conflicts\n");
       for(int ci = 0 ; ci < ncols ;ci++){
          printf("\nE %d\n", ci);
          int size = VInt_size(conflicts[ci]);
          printf("D: ");
          for(int cie = 0 ; cie < size ; cie++){
              printf(" %d ", VInt_get(conflicts[ci], cie));
          }
       }*/ // getchar();


    _time = ( (double)  (omp_get_wtime()-tinit) );
    if(VERBOSE==3) printf("\ntime finding conflicts: %f. ", _time);

    _time = ( (double) *timeLeft - (double) (omp_get_wtime()-tinit) );

    CGraph *cgraph = build_cgraph_conflicts(conf, nConf, ncols,_time);
    char recomputeDegree = 0;
    for(int j = 0 ; j < nJobs ; j++)
    {
        if(nVJobs[j]==0) continue;
        cgraph_add_clique(cgraph,vectorJobs[j],nVJobs[j]);
        recomputeDegree = 1;
    }

    if(recomputeDegree)
        cgraph_recompute_degree(cgraph);

    for(int t = 0 ; t < nTimes; t++)
    {
        free(nVTime[t]);
        for(int r = 0 ; r < nResources; r++)
            free(vectorTimes[t][r]);
        free(vectorTimes[t]);
    }
    for(int j = 0 ; j < nJobs ; j++)
        free(vectorJobs[j]);

    for(int p = 0 ; p < nProj ; p++)
        free(vectorProjects[p]);
    free(vectorJobs);
    free(vectorProjects);
    free(vectorTimes);
    free(nVJobs);
    free(nVProj);
    free(nVTime);

    for(int l = 0 ; l < ncols ; l++)
    {
        free(varInConfs[l]);
        VInt_free(&conflicts[l]);
        if(nConf[l] == 0) continue;
        free(conf[l]);
    }

    free(conflicts);
    free(varInConfs);
    free(conf);
    free(nConf);

//    cgraph_preprocess(cgraph, )
    int nCut=0;
    CliqueSeparation *sep = clq_sep_create(cgraph);

    int nRows = numrows*PERCMAXCLIQUE;
    if(nRows<maxcuts) nRows = maxcuts;
    clq_sep_set_maxcliques(sep,  nRows);
    clq_sep_set_max_it_bk(sep, 1000);

    _time =  ( (double) *timeLeft - (omp_get_wtime()-tinit) );
   // printf("before clq_sep_set_extend_method %f. \n", _time);fflush(stdout);
    clq_sep_set_extend_method(sep,0);

    _time = ( (double) *timeLeft - (omp_get_wtime()-tinit) );
 //   printf("after clq_sep_set_extend_method %f.\n", _time); fflush(stdout);

    if(_time <= clq_sep_get_max_time_bk(sep))
        clq_sep_set_max_time_bk(sep,_time);
   // clq_sep_params_print(sep);


  //  printf("before clq_sep_separate %f. \n", _time); fflush(stdout);
    clq_sep_separate( sep, xf);
  //  printf("after clq_sep_separate %f. \n", _time);  fflush(stdout);

    const CliqueSet *clqSet = clq_sep_get_cliques( sep );
    int nCliques = clq_get_number_of_cliques(clqSet);

    CliqueEliteSet *clqES;
    int nCliquesElite = 0;
    if(nCliques > 0) {
        if(maxcuts<1) {

            maxcuts =  maxcuts*numrows;
            clqES = clqES_create( clqSet, ncols,xf, maxcuts, 0.75, ncols );
        } else
            clqES = clqES_create( clqSet, ncols,xf, maxcuts, 0.75,ncols);
        nCliquesElite = clqES_size(clqES);
    }


    int j=0;
    i = 0;


    if(VERBOSE==2) printf("maxcuts %f. ", maxcuts);            fflush(stdout);

    ALLOCATE_VECTOR(ccl->cutElem,VecInt*,maxcuts);
    ALLOCATE_VECTOR(ccl->cutCoef,VecDbl*,maxcuts);

    while(j < nCliquesElite && j < nCliques && j < maxcuts && j < nRows) {

        double vi = 0.0;

        const int* iset = clqES_get_clique_elements(clqES,j);
        assert(iset);
        int size = clqES_get_clique_size(clqES,j);
        if(size==0) {
            fflush(stdout);
            j++;
            continue;
        }
        ccl->cutElem[i] = VInt_create();
        ccl->cutCoef[i] = VDbl_create();
        ccl->nAllocations++;
        for(int sel = 0 ; sel<size ; sel++) {
            VInt_pushBack(ccl->cutElem[i],iset[sel]);
            VDbl_pushBack(ccl->cutCoef[i],1.0);
            vi  += xf[iset[sel]];
        }

        char nameCut[STR_SIZE];
        sprintf( nameCut, "cutClique(%d)#%d", i,numrows+nCut++);

        int *cutIdx = VInt_getPtr(ccl->cutElem[i]);
        double *cutCoef = VDbl_getPtr(ccl->cutCoef[i]);
        VStr_pushBack(ccl->cutname, nameCut);
        VDbl_pushBack(ccl->cutrhs, 1.0);
        VInt_pushBack(ccl->cutsense, 0);
        VInt_pushBack(ccl->cutdominated,0);
        VInt_pushBack(ccl->cutnelem, size);
        VDbl_pushBack(ccl->cutviolation, vi-1);

        CutP_quick_sort_vec(cutIdx,cutCoef, size);
        lp_add_cut_grb( cbdata,  size, cutIdx, cutCoef, GRB_LESS_EQUAL, 1.0 );

        j++;
        i++;
    }

    /*  double timeinitdomresource = omp_get_wtime();
      int dom = 0;
      for(int ncdA = 0; ncdA < i ; ncdA++){
          for(int ncdB = ncdA+1; ncdB < i ; ncdB++){
              if(VInt_get(ccl->cutdominated,ncdA) == 1) break;
              if(VInt_get(ccl->cutdominated,ncdB) == 1) continue;
              CutP_dominanceBetweenTwo( lp, ncdA, VInt_getPtr(ccl->cutElem[ncdA]), VDbl_getPtr(ccl->cutCoef[ncdA]), VDbl_get(ccl->cutrhs,ncdA), VInt_get(ccl->cutsense,ncdA), VInt_get(ccl->cutnelem,ncdA), ncdB, VInt_getPtr(ccl->cutElem[ncdB]), VDbl_getPtr(ccl->cutCoef[ncdB]), VDbl_get(ccl->cutrhs,ncdB), VInt_get(ccl->cutsense,ncdB), VInt_get(ccl->cutnelem,ncdB),  ccl->cutdominated);
           //   getchar();
          }
      }
      printf("time compute dominance clique %f\n",omp_get_wtime()-timeinitdomresource);
    */

    clq_sep_free(&sep);
    if(nCliques>0)
        clqES_free(&clqES);


    _time = ( (double) *timeLeft - (omp_get_wtime()-tinit) );
    *timeLeft =_time;
    printf("\nCC added in CALLBACK: %d\n", nCut); fflush(stdout);//getchar();
    *nCuts += nCut;

    free(xf);
    free(rdc);
    CutCL_free(&ccl);
    cgraph_free( &cgraph );

    free(xIdx[0][0]);
    free(xIdx[0]);
    free(xIdx);
    free(x[0][0]);
    free(x[0]);
    free(x);

    //double endPC = om
    //double endPC = omp_get_wtime()- iniPC;
    //printf("timePC %f \n", endPC); fflush(stdout);


    //  cpr->timeseparation = (double) omp_get_wtime()-tinit;

}


