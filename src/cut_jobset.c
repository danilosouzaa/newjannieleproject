/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * mip_compact: routines for the creation of compact MIP formulations
 *              for the whole problem or for just one project
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */
#include "cut_jobset.h"
#include "build_cgraph.h"
#include "bron_kerbosch.h"
#include "lp.h"
#define VERBOSE  1


CutJS *CutJS_create(const Instance *inst)
{

    CutJS *cjs;
    ALLOCATE_INI(cjs,CutJS);

    cjs->nAllocations = 0;
    cjs->cutrhs = VDbl_create();
    cjs->cutviolation = VDbl_create();
    cjs->cutnelem = VInt_create();
    cjs->cutsense = VInt_create();
    cjs->cutname = VStr_create(256);
    cjs->cutdominated = VInt_create();

    cjs->inst = inst;

    return cjs;
}


void CutJS_free( CutJS **_cutJS )
{

    CutJS *cutJS = *_cutJS;

    for(int nr = 0 ; nr < cutJS->nAllocations; nr++) {
        VInt_free(&cutJS->cutElem[nr]);
        VDbl_free(&cutJS->cutCoef[nr]);
    }

    free(cutJS->cutElem);
    free(cutJS->cutCoef);

    VStr_free(&cutJS->cutname);
    VDbl_free(&cutJS->cutrhs);
    VDbl_free(&cutJS->cutviolation);
    VInt_free(&cutJS->cutnelem);
    VInt_free(&cutJS->cutdominated);
    VInt_free(&cutJS->cutsense);

    free( cutJS );
    *_cutJS = NULL;

}

void CutJS_add_cuts_jobset_parallel( CutJS *ccjs, LinearProgram *lp, const int *origCols, LinearProgram *origLP, const Instance *inst, double timeLeft, int continuous, int lifting, double maxcuts, int nround )
{

 //    double tinit = omp_get_wtime();
    assert(inst);
//    int nCut = lp_rows(lp);

    const double *xf = lp_x(lp);
    char name[256];
    int i;

#define MAX_IDX 16
    char prefix[STR_SIZE];
    int idx[MAX_IDX];
    int maxT = 0;



    for ( i=0 ; (i<lp_cols(lp)) ; ++i ) {
        lp_col_name( origLP, origCols[i], name );
        parseName( name, prefix, idx );
        if (prefix[0]=='x') {
//            int j = idx[0];
            const int t = idx[2];
//            const int m = idx[1];
            maxT = MAX( maxT, t );
        }
    }

    int nTimes = maxT+1;


    double _time;
    double tinit = omp_get_wtime();

    LinearProgram *mipJobSet = lp_create();
    lp_set_print_messages(mipJobSet,1);
    lp_set_max_nodes(mipJobSet,100);
    lp_set_mip_emphasis( mipJobSet, LP_ME_FEASIBILITY );

    VecStr *namesY = VStr_create( STR_SIZE );
    VecStr *namesZ = VStr_create( STR_SIZE );
    VecStr *namesX = VStr_create( STR_SIZE );
    VecStr *namesW = VStr_create( STR_SIZE );

    int n=5;
    int m=10;
    char namer[STR_SIZE];

    //enphase heuristica e limite de 100 nós para rodar
    double valueobjy[Inst_nJobs(inst)];
    double valueobjz[nTimes];
    double valueobjx[Inst_nJobs(inst)*nTimes];
    double valueobjw[Inst_nJobs(inst)*Inst_nJobs(inst)];

    for(int t = 0 ; t < nTimes ; t++){
        sprintf(name, "z(%d)", t);
        VStr_pushBack( namesZ, name );
        valueobjz[t]=0;
    }
    lp_add_bin_cols( mipJobSet, VStr_size(namesZ),valueobjz, VStr_ptr(namesZ) );


    double valuetj[ Inst_nJobs(inst)][nTimes];

     for(int j = 0 ; j < Inst_nJobs(inst) ; j++) {
        const Job *job = Inst_job(inst,j);
        for(int t = 0 ; t < nTimes ; t++){
            valuetj[j][t] = 0;
            for(int m = 0 ; m < Job_nModes(job) ; m++){
                const Mode *mode = Job_mode(job,m);

                double totalvalues =0;
                for(int tm = t-Mode_duration(mode)+1 ; tm < t ; tm++){
                    double value = 0;
                    char vnamex[STR_SIZE];
                    sprintf( vnamex, "x(%d,%d,%d)", j, m,tm);
                    int elem = lp_col_index(lp, vnamex);
                    if(elem!=-1){
                        value = xf[elem];
                    }else{
                        continue;
                    }

                    int user =0;
                    for(int r = 0 ; r < Mode_nResR(mode) ; r++){
                        user += Mode_useResR(mode,r);
                    }
                    totalvalues += (double)value*(double)user;
                }
                valuetj[j][t] +=  totalvalues;
            }
        }
    }

    int idxw =0,idxx=0;
    for(int j = 0 ; j < Inst_nJobs(inst) ; j++) {
        const Job *job = Inst_job(inst,j);
        char name[STR_SIZE];
        sprintf(name, "y(%d)", j);
        VStr_pushBack( namesY, name );
        valueobjy[j]=0;
        for(int t = 0 ; t < nTimes ; t++){
            char name[STR_SIZE];
            sprintf(name, "x(%d,%d)",j,t);
            VStr_pushBack( namesX, name );
            valueobjx[idxx] = -1*valuetj[j][t];
            //printf("%s, %f\n",name,  valueobjx[idxx]); fflush(stdout);
            idxx++;
        }
        for(int j2 = 0 ; j2 < Job_nSucc(job); j2++) {
            int idj2 = Job_succ(job,j2);
            char name[STR_SIZE];
            sprintf(name, "w(%d,%d)", j,idj2);
            VStr_pushBack( namesW, name );
            valueobjw[idxw]=-1*EPS;
            idxw++;

        }
    }
    lp_add_bin_cols( mipJobSet, VStr_size(namesY),valueobjy, VStr_ptr(namesY) );
    lp_add_bin_cols( mipJobSet, VStr_size(namesW),valueobjw, VStr_ptr(namesW) );
    lp_add_bin_cols( mipJobSet, VStr_size(namesX),valueobjx, VStr_ptr(namesX) );



    for(int t = 0 ; t < nTimes ; t++) {
        for(int t2 = 0 ; t2 < nTimes ; t2++) {
            if((t2-t)<=m)continue;
            sprintf( name, "z(%d)", t);
            int idxZ1 = lp_col_index(mipJobSet, name);
            sprintf( name, "z(%d)", t2);
            int idxZ2 = lp_col_index(mipJobSet, name);

            int idxZ[2];
            idxZ[0] = idxZ1;
            idxZ[1] = idxZ2;
            double coef[2];
            coef[0] = 1.0;
            coef[1] = 1.0;

            char namer[STR_SIZE];
            sprintf(namer, "timez(%d,%d)", t,t2);
            lp_add_row( mipJobSet, 2, idxZ, coef, namer, 'L', 1.0 );
        }
    }


    int idxY[Inst_nJobs(inst)];
    double coefY[Inst_nJobs(inst)];
    int idxYX[2];
    double coefYX[2];
    int idxYZ[2];
    double coefYZ[2];

    int idxWY[2];
    double coefWY[2];

    int idxWYJ2[2];
    double coefWYJ2[2];

    for(int j = 0 ; j < Inst_nJobs(inst) ; j++) {

        char name[STR_SIZE];
        sprintf( name, "y(%d)", j);
        int id = lp_col_index(mipJobSet, name);
        idxY[j] = id;
        coefY[j] = 1.0;
        const Job *job = Inst_job(inst,j);

        for(int t = 0 ; t < nTimes ; t++) {
            char name[STR_SIZE];
            sprintf( name, "x(%d,%d)", j,t);
            int idx = lp_col_index(mipJobSet, name);
            idxYX[0] =  idx;
            idxYX[1] =  id;
            coefYX[0] = 1.0;
            coefYX[1] = -1.0;

            sprintf(namer, "constxy(%d,%d)",j,t);
            lp_add_row( mipJobSet, 2, idxYX, coefYX, namer, 'L', 0 );

            sprintf( name, "z(%d)", t);
            int idz = lp_col_index(mipJobSet, name);
            idxYZ[0] =  idx;
            idxYZ[1] =  idz;
            coefYZ[0] = 1.0;
            coefYZ[1] = -1.0;
            sprintf(namer, "constxz(%d,%d)",j,t);
            lp_add_row( mipJobSet, 2, idxYZ, coefYZ, namer, 'L', 0 );
        }


        for(int j2 = 0 ; j2 < Job_nSucc(job) ; j2++) {
            int idj2 = Job_succ(job,j2);
            sprintf( name, "w(%d,%d)", j,idj2);
            int idw = lp_col_index(mipJobSet, name);
            idxWY[0] = idw;
            coefWY[0] = 1.0;
            idxWY[1] = id;
            coefWY[1] = -1.0;
            sprintf(namer, "constwyj1(%d,%d,%d)",j,idj2,j);
            lp_add_row( mipJobSet, 2, idxWY, coefWY, namer, 'L', 0 );

            sprintf( name, "y(%d)", idj2);
            int idyj2 = lp_col_index(mipJobSet, name);
            idxWYJ2[0] = idw;
            coefWYJ2[0] = 1.0;
            idxWYJ2[1] = idyj2;
            coefWYJ2[1] = -1.0;
            sprintf(namer, "constwyj2(%d,%d,%d)",j,idj2,idj2);
            lp_add_row( mipJobSet, 2, idxWYJ2, coefWYJ2, namer, 'L', 0 );

        }
    }
    sprintf(namer, "maxnj");
    lp_add_row( mipJobSet, Inst_nJobs(inst), idxY, coefY, namer, 'L', n );



    printf("printing model") ;
    lp_write_lp(mipJobSet, "testeFistLP.lp");
    //lp_set_concurrentMIP(mipJobSet,1);
    //lp_set_method(mipJobSet,4);
    //lp_set_seed(mipJobSet,100000);
    lp_optimize(mipJobSet);
    lp_write_sol(mipJobSet, "testeFirstSol.sol");
  //  getchar();

    int id[4];
    double *xfl = lp_x(mipJobSet);
    int jobSet[Inst_nJobs(inst)];
    int jobSetTime[Inst_nJobs(inst)];
    int cont=0;
    int selmaxt  = 0;
    int selmint = INT_MAX;
    for(int l = 0 ; l < lp_cols(mipJobSet); l++){
        if (fabs(xfl[l])>1e-5 ) {
            char vnamey[STR_SIZE];
            lp_col_name(mipJobSet,l,vnamey);
            if (tolower(vnamey[0])=='x') {
                parseName( vnamey, prefix, id );
                if (prefix[0]=='x') {
                   jobSet[cont] = id[0];
                   jobSetTime[cont] = id[1];
                   cont++;
                   if(id[1]>selmaxt){
                        selmaxt = id[1];
                   }
                   if(id[1]<selmint){
                        selmint = id[1];
                   }
                 }
            }
        }
    }

    lp_free(&mipJobSet);


    double ***coefR;
    int **contelemcoefrt;
    ALLOCATE_VECTOR_INI(coefR,double**,Inst_nResR(inst));
    ALLOCATE_VECTOR_INI(contelemcoefrt,int*,Inst_nResR(inst));
    for(int res = 0 ; res < Inst_nResR(inst); res++){
        ALLOCATE_VECTOR_INI(coefR[res],double*,nTimes);
        ALLOCATE_VECTOR_INI(contelemcoefrt[res],int,nTimes);
        for(int td = 0; td < nTimes ;td++){
            ALLOCATE_VECTOR_INI(coefR[res][td],double,cont*Inst_nMaxModes(inst));
        }
    }

    double **coefNR;
    int *contelemcoefRN;
    ALLOCATE_VECTOR_INI(coefNR,double*,Inst_nResN(inst));
    ALLOCATE_VECTOR_INI(contelemcoefRN,int,Inst_nResN(inst));
    for(int res = 0 ; res < Inst_nResN(inst); res++){
        ALLOCATE_VECTOR_INI(coefNR[res],double,cont*Inst_nMaxModes(inst)*(selmaxt-selmint));
    }


    LinearProgram *mipCut = lp_create();
    lp_set_print_messages(mipCut,1);


    VecStr *namesVars = VStr_create( STR_SIZE );
    VecInt *idxVars = VInt_create();
//    VecDbl *coefVars = VDbl_create();
    VecDbl *objVars = VDbl_create();

    VecInt *cutsElem= VInt_create();
    VecDbl *coefElem= VDbl_create();

    // maximizar os x * pelo tempo (coeficiente do corte) <= rhs
    //por tempo os recursos

    int cutnelem = 0;
    for(int j = 0; j < cont ; j++ ){
        const Job* job = Inst_job(inst,jobSet[j]);
        for(int m = 0 ; m < Job_nModes(job) ; m++){
            const Mode *mode = Job_mode(job,m);
            //for(int i = s ; i < t ; i++){
                char vnamex[STR_SIZE];
                sprintf( vnamex, "x(%d,%d,%d)", jobSet[j], m,jobSetTime[j]);
               // printf( "x(%d,%d,%d)", jobSet[j], m,i); fflush(stdout);
                int elemlp = lp_col_index(lp, vnamex);
                if(elemlp==-1) continue;
                VInt_pushBack(cutsElem,elemlp);
                VDbl_pushBack(coefElem,jobSetTime[j]*xf[elemlp]);

                VStr_pushBack(namesVars,vnamex);
                VInt_pushBack(idxVars,cutnelem);
                VDbl_pushBack(objVars,-1*(jobSetTime[j]*xf[elemlp]));
                for(int res = 0 ; res < Mode_nResR(mode); res++){
                    for(int ttd = jobSetTime[j] ; ttd< jobSetTime[j]+Mode_duration(mode)+1 ; ttd++){
                        int r = Mode_idxResR(mode,res);
                        coefR[r][ttd][cutnelem] = Mode_useResR(mode,res);
                        contelemcoefrt[r][ttd]++;
                    }
                }

                for(int res = 0 ; res < Mode_nResN(mode); res++){
                    int nr = Mode_idxResN(mode,res);
                    coefNR[nr][cutnelem] = Mode_useResN(mode,res);
                    contelemcoefRN[nr]++;
                }
                cutnelem++;
            //}
        }
    }

    lp_add_bin_cols( mipCut, VStr_size(namesVars), VDbl_getPtr(objVars), VStr_ptr(namesVars) );

    for(int res = 0 ; res < Inst_nResR(inst); res++){
        for(int td = selmint; td <= selmaxt ;td++){
            if(contelemcoefrt[res][td]==0) continue;
            char nameR[STR_SIZE];
            sprintf( nameR, "cResR(%d,%d)", res,td);
            lp_add_row(mipCut,cutnelem,VInt_getPtr(idxVars),coefR[res][td],nameR,'L',Inst_capResR(inst,res));
        }
    }
    for(int res = 0 ; res < Inst_nResN(inst); res++){
        if(contelemcoefRN[res]==0) continue;
       char nameR[STR_SIZE];
       sprintf( nameR, "cResNR(%d)", res);
       lp_add_row(mipCut,cutnelem,VInt_getPtr(idxVars),coefNR[res],nameR,'L',Inst_capResN(inst,res));
    }

    VecInt** conflicts;
    ALLOCATE_VECTOR(conflicts,VecInt*,cutnelem);

    int **varInConfs;
    ALLOCATE_VECTOR_INI(varInConfs, int*, cutnelem);
    int nProj = Inst_nProjects(inst);
    int *nVProj;
    ALLOCATE_VECTOR_INI(nVProj,int,nProj);
    IntTriple **vectorProjects;
    ALLOCATE_VECTOR(vectorProjects,IntTriple*,nProj);
    for(int i = 0 ; i < nProj ; i++)
        ALLOCATE_VECTOR_INI(vectorProjects[i],IntTriple,cutnelem);

    VecInt **sameJob;
    ALLOCATE_VECTOR_INI(sameJob,VecInt*,Inst_nJobs(inst));
    for(int i = 0; i < Inst_nJobs(inst) ; i++)
        sameJob[i] = VInt_create();

    for(int l = 0; l < cutnelem; l++) {
        conflicts[l] = VInt_create();
        ALLOCATE_VECTOR_INI(varInConfs[l], int, cutnelem);

        int j,m,t;
        lp_col_name(mipCut,VInt_get(idxVars,l),name);
        parseName( name, prefix, id );
        j = id[0];
        m = id[1];
        const Job *job = Inst_job(inst,j);
//        const Mode *mode = Job_mode(job,m);
        t = id[2];

        //incluir restrição da maximização // valor fracionário do original
        VInt_pushBack(sameJob[j], l);
        int plp = Job_project(job);

        vectorProjects[plp][nVProj[plp]].idx = l;
        vectorProjects[plp][nVProj[plp]].j = j;
        vectorProjects[plp][nVProj[plp]].m = m;
        vectorProjects[plp][nVProj[plp]].t = t;
        nVProj[plp]++;

    }



    for(int i = 0; i < Inst_nJobs(inst) ; i++) {
        int size = VInt_size(sameJob[i]);
        if(size<=0) continue;
        //if(jobsvar[] != 1) continue;
//        int nme=0;
        double coefsj[size];
        // printf("\nElem for job %d\n", i);
        for(int e = 0 ; e < size ; e++ )
            coefsj[e] = 1.0;
        double rhsm = 1;
        char nameR[STR_SIZE];
        sprintf( nameR, "cSelmode(%d)", i);
        lp_add_row(mipCut,size,VInt_getPtr(sameJob[i]),coefsj,nameR,'L',rhsm);

    }

    for(int p = 0 ; p < nProj; p++) {
        //  if(cont >= lp_rows(lp)) break;
        for(int elem = 0 ; elem < nVProj[p] ; elem++) {
            //if(cont >= lp_rows(lp)) break;
            int l = vectorProjects[p][elem].idx;
            int j = vectorProjects[p][elem].j, m = vectorProjects[p][elem].m, t =  vectorProjects[p][elem].t;
        //        printf("\n E (%d,%d,%d) %d ", j, m, t, l);
            const Job *job = Inst_job(inst,j);
            const Mode *mode = Job_mode(job,m);
            if(Mode_duration(mode)==0)continue;

            for(int elem2 = 0 ; elem2 < nVProj[p] ; elem2++) {
                //  if(cont >= lp_rows(lp)) break;
                if(elem == elem2 ) continue;
                int j2 = vectorProjects[p][elem2].j, m2 = vectorProjects[p][elem2].m, t2 = vectorProjects[p][elem2].t;
                if(j==j2 && m==m2 && t ==t2) continue;

                if(j==j2) {
                    int l2 = vectorProjects[p][elem2].idx;
                    if(varInConfs[l][l2]==1) continue;
                    VInt_pushBack(conflicts[l],l2);
              //          printf("Same Job (%d,%d,%d) %d ", j2, m2, t2, l2);
                    varInConfs[l][l2] = 1;
                    continue;
                }


                const Job *job2 = Inst_job(inst,j2);
                const Mode *mode2 = Job_mode(job2,m2);
                if(Mode_duration(mode2)==0)continue;

                int durationE = Mode_duration(mode);
                int timeEndPred = t+durationE;

                if( Job_hasIndSucc(job,j2) ) {
                    int value = Inst_getMaxDIJ(inst,j,j2)- Job_minDuration(job);
                    int winTime = value ;
                    if(timeEndPred>t2 || ((t2-timeEndPred) < winTime)) {
                        //     printf("timeEndPred = t %d + durationE %d \n", t, durationE );
                   //     printf("\ntimeEndPred %d > t2 %d || ((t2 %d - timeEndPred %d) < winTime %d)", timeEndPred,t2,t2,timeEndPred, winTime);
                        int l2 = vectorProjects[p][elem2].idx;
                        if(varInConfs[l][l2]==1) continue;
                        //     printf("prec D(%d,%d,%d) %d ", j2, m2, t2, l2);
                        VInt_pushBack(conflicts[l],l2);
                        varInConfs[l][l2] = 1;
                        continue;
                    }
                }
            }

            for(int p2 = 0 ; p2 < nProj ; p2++) {
                //  if(cont >= lp_rows(lp)) break;
                if(p==p2) continue;
                for(int elem2 = 0 ; elem2 < nVProj[p2] ; elem2++) {
                    //  if(cont >= lp_rows(lp)) break;
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
                    if(part1 + part2 > Inst_getSumTPD(inst)) {
                        //printf("\npart1 %d (t%d-cp1 %d) + part2 %d (t2%d-cp2 %d) > Inst_getSumTPD(inst) %d  releaseDatep1 %d + criticalPathp1 %d  DJMj1 %d releaseDatep2 %d + criticalPathp %d DIM2 %d", part1, t, cp1, part2, t2, cp2, Inst_getSumTPD(inst), Project_releaseDate(project),Project_criticalPath(project) , Inst_getMaxDIM(inst,job,Mode_index(mode)), Project_releaseDate(project2),Project_criticalPath(project2),Inst_getMaxDIM(inst,job2,Mode_index(mode2)));
                        int l2 = vectorProjects[p2][elem2].idx;
                        if(varInConfs[l][l2]==1) continue;
                        VInt_pushBack(conflicts[l],l2);
                         //        printf("tpd p(%d,%d,%d) %d ", j2, m2, t2, l2);
                        varInConfs[l][l2] = 1;

                    }
                }
            }
        }
    }


    int **conf;
    ALLOCATE_VECTOR(conf, int*, cutnelem);

    int *nConf;
    ALLOCATE_VECTOR_INI(nConf, int, cutnelem);
    int weight[cutnelem];
    //printf("find conflicts clique: end\n");fflush(stdout);
    for(int l = 0 ; l < cutnelem ; l++) {
        int sizeConflictsJ = VInt_size(conflicts[l]);
        nConf[l] = sizeConflictsJ;
        // printf("\nconflicts[%d] size %d: ", l, sizeConflictsJ);
        if(sizeConflictsJ != 0) {
            int *confJ = VInt_getPtr(conflicts[l]);
            double *coef;
            ALLOCATE_VECTOR_INI(coef, double, sizeConflictsJ);
            FILL(coef,0,sizeConflictsJ,1.0);
            ALLOCATE_VECTOR(conf[l], int, sizeConflictsJ);
            for(int o = 0 ; o <  sizeConflictsJ; o++) {

                //printf("confJ[o]%d", confJ[o]);
                conf[l][o]= confJ[o];
                //printf("conf[%d][%d] %d ", l,o, conf[l][o]);
                //getchar();
            }
            CutP_quick_sort_vec(conf[l],coef,sizeConflictsJ);
            if(xf[VInt_get(cutsElem,l)]<=0.00001)
                weight[l] = 1;
            else
                weight[l] = xf[VInt_get(cutsElem,l)]*1000.0;
            free(coef);
        } else
            weight[l] = 0.0;
    }

    _time = ( (double) timeLeft - (omp_get_wtime()-tinit) );
    CGraph *cgraph = build_cgraph_conflicts(conf, nConf, cutnelem, _time);

    cgraph_set_weight( cgraph, weight );
    BronKerbosch *bk =  bk_create(cgraph);
//    _time = ( (double) timeLeft - (omp_get_wtime()-tinit) );
  //  bk_set_min_weight(bk, weight);
    bk_set_timelimit(bk, _time);
    bk_run(bk);

    CliqueSet *clqSet = bk_get_clq_set(bk);
    int nCliqueSet =  clq_get_number_of_cliques(clqSet);

    int totalConf = 0;
    for(int nc = 0 ; nc < nCliqueSet ; nc++) {
        int* iset = clq_set_clique_elements(clqSet, nc);
        int size = clq_set_clique_size(clqSet, nc);
        double coef[size];
        FILL(coef,0,size,1);
        char nameR[STR_SIZE];
        sprintf( nameR, "cConf(%d)", nc);
        lp_add_row(mipCut,size,iset,coef,nameR,'L',1.0);
        totalConf++;
    }

    //lp_set_concurrentMIP(mipCut,1);
    //lp_set_method(mipCut,4);
    //lp_set_seed(mipCut,100000);
    int st = lp_optimize(mipCut);
    double newrhs = 0;
    double slack = 0.0;
    if(st==LP_OPTIMAL)
    {
            newrhs =lp_obj_value(mipCut)*-1;
            printf("printing mipCut\n");
            fflush(stdout);
            lp_write_lp(mipCut,"mipCut.lp");
            lp_write_sol(mipCut,"mipCut.sol");//getchar();

                /* adding cuts*/

            ALLOCATE_VECTOR(ccjs->cutElem,VecInt*,1);
            ALLOCATE_VECTOR(ccjs->cutCoef,VecDbl*,1);

            ccjs->cutElem[ccjs->nAllocations] = VInt_create();
            ccjs->cutCoef[ccjs->nAllocations] = VDbl_create();

            //VecInt *idx = VInt_create();
            double *xflm = lp_x(mipCut);
            int contelem = 0;
            for(int l = 0 ; l < lp_cols(mipCut); l++){
                if (fabs(xflm[l])>1e-5 ) {
                    char vnamex[STR_SIZE];
                    lp_col_name(mipCut,l,vnamex);
                    if (tolower(vnamex[0])=='x') {
                        parseName( vnamex, prefix, id );
                        if (prefix[0]=='x') {
                            int elem = lp_col_index(lp,vnamex);
                            VInt_pushBack(ccjs->cutElem[ccjs->nAllocations],elem);
                            VDbl_pushBack(ccjs->cutCoef[ccjs->nAllocations],id[2]*xf[elem]);
                            slack += id[2]*xf[elem]*xflm[l];
                            contelem++;
                            printf( " %s :  %f * %d, slack %f \n", vnamex, id[2]*xf[elem],elem, slack);
                        }
                    }
                }
            }
            char nameCut[STR_SIZE];
            sprintf(nameCut, "liftJobSet(%d)", lp_cols(lp) );
            printf("liftJobSet(%d)", lp_cols(lp) );
            VStr_pushBack(ccjs->cutname,nameCut);
            slack = slack-newrhs;
            int *cutIdx = VInt_getPtr(ccjs->cutElem[ccjs->nAllocations]);
            double *cutCoef = VDbl_getPtr(ccjs->cutCoef[ccjs->nAllocations]);
            VDbl_pushBack(ccjs->cutrhs,newrhs);
            VInt_pushBack(ccjs->cutsense,0);
            VInt_pushBack(ccjs->cutdominated,0);
            VInt_pushBack(ccjs->cutnelem,contelem);
            VDbl_pushBack(ccjs->cutviolation,slack);
            CutP_quick_sort_vec(cutIdx, cutCoef, contelem);
            ccjs->nAllocations++;
            printf( " rhs%f, sense %d, nelem %d, viol %f \n",newrhs, 0,contelem,slack);
          //  getchar();
    }

    //limite de nós e lower_bound valido arredondado para cima se nao achar o ótimo.

    for(int i = 0; i < Inst_nJobs(inst) ; i++)
        VInt_free(&sameJob[i]);
    free(sameJob);

//gurobi
//barrier

//precisao quad
//setar o quad -1
//numericfocus 3 ok

    VStr_free(&namesY);
    VStr_free(&namesZ);
    VStr_free(&namesX);
    VStr_free(&namesW);

    for(int l = 0 ; l < cutnelem ; l++) {
        if( nConf[l]!=0)  free(conf[l]);
        VInt_free(&conflicts[l]);
        free(varInConfs[l]);
    }

    free(conf);
    free(nConf);
    free(conflicts);
    free(varInConfs);


    for(int res = 0 ; res < Inst_nResR(inst); res++){

        for(int td = 0; td < nTimes ;td++){
            free(coefR[res][td]);
        }
        free(coefR[res]);
        free(contelemcoefrt[res]);
    }
    free(coefR);
    free(contelemcoefrt);

    for(int res = 0 ; res < Inst_nResN(inst); res++)
        free(coefNR[res]);
    free(coefNR);


    for(int p = 0 ; p < nProj ; p++)
        free(vectorProjects[p]);
    free(vectorProjects);
    free(nVProj);
    cgraph_free(&cgraph);
    bk_free(bk);
    lp_free(&mipCut);
}
