
/* GOAL PSP Solver II
 * Resource-Constrained Project Scheduling Problems (RCPSP)
 *
 * Develop as part of the D.Sc. thesis of Araujo, Janniele A. S., with collaboration
 *                                   of Santos, H.G.
 */

#include "cut_cg.h"
#define VERBOSE  1
#define PERCMAXCLIQUE  0.5

/*create and initialize the structure of cutCG*/
CutCG *CutCG_create(const Instance *inst, int nRow, int nElem, int nCols)
{

    CutCG *ccg;
    ALLOCATE_INI(ccg,CutCG);

    ccg->xfElemPP = VDbl_create();
    ccg->idxElemPP = VInt_create();
    ccg->idxelementsvec = VInt_create();
    ccg->rowrhs = VDbl_create();
    ccg->rowtype = VInt_create();
    ccg->rowtimeslotmin = VInt_create();
    ccg->rowtimeslotmax = VInt_create();

    ALLOCATE_VECTOR_INI(ccg->rowElem, VecInt*, nRow);
    ALLOCATE_VECTOR_INI(ccg->rowCoef, VecDbl*, nRow);

    for(int r = 0 ; r < nRow  ; r++) {
        ccg->rowElem[r] = VInt_create();
        ccg->rowCoef[r] = VDbl_create();
    }

    ccg->cont = nRow;
    ccg->contElem = nElem;
    ccg->nCols = nCols;

    ALLOCATE_VECTOR_INI(ccg->elements,int,nCols);
    ALLOCATE_VECTOR_INI(ccg->idxelements,int,nCols);
    ALLOCATE_VECTOR(ccg->elemrow,VecInt*,nCols);
    ALLOCATE_VECTOR(ccg->coefelemrows,VecDbl*,nCols);
    for ( int c=0 ; (c<nCols) ; ++c ) {
        ccg->elemrow[c] = VInt_create();
        ccg->coefelemrows[c] = VDbl_create();
    }

    ccg->rname = VStr_create(256);
    ccg->nameelements = VStr_create(256);
    ccg->rrhs = VDbl_create();

    ccg->nAllocations = 0;//1
    ccg->cutrhs = VDbl_create();
    ccg->cutnelem = VInt_create();
    ccg->cutviolation = VDbl_create();
    ccg->cutsense = VInt_create();
    ccg->cutname = VStr_create(256);
    ccg->cutdominated = VInt_create();

    ccg->contUR=0; ccg->contURN = 0; ccg->contUM =0; ccg->contUCI = 0; ccg->contUCP = 0; ccg->nJumps =0;

    return ccg;
}

/*create MIP to CG*/
LinearProgram * lp_create_cgsep(const double *xf, int nr, int nelem, int *idxelements, VecStr *nameelements, VecInt **elemrow, VecDbl **coefelemrows, VecDbl *rrhs, VecStr *rname, double timeleft)
{
    double _time;
    double startT = omp_get_wtime();


    LinearProgram *mipCGSep = lp_create();
    lp_set_print_messages(mipCGSep,0);
    if(VERBOSE==3) lp_set_print_messages(mipCGSep,1);
    lp_set_mip_emphasis( mipCGSep, LP_ME_FEASIBILITY );
    _time = ( (double) timeleft - (omp_get_wtime()-startT) );
    //printf("time %f", _time);fflush(stdout);
    lp_set_max_seconds(mipCGSep,(int)_time);


    //lp_set_max_nodes(mipCGSep,10000);//GRB_DBL_PAR_NODELIMIT
    //printf("antes callback\n"); fflush(stdout);
#ifdef GRB
    getcallback(mipCGSep,0);
#endif // GRB
#ifdef CBC
    getcallback(mipCGSep,0);
#endif //

    //printf("passou callback\n"); fflush(stdout);

    VecStr *names = VStr_create( STR_SIZE );
    VecStr *namesU = VStr_create( STR_SIZE );
    VecStr *namesF = VStr_create( STR_SIZE );
    char name[STR_SIZE];
    sprintf(name, "aRHS");
    VStr_pushBack( names, name );
    sprintf(name, "fRHS");
    VStr_pushBack( namesF, name );

    for(int n = 0 ; n < nelem ; n++) {
        //        int e = idxelements[n];
        char namevar[STR_SIZE];
        sprintf(namevar, "%s", VStr_get(nameelements,n));
        char cname[STR_SIZE];
        sprintf(cname, "a%s", namevar);
        VStr_pushBack( names, cname );
        char fname[STR_SIZE];
        sprintf(fname, "f%s", namevar);
        VStr_pushBack( namesF, fname );
    }

    int nVars = VStr_size( names );

    double *obj, *lb, *ub;//obj[nVars], lb[nVars], ub[nVars];
    ALLOCATE_VECTOR_INI(obj, double, nVars);
    ALLOCATE_VECTOR_INI(lb, double, nVars);
    ALLOCATE_VECTOR_INI(ub, double, nVars);
    char *integer;
    ALLOCATE_VECTOR_INI(integer, double, nVars);

    FILL( lb, 0, nVars, -INT_MAX );
    FILL( ub, 0, nVars, INT_MAX );

    FILL( integer, 0, nVars, 1.0 );

    int var = 0;
    obj[var] = 1.0;
    for(int n = 0 ; n < nelem ; n++) {
        int e = idxelements[n];
        //    int e = VInt_get(idxelements,n);
        var++;
        obj[var] = -xf[e];
    }

    lp_add_cols( mipCGSep, VStr_size(names), obj, lb, ub, integer, VStr_ptr(names) );

    for(int i = 0 ; i < nr ; i++ ) {
        char namevar[STR_SIZE];
        sprintf(namevar, "%s",VStr_get(rname,i));
        char name[STR_SIZE];
        sprintf(name, "u%s", namevar);
        // printf( "u%s", namevar); fflush(stdout);
        VStr_pushBack( namesU, name );
    }

    int nVarsF = VStr_size( namesF );
    int nVarsU = VStr_size( namesU );
    double *objU, *lbU, *ubU, *objF, *lbF, *ubF;//obj[nVars], lb[nVars], ub[nVars];
    ALLOCATE_VECTOR_INI(objF, double, nVarsF);
    ALLOCATE_VECTOR_INI(lbF, double, nVarsF);
    ALLOCATE_VECTOR_INI(ubF, double, nVarsF);
    ALLOCATE_VECTOR_INI(objU, double, nVarsU);
    ALLOCATE_VECTOR_INI(lbU, double, nVarsU);
    ALLOCATE_VECTOR_INI(ubU, double, nVarsU);

    char *integerU;
    ALLOCATE_VECTOR_INI(integerU, double, nVarsU);
    char *integerF;
    ALLOCATE_VECTOR_INI(integerF, double, nVarsF);

    FILL( objU, 0, nVarsU, 0.0001);//EPS );
    FILL( objF, 0, nVarsF, 0.0 );
    FILL( lbU, 0, nVarsU, 0.0 );
    //FILL( ubU, 0, nVarsU, 1 );
    FILL( ubU, 0, nVarsU, 0.99 );
    FILL( lbF, 0, nVarsF, 0.0 );
    FILL( ubF, 0, nVarsF, 0.99 );
    FILL( integerU, 0, nVarsU, 0.0 );
    FILL( integerF, 0, nVarsF, 0.0 );

    lp_add_cols( mipCGSep, VStr_size(namesU), objU, lbU, ubU, integerU, VStr_ptr(namesU) );
    lp_add_cols( mipCGSep, VStr_size(namesF), objF, lbF, ubF, integerF, VStr_ptr(namesF) );


    for(int n = 0 ; n < nelem ; n++) {
        int e = idxelements[n];
        //int e = VInt_get(idxelements,n);
        int size = VInt_size(elemrow[e]);
        int idx[size+2];
        double coef[size+2];
        int ne = 0;
        for(int rj = 0 ; rj < size ; rj ++) {
            double coeferj = VDbl_get(coefelemrows[e],rj);
            int i = VInt_get(elemrow[e], rj);
            //  printf("ne %d i %d coeferj %f\n", ne, i, coeferj); fflush(stdout);
            char namevar[STR_SIZE];
            sprintf(namevar, "%s",VStr_get(rname,i));
            char name[STR_SIZE];
            sprintf(name, "u%s", namevar);
            int ui = lp_col_index(mipCGSep,name);
            idx[ne] = ui;
            coef[ne] = coeferj;
            ne++;
        }
        char namevar[STR_SIZE];
        sprintf(namevar,"%s",VStr_get(nameelements,n));
        sprintf(name, "a%s", namevar);
        int aj = lp_col_index(mipCGSep,name);
        idx[ne] = aj;
        coef[ne] = -1;
        //  printf("ne %d aj %d coeferj -1\n", ne,  aj);
        ne++;
        sprintf(name, "f%s",namevar);
        int fj = lp_col_index(mipCGSep,name);
        idx[ne] = fj;
        coef[ne] = -1;
        //   printf("ne %d fj %d coeferj -1\n", ne,  fj);
        ne++;
        char namer[STR_SIZE];
        sprintf(namer, "ConstraintURE%s#%d", namevar, lp_rows(mipCGSep)+1);
        lp_add_row( mipCGSep, ne, idx, coef, namer, 'E', 0.0 );

    }

    int idx[nr+2];
    double coef[nr+2];

    int ni=0;
    for(int i = 0 ; i < nr ; i ++) {
        char namevar[STR_SIZE];
        sprintf(namevar, "%s",VStr_get(rname,i));
        //        printf("namevar %s\n", namevar); fflush(stdout);
        char name[STR_SIZE];
        sprintf(name, "u%s", namevar);
        int ui = lp_col_index(mipCGSep,name);
        double rhsi = VDbl_get(rrhs,i);
        idx[ni] = ui;
        coef[ni] =  rhsi;
        //  printf("nr %d ui %d rhsi %f\n", nr, ui, rhsi);
        ni++;
    }
    sprintf(name, "aRHS");
    int a0 = lp_col_index(mipCGSep,name);
    idx[ni] = a0;
    coef[ni] = -1;
    // printf("nr %d aj %d f0 -1 \n", nr, a0);
    ni++;
    sprintf(name, "fRHS");
    int f0 = lp_col_index(mipCGSep,name);
    idx[ni] = f0;
    coef[ni] = -1;
    // printf("nr %d fj %d f0 -1 \n", nr, f0k);
    ni++;
    char namer[STR_SIZE];
    sprintf(namer, "ConstraintUR#%d", lp_rows(mipCGSep)+1);
    lp_add_row( mipCGSep, ni, idx, coef, namer, 'E', 0.0 );

    int idx3[ nelem];
    double coef3[ nelem];
    int ne3=0;
    for(int n = 0 ; n < nelem ; n++) {
        int e = idxelements[n];
        //int e = VInt_get(idxelements,n);
        if(fabs( xf[e])<1e-5)continue;
        int size = VInt_size(elemrow[e]);
        if(size==0) continue;
        char namevar[STR_SIZE];
        sprintf(namevar, "%s", VStr_get(nameelements,n));
        sprintf(name, "a%s", namevar);
        int aj3 = lp_col_index(mipCGSep,name);
        idx3[ ne3] = aj3;
        coef3[ ne3] = xf[e];
        ne3++;
    }
    sprintf(name, "aRHS");
    int a03 = lp_col_index(mipCGSep,name);
    idx3[ ne3] = a03;
    coef3[ ne3] = -1;
    ne3++;

    char namer3[STR_SIZE];
    sprintf(namer3, "ConstraintURViol#%d", lp_rows(mipCGSep)+1);
    lp_add_row( mipCGSep, ne3, idx3, coef3, namer3, 'G', 0.0002);//+EPS ); //+EPS


    free(obj);
    free(lb);
    free(ub);
    free(integer);
    VStr_free( &names );
    free(objF);
    free(lbF);
    free(ubF);
    free(integerF);
    VStr_free( &namesF );
    free(objU);
    free(lbU);
    free(ubU);
    free(integerU);
    VStr_free( &namesU );

    return mipCGSep;

}


/*identify the set of important rows to compose the CG procedure with original constraints about renewable resources*/
CutCG *CutCG_create_and_identify_rows(  LinearProgram *lp, const int *origCols, LinearProgram *origLP, const Instance *inst, double timeLeft, int continuous,  int **maxTJM,  double maxcuts, int mininstant, int maxinstant, int jump, int nround)
{

    assert(lp);
    assert(inst);
    assert(origLP);

    //    int nround = Res_getRound(res);

    const double *xf = lp_x(lp);
    //double maxRC = lp_get_max_reduced_cost(origLP);
    //const double *rdc = lp_reduced_cost(lp);
    //  const double *rdc2 = lp_reduced_cost(origLP);
    char name[256];

    /* maximum number of indexes in a var name */
#define MAX_IDX 16
    char prefix[STR_SIZE];
    int id[MAX_IDX];
    int maxT = 0;
    int ncols =lp_cols(lp);

    int nJobs = Inst_nJobs(inst);
    int nModes = Inst_nMaxModes(inst);
    // checking maximum time  to all variables //where there is some allocation ** if
    for ( int i=0 ; (i<lp_cols(lp)) ; ++i ) {
        if (fabs(xf[i])>1e-5 ) {
            // if (fabs(xf[i])>1e-5 || maxRC != -1.0) {
            //   if(maxRC != -1.0)
            //     if(rdc[i] > maxRC) continue;

            lp_col_name( origLP, origCols[i], name );
            if (tolower(name[0])=='x') {
                parseName( name, prefix, id );
                if (prefix[0]=='x') {
                    //    int j = id[0];
                    //   int m = id[1];
                    int t = id[2];
                    maxT = MAX( maxT, t );
                    //      printf( "x(%d,%d,%d) xf: %g xfor %g, rdc: %f \n", j, m, t, xf[i], xf[origCols[i]], rdc[i]);
                }
            }
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

    // filling x from fractional solution
    for ( int i=0 ; (i<lp_cols(lp)) ; ++i ) {
        if (fabs(xf[i])>1e-5 ) {
            //if (fabs(xf[i])>1e-5 || maxRC != -1.0) {
            //  if(maxRC != -1.0)
            //    if(rdc[i] > maxRC) continue;

            lp_col_name( origLP, origCols[i], name );
            if (tolower(name[0]=='x')) {
                parseName( name, prefix, id );
                if (prefix[0]=='x')
                    x[id[0]][id[1]][id[2]] = xf[i];

            }
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
    for ( int i=0 ; (i<lp_cols(lp)) ; ++i ) {
        if (fabs(xf[i])>1e-5 ) {
            //if (fabs(xf[i])>1e-5 || maxRC != -1.0) {
            // if(maxRC != -1.0)
            //      if(rdc[i] > maxRC) continue;

            lp_col_name( origLP, origCols[i], name );
            if (tolower(name[0])=='x') {
                parseName( name, prefix, id );
                if (prefix[0]=='x') {
                    int jlp = id[0];
                    int mlp =id[1];
                    int tlp = id[2];
                    xIdx[jlp][mlp][tlp] = i;
                    //     printf("xIdx[%d][%d][%d] = %d ", jlp,mlp,tlp,i);
                }
            }
        }
    }

    //  getchar();
    //CHOSE CONSTRAINTS

    //    printf("iniciou alocacao memory\n");// getchar();
    //    int lpr = lp_rows(lp);

    VecStr *rname = VStr_create(256);
    VecStr *nameelements = VStr_create(256);
    VecDbl * rrhs = VDbl_create();
    VecInt * rtype = VInt_create();
    VecInt * rtimeslotmin = VInt_create();
    VecInt * rtimeslotmax = VInt_create();

    /*VecChar * rsense = VChar_create();
    VecInt ** ridx;
    ALLOCATE_VECTOR(ridx,VecInt*,lp_rows(lp));
    VecDbl ** rcoef;
    ALLOCATE_VECTOR(rcoef,VecDbl*,lp_rows(lp));

    for ( int r=0 ; (r<lp_rows(lp)) ; ++r ) {
        ridx[r] = VInt_create();
        rcoef[r] = VDbl_create();
    }
    */
//    CutP_create_job_set(inst,lp, maxT);

    //    int contUR=0, contURN = 0, contUM =0, contUCI = 0, contUCP = 0;
    int nJumps =0;
    int *elements;
    ALLOCATE_VECTOR_INI(elements,int,lp_cols(lp));
    int *idxelements;
    ALLOCATE_VECTOR_INI(idxelements,int,lp_cols(lp));
    VecInt **elemrow;
    ALLOCATE_VECTOR(elemrow,VecInt*,lp_cols(lp));
    VecDbl **coefelemrows;
    ALLOCATE_VECTOR(coefelemrows,VecDbl*,lp_cols(lp));
    for ( int c=0 ; (c<lp_cols(lp)) ; ++c ) {
        elemrow[c] = VInt_create();
        coefelemrows[c] = VDbl_create();
    }

    // printf("passou alocacao memory\n"); fflush(stdout); //getchar();

    int cont = 0, contElem = 0;
    int mininterval = mininstant, maxinterval = maxinstant;
    double maxvalue = 0;
    double sumvalue = 0;
    double perc = 1.25;

    //if(maxinstant >= maxT) maxinstant = maxT;
    for ( int t=mininstant; (t<maxinstant && maxinstant <= maxT) ; ++t ) {
        for ( int j=0; (j<nJobs) ; ++j ) {
            for ( int m=0 ; (m<nModes) ; ++m ) {
                int ind = xIdx[j][m][t];
                if( ind == -1) continue;
                if (fabs(xf[ind])>1e-5 ) {
                    //if (fabs(xf[ind])>1e-5 || maxRC != -1.0) {
                    //  if(maxRC != -1.0)
                    //    if(rdc[ind] > maxRC) continue;

                    double intcloser = ROUND(xf[ind]);
                    double value = fabs(intcloser-xf[ind]);
                    //printf("int closer %f, xf[%d] %f, value %f\n", intcloser, ind, xf[ind], value);
                    sumvalue +=value;
                }
            }
        }
        if(t==maxinstant-1) {

            if(sumvalue >= (double)(maxvalue)*perc ) {
                mininterval = mininstant;
                maxinterval = maxinstant;
                maxvalue = sumvalue;
                nJumps++;
                //   printf("sumvalue %f, maxvalue %f\n",  sumvalue, maxvalue );
            }
            //   printf("maxT %d, min interval %d, max interval %d, min instant %d, max instant %d, jump %d , sumvalue %f, maxvalue %f\n",maxT,mininterval, maxinterval, mininstant, maxinstant, jump, sumvalue, maxvalue );
            //if(VERBOSE==3)
            sumvalue = 0;
            mininstant += jump;
            maxinstant += jump;
            if(maxinstant >= maxT) maxinstant = maxT;
            if(jump !=0 )
                t = mininstant;
        }
    }

    for ( int r=0 ; (r<lp_rows(lp)) ; ++r ) {
        char nameR[256];
        lp_row_name(lp, r, nameR);
        // printf("nameR %s \n", nameR);
        // getchar();
        if(strncmp(nameR,"resR",4)== 0 ){// || strncmp(nameR,"cutE",4)== 0){
            parseName( nameR, prefix, id );
            if(id[1] < mininterval || id[1] > maxinterval ) continue;
            //printf("mininstant %d maxinstant %d %s(%d,%d)\n", mininterval, maxinterval,nameR, id[0], id[1]);
            int *idx;
            ALLOCATE_VECTOR_INI(idx,int,lp_cols(lp));
            double *coef;
            ALLOCATE_VECTOR_INI(coef,double,lp_cols(lp));
            int size = lp_row(lp,r,idx,coef);
            int aux = 0;
            if(size<=0) {
                free(idx);
                free(coef);
                continue;
            }
            for(int i  = 0 ; i < size ; i++) {
                char nameE[256];
                lp_col_name(lp,idx[i],nameE);
                parseName( nameE, prefix, id );
                int ind = xIdx[id[0]][id[1]][id[2]];
                if(ind != -1) {

                    if(elements[ind]==0 ) {
                        elements[ind] = 1;
                        idxelements[contElem] = ind;
                        VStr_pushBack(nameelements, nameE);
                        contElem++;
                    }
                    //VInt_pushBack(ridx[cont],idx[i]);
                    //VDbl_pushBack(rcoef[cont],coef[i]);
                    VInt_pushBack(elemrow[ind], cont);
                    VDbl_pushBack(coefelemrows[ind], coef[i]);

                    aux++;
                }
            }
            if(aux>0) {
                double rhs = lp_rhs(lp,r);
                VDbl_pushBack(rrhs,rhs);
                if(strncmp(nameR,"resR",4)== 0) VInt_pushBack(rtype,RES_RR);
                //if(strncmp(nameR,"cutRR",5)== 0) VInt_pushBack(rtype,LPC_RR);
                //if(strncmp(nameR,"cutE",4)== 0) VInt_pushBack(rtype,LPC_ENERGETIC);

                VInt_pushBack(rtimeslotmin,mininterval);
                VInt_pushBack(rtimeslotmax,maxinterval);
                //char sense = lp_sense(lp,r);
                //VChar_pushBack(rsense,sense);
                VStr_pushBack(rname,nameR);
                cont++;
                //    printf("%s",nameR );
            }
            free(idx);
            free(coef);
        }
    }

    CutCG *ccg  = CutCG_create(inst, cont, contElem, ncols);

    for(int r = 0 ; r < cont  ; r++) {
        double rhs = VDbl_get(rrhs,r);
        VDbl_pushBack(ccg->rowrhs,rhs);
        // VDbl_pushBack(ccg->rrhs,rhs);
        int type = VInt_get(rtype,r);
        VInt_pushBack(ccg->rowtype,type);
        int timeslotmin = VInt_get(rtimeslotmin,r);
        int timeslotmax = VInt_get(rtimeslotmax,r);
        VInt_pushBack(ccg->rowtimeslotmin,timeslotmin);
        VInt_pushBack(ccg->rowtimeslotmax,timeslotmax);
        VStr_pushBack(ccg->rname, VStr_get(rname,r));
    }

    for(int e = 0 ; e < contElem ; e++) {
        int idxel = idxelements[e];
        ccg->idxelements[e] = idxelements[e];
        VDbl_pushBack(ccg->xfElemPP,xf[idxel]);
        VInt_pushBack(ccg->idxElemPP,e);
        VInt_pushBack(ccg->idxelementsvec,idxel);

        VStr_pushBack(ccg->nameelements,VStr_get(nameelements,e));
    }


    for(int el = 0 ; el < contElem ; el++) {
        int elem = idxelements[el];
        int elemPP = VInt_get(ccg->idxElemPP,el);
        int nRows = VInt_size(elemrow[elem]);
        for(int i  = 0 ; i < nRows ; i++ ) {
            int row = VInt_get(elemrow[elem], i);
            VInt_pushBack(ccg->elemrow[elem],row);
            double coef = VDbl_get(coefelemrows[elem], i);
            VDbl_pushBack(ccg->coefelemrows[elem], coef);
            VInt_pushBack(ccg->rowElem[row],elemPP);
            VDbl_pushBack(ccg->rowCoef[row],coef);
        }
    }


    free(elements);
    free(idxelements);
    VStr_free(&rname);
    VDbl_free(&rrhs);
    VInt_free(&rtype);
    VInt_free(&rtimeslotmin);
    VInt_free(&rtimeslotmax);

    for ( int c=0 ; (c<lp_cols(lp)) ; ++c ) {
        VInt_free(&elemrow[c]);
         VDbl_free(&coefelemrows[c]);
    }

    free(coefelemrows);
    free(elemrow);
    VStr_free(&nameelements);

    free(xIdx[0][0]);
    free(xIdx[0]);
    free(xIdx);
    free(x[0][0]);
    free(x[0]);
    free(x);


    return ccg;

}

/*identify the set of important rows to compose the CG procedure with all constraints (RR, cutRR, NR, Job Selection)*/
CutCG *CutCG_create_and_identify_rows_all(  LinearProgram *lp, const int *origCols, LinearProgram *origLP, const Instance *inst, double timeLeft, int mininstant, int maxinstant, int jump, int nround)
{

  //  double _time;
//    double tinit = omp_get_wtime();

    assert(lp);
    assert(inst);
    assert(origLP);

    //    int nround = Res_getRound(res);

    const double *xf = lp_x(lp);

    //double maxRC = lp_get_max_reduced_cost(origLP);
    //const double *rdc = lp_reduced_cost(lp);
    //  const double *rdc2 = lp_reduced_cost(origLP);
    char name[256];

    /* maximum number of indexes in a var name */
#define MAX_IDX 16
    char prefix[STR_SIZE];
    int id[MAX_IDX];
    int maxT = 0;
    int ncols =lp_cols(lp);

    int nJobs = Inst_nJobs(inst);
    int nModes = Inst_nMaxModes(inst);
    // checking maximum time  to all variables //where there is some allocation ** if
    for ( int i=0 ; (i<lp_cols(lp)) ; ++i ) {
        if (fabs(xf[i])>1e-5 ) {
            // if (fabs(xf[i])>1e-5 || maxRC != -1.0) {
            //   if(maxRC != -1.0)
            //     if(rdc[i] > maxRC) continue;

            lp_col_name( origLP, origCols[i], name );
            if (tolower(name[0])=='x') {
                parseName( name, prefix, id );
                if (prefix[0]=='x') {
                    //    int j = id[0];
                    //   int m = id[1];
                    int t = id[2];
                    maxT = MAX( maxT, t );
                    //      printf( "x(%d,%d,%d) xf: %g xfor %g, rdc: %f \n", j, m, t, xf[i], xf[origCols[i]], rdc[i]);
                }
            }
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

    // filling x from fractional solution
    for ( int i=0 ; (i<lp_cols(lp)) ; ++i ) {
        if (fabs(xf[i])>1e-5 ) {
            //if (fabs(xf[i])>1e-5 || maxRC != -1.0) {
            //  if(maxRC != -1.0)
            //    if(rdc[i] > maxRC) continue;

            lp_col_name( origLP, origCols[i], name );
            if (tolower(name[0]=='x')) {
                parseName( name, prefix, id );
                if (prefix[0]=='x')
                    x[id[0]][id[1]][id[2]] = xf[i];

            }
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
    for ( int i=0 ; (i<lp_cols(lp)) ; ++i ) {
        if (fabs(xf[i])>1e-5 ) {
            //if (fabs(xf[i])>1e-5 || maxRC != -1.0) {
            // if(maxRC != -1.0)
            //      if(rdc[i] > maxRC) continue;

            lp_col_name( origLP, origCols[i], name );
            if (tolower(name[0])=='x') {
                parseName( name, prefix, id );
                if (prefix[0]=='x') {
                    int jlp = id[0];
                    int mlp =id[1];
                    int tlp = id[2];
                    xIdx[jlp][mlp][tlp] = i;
                    //     printf("xIdx[%d][%d][%d] = %d ", jlp,mlp,tlp,i);
                }
            }
        }
    }

    //  getchar();
    //CHOSE CONSTRAINTS

    //    printf("iniciou alocacao memory\n");// getchar();
    //    int lpr = lp_rows(lp);

    VecStr *rname = VStr_create(256);
    VecStr *nameelements = VStr_create(256);
    VecDbl * rrhs = VDbl_create();
    VecInt * rtype = VInt_create();
    VecInt * rtimeslotmin = VInt_create();
    VecInt * rtimeslotmax = VInt_create();

    /*VecChar * rsense = VChar_create();
    VecInt ** ridx;
    ALLOCATE_VECTOR(ridx,VecInt*,lp_rows(lp));
    VecDbl ** rcoef;
    ALLOCATE_VECTOR(rcoef,VecDbl*,lp_rows(lp));

    for ( int r=0 ; (r<lp_rows(lp)) ; ++r ) {
        ridx[r] = VInt_create();
        rcoef[r] = VDbl_create();
    }
    */
//    CutP_create_job_set(inst,lp, maxT);

    //    int contUR=0, contURN = 0, contUM =0, contUCI = 0, contUCP = 0;
    int nJumps =0;
    int *elements;
    ALLOCATE_VECTOR_INI(elements,int,lp_cols(lp));
    int *idxelements;
    ALLOCATE_VECTOR_INI(idxelements,int,lp_cols(lp));
    VecInt **elemrow;
    ALLOCATE_VECTOR(elemrow,VecInt*,lp_cols(lp));
    VecDbl **coefelemrows;
    ALLOCATE_VECTOR(coefelemrows,VecDbl*,lp_cols(lp));
    for ( int c=0 ; (c<lp_cols(lp)) ; ++c ) {
        elemrow[c] = VInt_create();
        coefelemrows[c] = VDbl_create();
    }

    // printf("passou alocacao memory\n"); fflush(stdout); //getchar();

    int cont = 0, contElem = 0;
    int mininterval = mininstant, maxinterval = maxinstant;
    double maxvalue = 0;
    double sumvalue = 0;
    //if(maxinstant >= maxT) maxinstant = maxT;
    for ( int t=mininstant; (t<maxinstant && maxinstant <= maxT) ; ++t ) {
        for ( int j=0; (j<nJobs) ; ++j ) {
            for ( int m=0 ; (m<nModes) ; ++m ) {
                int ind = xIdx[j][m][t];
                if( ind == -1) continue;
                if (fabs(xf[ind])>1e-5 ) {
                    //if (fabs(xf[ind])>1e-5 || maxRC != -1.0) {
                    //  if(maxRC != -1.0)
                    //    if(rdc[ind] > maxRC) continue;

                    double intcloser = ROUND(xf[ind]);
                    double value = fabs(intcloser-xf[ind]);
                    //printf("int closer %f, xf[%d] %f, value %f\n", intcloser, ind, xf[ind], value);
                    sumvalue +=value;
                }
            }
        }
        if(t==maxinstant-1) {

            if(sumvalue >= (maxvalue)*1.25 ) {
                mininterval = mininstant;
                maxinterval = maxinstant;
                maxvalue = sumvalue;
                nJumps++;
                //   printf("sumvalue %f, maxvalue %f\n",  sumvalue, maxvalue );
            }
            //   printf("maxT %d, min interval %d, max interval %d, min instant %d, max instant %d, jump %d , sumvalue %f, maxvalue %f\n",maxT,mininterval, maxinterval, mininstant, maxinstant, jump, sumvalue, maxvalue );
            //if(VERBOSE==3)
            sumvalue = 0;
            mininstant += jump;
            maxinstant += jump;
            if(maxinstant >= maxT) maxinstant = maxT;
            if(jump !=0 )
                t = mininstant;
        }
    }

    for ( int r=0 ; (r<lp_rows(lp)) ; ++r ) {
        char nameR[256];
        lp_row_name(lp, r, nameR);
        // printf("nameR %s \n", nameR);
        // getchar();
        if(strncmp(nameR,"resR",4)== 0 || strncmp(nameR,"cutRR",5)== 0){// || strncmp(nameR,"cutE",4)== 0){
            parseName( nameR, prefix, id );
            if(id[1] < mininterval || id[1] > maxinterval ) continue;
            //printf("mininstant %d maxinstant %d %s(%d,%d)\n", mininterval, maxinterval,nameR, id[0], id[1]);
            int *idx;
            ALLOCATE_VECTOR_INI(idx,int,lp_cols(lp));
            double *coef;
            ALLOCATE_VECTOR_INI(coef,double,lp_cols(lp));
            int size = lp_row(lp,r,idx,coef);
            int aux = 0;
            if(size<=0) {
                free(idx);
                free(coef);
                continue;
            }
            for(int i  = 0 ; i < size ; i++) {
                char nameE[256];
                lp_col_name(lp,idx[i],nameE);
                parseName( nameE, prefix, id );
                int ind = xIdx[id[0]][id[1]][id[2]];
                if(ind != -1) {

                    if(elements[ind]==0 ) {
                        elements[ind] = 1;
                        idxelements[contElem] = ind;
                        VStr_pushBack(nameelements, nameE);
                        contElem++;
                    }
                    //VInt_pushBack(ridx[cont],idx[i]);
                    //VDbl_pushBack(rcoef[cont],coef[i]);
                    VInt_pushBack(elemrow[ind], cont);
                    VDbl_pushBack(coefelemrows[ind], coef[i]);

                    aux++;
                }
            }
            if(aux>0) {
                double rhs = lp_rhs(lp,r);
                VDbl_pushBack(rrhs,rhs);
                if(strncmp(nameR,"resR",4)== 0) VInt_pushBack(rtype,RES_RR);
                if(strncmp(nameR,"cutRR",5)== 0) VInt_pushBack(rtype,LPC_RR);
                //if(strncmp(nameR,"cutE",4)== 0) VInt_pushBack(rtype,LPC_ENERGETIC);

                VInt_pushBack(rtimeslotmin,mininterval);
                VInt_pushBack(rtimeslotmax,maxinterval);
                //char sense = lp_sense(lp,r);
                //VChar_pushBack(rsense,sense);
                VStr_pushBack(rname,nameR);
                cont++;
                //    printf("%s",nameR );
            }
            free(idx);
            free(coef);
        }
    }
//RES_RR,LPC_RR,RES_NR,RES_MODE
    VecInt **sameJob;
    ALLOCATE_VECTOR_INI(sameJob,VecInt*,Inst_nJobs(inst));
    for(int i = 0; i < Inst_nJobs(inst) ; i++)
        sameJob[i] = VInt_create();
    VecInt **sameResN,**sameResNCoef;
    ALLOCATE_VECTOR_INI(sameResN,VecInt*,Inst_nResN(inst));
    ALLOCATE_VECTOR_INI(sameResNCoef,VecInt*,Inst_nResN(inst));
    for(int i = 0; i < Inst_nResN(inst); i++) {
        sameResN[i] = VInt_create();
        sameResNCoef[i] = VInt_create();
    }


    for(int l = 0; l < ncols && cont > 0; l++) {

        if (fabs(xf[l])>1e-5 ) {
            //if (fabs(xf[l])>1e-5 || maxRC != -1.0) {
            //  if(maxRC != -1.0)
            //    if(rdc[l] > maxRC) continue;

            lp_col_name( lp,l, name );

            if (tolower(name[0])=='x') {
                parseName( name, prefix, id );
                if (prefix[0]=='x') {

                    int jlp = id[0];
                    int mlp =id[1];
                    int tlp = id[2];

                    int ind = xIdx[jlp][mlp][tlp];
                    if(ind == -1) continue;
                    //  printf("ind %d xIdx[%d][%d][%d];\n",ind,jlp,mlp,tlp);
                    if(elements[ind]==0) continue;

                    VInt_pushBack(sameJob[jlp], ind);

                    const Job *job = Inst_job(inst,jlp);
                    const Mode *mode = Job_mode(job,mlp);
                    if(Mode_duration(mode)==0)continue;


                    for(int mmm = 0 ; mmm < Mode_nResN(mode); mmm++) {
                        int idxResN = Mode_idxResN(mode,mmm);
                        VInt_pushBack(sameResN[idxResN], ind);
                        VInt_pushBack(sameResNCoef[idxResN], Mode_useResN(mode,mmm));
                    }
                }
            }
        }
    }


    for(int i = 0; i < Inst_nJobs(inst) ; i++) {
        int size = VInt_size(sameJob[i]);
        if(size<=0) continue;
        int nme=0;
        // printf("\nElem for job %d\n", i);
        for(int e = 0 ; e < size ; e++ ) {
            int elem = VInt_get(sameJob[i],e);
            //    printf("elem %d ", elem);
            VInt_pushBack(elemrow[elem], cont);
            VDbl_pushBack(coefelemrows[elem], 1);
            nme++;
        }
        if(nme>0) {
            double rhs = 1;
            VDbl_pushBack(rrhs,rhs);
            VInt_pushBack(rtype,RES_MODE);
            VInt_pushBack(rtimeslotmin,mininterval);
            VInt_pushBack(rtimeslotmax,maxinterval);
            char nameCC[256];
            sprintf(nameCC,"modeSel(%d)",i);
            // printf(" %s type %d \n",nameCC,  RES_MODE);
            VStr_pushBack(rname,nameCC);
            cont++;
        }
    }

    for(int r = 0; r < Inst_nResN(inst); r++) {
        if(cont >= lp_rows(lp)) break;
        int size = VInt_size(sameResN[r]);
        if(size<=0) continue;
        int nre =0 ;
        //printf("resource %d ", r);
        for(int elem = 0 ; elem < size ; elem++) {
            int l = VInt_get(sameResN[r],elem);
            double coef= (double) VInt_get(sameResNCoef[r],elem);
            //printf("%f*%d  ", coef, l);
            VInt_pushBack(elemrow[l], cont);
            VDbl_pushBack(coefelemrows[l], coef);
            nre++;
        }
        if(nre>0) {
            double rhs = (double) Inst_capResN(inst,r);
            //   printf("rhs %f\n", rhs);
            VDbl_pushBack(rrhs,rhs);
            VInt_pushBack(rtype,RES_NR);
            VInt_pushBack(rtimeslotmin,mininterval);
            VInt_pushBack(rtimeslotmax,maxinterval);

            char nameCC[256];
            sprintf(nameCC,"resN(%d)",r);
            //printf("%s type %d",nameCC, RES_NR );
            VStr_pushBack(rname,nameCC);
            cont++;
            nre=0;
        }
    }


    CutCG *ccg  = CutCG_create(inst, cont, contElem, ncols);

    for(int r = 0 ; r < cont  ; r++) {
        double rhs = VDbl_get(rrhs,r);
        VDbl_pushBack(ccg->rowrhs,rhs);
        // VDbl_pushBack(ccg->rrhs,rhs);
        int type = VInt_get(rtype,r);
        VInt_pushBack(ccg->rowtype,type);
        int timeslotmin = VInt_get(rtimeslotmin,r);
        int timeslotmax = VInt_get(rtimeslotmax,r);
        VInt_pushBack(ccg->rowtimeslotmin,timeslotmin);
        VInt_pushBack(ccg->rowtimeslotmax,timeslotmax);
        VStr_pushBack(ccg->rname, VStr_get(rname,r));
    }

    for(int e = 0 ; e < contElem ; e++) {
        int idxel = idxelements[e];
        ccg->idxelements[e] = idxelements[e];
        VDbl_pushBack(ccg->xfElemPP,xf[idxel]);
        VInt_pushBack(ccg->idxElemPP,e);
        VInt_pushBack(ccg->idxelementsvec,idxel);

        VStr_pushBack(ccg->nameelements,VStr_get(nameelements,e));
    }

    for(int el = 0 ; el < contElem ; el++) {
        int elem = idxelements[el];
        int elemPP = VInt_get(ccg->idxElemPP,el);
        int nRows = VInt_size(elemrow[elem]);
        for(int i  = 0 ; i < nRows ; i++ ) {
            int row = VInt_get(elemrow[elem], i);
            VInt_pushBack(ccg->elemrow[elem],row);
            double coef = VDbl_get(coefelemrows[elem], i);
            VDbl_pushBack(ccg->coefelemrows[elem], coef);
            VInt_pushBack(ccg->rowElem[row],elemPP);
            VDbl_pushBack(ccg->rowCoef[row],coef);
        }
    }

    free(elements);
    free(idxelements);
    VStr_free(&rname);
    VDbl_free(&rrhs);
    VInt_free(&rtype);
    VInt_free(&rtimeslotmin);
    VInt_free(&rtimeslotmax);
    for(int i = 0; i < Inst_nResN(inst); i++) {
        VInt_free(&sameResN[i]);
        VInt_free(&sameResNCoef[i]);
    }
    for(int j = 0 ; j < Inst_nJobs(inst) ; j++)
        VInt_free(&sameJob[j]);

    for ( int c=0 ; (c<lp_cols(lp)) ; ++c ) {
        VInt_free(&elemrow[c]);
        VDbl_free(&coefelemrows[c]);
    }

    free(sameJob);
    free(sameResN);
    free(sameResNCoef);
    free(coefelemrows);
    free(elemrow);

    VStr_free(&nameelements);

    free(xIdx[0][0]);
    free(xIdx[0]);
    free(xIdx);
    free(x[0][0]);
    free(x[0]);
    free(x);


    return ccg;

}

/*identify the set of important rows to compose the CG procedure with all constraints (RR, cutRR, NR, Job Selection,) and conflicts about precedence relationship and project delay*/
CutCG *CutCG_create_and_identify_rows_all_conflitcts(  LinearProgram *lp, const int *origCols, LinearProgram *origLP, const Instance *inst, double timeLeft, int mininstant, int maxinstant, int jump, int nround, int *horizon)
{

    double _time;
    double tinit = omp_get_wtime();

    assert(lp);
    assert(inst);
    assert(origLP);

    //    int nround = Res_getRound(res);

    const double *xf = lp_x(lp);
    //double maxRC = lp_get_max_reduced_cost(origLP);
    //const double *rdc = lp_reduced_cost(lp);
    //  const double *rdc2 = lp_reduced_cost(origLP);
    char name[256];

    /* maximum number of indexes in a var name */
#define MAX_IDX 16
    char prefix[STR_SIZE];
    int id[MAX_IDX];
    int maxT = 0;
    int maxHorizon=0;
    int ncols =lp_cols(lp);

    int nJobs = Inst_nJobs(inst);
    int nModes = Inst_nMaxModes(inst);
    // checking maximum time  to all variables //where there is some allocation ** if
    for ( int i=0 ; (i<lp_cols(lp)) ; ++i ) {
        if (fabs(xf[i])>EPS ) {
            // if (fabs(xf[i])>1e-5 || maxRC != -1.0) {
            //   if(maxRC != -1.0)
            //     if(rdc[i] > maxRC) continue;

            lp_col_name( origLP, origCols[i], name );
            if (tolower(name[0])=='x') {
                parseName( name, prefix, id );
                if (prefix[0]=='x') {
                    int j = id[0];
                    int m = id[1];
                    int t = id[2];
                    maxT = MAX( maxT, t );
                    const Job * job = Inst_job(inst,j);
                    const Mode *mode = Job_mode(job, m);
                    maxHorizon = MAX(maxHorizon,t+Mode_duration(mode)+1);
                    //      printf( "x(%d,%d,%d) xf: %g xfor %g, rdc: %f \n", j, m, t, xf[i], xf[origCols[i]], rdc[i]);
                }
            }
        }
    }
    int nTimes = maxT+1;
    *horizon = maxHorizon;
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

    // filling x from fractional solution
    for ( int i=0 ; (i<lp_cols(lp)) ; ++i ) {
        if (fabs(xf[i])>EPS ) {
            //if (fabs(xf[i])>1e-5 || maxRC != -1.0) {
            //  if(maxRC != -1.0)
            //    if(rdc[i] > maxRC) continue;

            lp_col_name( origLP, origCols[i], name );
            if (tolower(name[0]=='x')) {
                parseName( name, prefix, id );
                if (prefix[0]=='x')
                    x[id[0]][id[1]][id[2]] = xf[i];

            }
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
    for ( int i=0 ; (i<lp_cols(lp)) ; ++i ) {
        if (fabs(xf[i])>EPS ) {
            //if (fabs(xf[i])>1e-5 || maxRC != -1.0) {
            // if(maxRC != -1.0)
            //      if(rdc[i] > maxRC) continue;

            lp_col_name( origLP, origCols[i], name );
            if (tolower(name[0])=='x') {
                parseName( name, prefix, id );
                if (prefix[0]=='x') {
                    int jlp = id[0];
                    int mlp =id[1];
                    int tlp = id[2];
                    xIdx[jlp][mlp][tlp] = i;
                    //     printf("xIdx[%d][%d][%d] = %d ", jlp,mlp,tlp,i);
                }
            }
        }
    }

    //  getchar();
    //CHOSE CONSTRAINTS

    //    printf("iniciou alocacao memory\n");// getchar();
    //    int lpr = lp_rows(lp);

    VecStr *rname = VStr_create(256);
    VecStr *nameelements = VStr_create(256);
    VecDbl * rrhs = VDbl_create();
    VecInt * rtype = VInt_create();
    VecInt * rtimeslotmin = VInt_create();
    VecInt * rtimeslotmax = VInt_create();

    /*VecChar * rsense = VChar_create();
    VecInt ** ridx;
    ALLOCATE_VECTOR(ridx,VecInt*,lp_rows(lp));
    VecDbl ** rcoef;
    ALLOCATE_VECTOR(rcoef,VecDbl*,lp_rows(lp));

    for ( int r=0 ; (r<lp_rows(lp)) ; ++r ) {
        ridx[r] = VInt_create();
        rcoef[r] = VDbl_create();
    }
    */
//    CutP_create_job_set(inst,lp, maxT);

    //    int contUR=0, contURN = 0, contUM =0, contUCI = 0, contUCP = 0;
    int nJumps =0;
    int *elements;
    ALLOCATE_VECTOR_INI(elements,int,lp_cols(lp));
    int *idxelements;
    ALLOCATE_VECTOR_INI(idxelements,int,lp_cols(lp));
    VecInt **elemrow;
    ALLOCATE_VECTOR(elemrow,VecInt*,lp_cols(lp));
    VecDbl **coefelemrows;
    ALLOCATE_VECTOR(coefelemrows,VecDbl*,lp_cols(lp));
    for ( int c=0 ; (c<lp_cols(lp)) ; ++c ) {
        elemrow[c] = VInt_create();
        coefelemrows[c] = VDbl_create();
    }

    // printf("passou alocacao memory\n"); fflush(stdout); //getchar();

    int cont = 0, contElem = 0;
    int mininterval = mininstant, maxinterval = maxinstant;
    double maxvalue = 0;
    double sumvalue = 0;
    //if(maxinstant >= maxT) maxinstant = maxT;
    for ( int t=mininstant; (t<maxinstant && maxinstant <= maxT) ; ++t ) {
        for ( int j=0; (j<nJobs) ; ++j ) {
            for ( int m=0 ; (m<nModes) ; ++m ) {
                int ind = xIdx[j][m][t];
                if( ind == -1) continue;
                if (fabs(xf[ind])>EPS ) {
                    //if (fabs(xf[ind])>1e-5 || maxRC != -1.0) {
                    //  if(maxRC != -1.0)
                    //    if(rdc[ind] > maxRC) continue;

                    double intcloser = ROUND(xf[ind]);
                    double value = fabs(intcloser-xf[ind]);
                    //printf("int closer %f, xf[%d] %f, value %f\n", intcloser, ind, xf[ind], value);
                    sumvalue +=value;
                }
            }
        }
        if(t==maxinstant-1) {

            if(sumvalue >= (maxvalue)*1.25 ) {
                mininterval = mininstant;
                maxinterval = maxinstant;
                maxvalue = sumvalue;
                nJumps++;
                //   printf("sumvalue %f, maxvalue %f\n",  sumvalue, maxvalue );
            }
            //   printf("maxT %d, min interval %d, max interval %d, min instant %d, max instant %d, jump %d , sumvalue %f, maxvalue %f\n",maxT,mininterval, maxinterval, mininstant, maxinstant, jump, sumvalue, maxvalue );
            //if(VERBOSE==3)
            sumvalue = 0;
            mininstant += jump;
            maxinstant += jump;
            if(maxinstant >= maxT) maxinstant = maxT;
            if(jump !=0 )
                t = mininstant;
        }
    }

    for ( int r=0 ; (r<lp_rows(lp)) ; ++r ) {
        char nameR[256];
        lp_row_name(lp, r, nameR);
        // printf("nameR %s \n", nameR);
        // getchar();
        if(strncmp(nameR,"resR",4)== 0 || strncmp(nameR,"cutRR",5)== 0){// || strncmp(nameR,"cutE",4)== 0){
            parseName( nameR, prefix, id );
            if(id[1] < mininterval || id[1] > maxinterval ) continue;
            //printf("mininstant %d maxinstant %d %s(%d,%d)\n", mininterval, maxinterval,nameR, id[0], id[1]);
            int *idx;
            ALLOCATE_VECTOR_INI(idx,int,lp_cols(lp));
            double *coef;
            ALLOCATE_VECTOR_INI(coef,double,lp_cols(lp));
            int size = lp_row(lp,r,idx,coef);
            int aux = 0;
            if(size<=0) {
                free(idx);
                free(coef);
                continue;
            }
            for(int i  = 0 ; i < size ; i++) {
                char nameE[256];
                lp_col_name(lp,idx[i],nameE);
                parseName( nameE, prefix, id );
                int ind = xIdx[id[0]][id[1]][id[2]];
                if(ind != -1) {

                    if(elements[ind]==0 ) {
                        elements[ind] = 1;
                        idxelements[contElem] = ind;
                        VStr_pushBack(nameelements, nameE);
                        contElem++;
                    }
                    //VInt_pushBack(ridx[cont],idx[i]);
                    //VDbl_pushBack(rcoef[cont],coef[i]);
                    VInt_pushBack(elemrow[ind], cont);
                    VDbl_pushBack(coefelemrows[ind], coef[i]);

                    aux++;
                }
            }
            if(aux>0) {
                double rhs = lp_rhs(lp,r);
                VDbl_pushBack(rrhs,rhs);
                if(strncmp(nameR,"resR",4)== 0) VInt_pushBack(rtype,RES_RR);
                if(strncmp(nameR,"cutRR",5)== 0) VInt_pushBack(rtype,LPC_RR);
                //if(strncmp(nameR,"cutE",4)== 0) VInt_pushBack(rtype,LPC_ENERGETIC);

                VInt_pushBack(rtimeslotmin,mininterval);
                VInt_pushBack(rtimeslotmax,maxinterval);
                //char sense = lp_sense(lp,r);
                //VChar_pushBack(rsense,sense);
                VStr_pushBack(rname,nameR);
                cont++;
                //    printf("%s",nameR );
            }
            free(idx);
            free(coef);
        }
    }


    VecInt** conflicts;
    ALLOCATE_VECTOR(conflicts,VecInt*,ncols);


    int **varInConfs;
    ALLOCATE_VECTOR_INI(varInConfs, int*, ncols);
    int nProj = Inst_nProjects(inst);
    int *nVProj;
    ALLOCATE_VECTOR_INI(nVProj,int,nProj);
    IntTriple **vectorProjects;
    ALLOCATE_VECTOR(vectorProjects,IntTriple*,nProj);
    for(int i = 0 ; i < nProj ; i++)
        ALLOCATE_VECTOR_INI(vectorProjects[i],IntTriple,ncols);


    VecInt **sameJob;
    ALLOCATE_VECTOR_INI(sameJob,VecInt*,Inst_nJobs(inst));
    for(int i = 0; i < Inst_nJobs(inst) ; i++)
        sameJob[i] = VInt_create();
    VecInt **sameResN,**sameResNCoef;
    ALLOCATE_VECTOR_INI(sameResN,VecInt*,Inst_nResN(inst));
    ALLOCATE_VECTOR_INI(sameResNCoef,VecInt*,Inst_nResN(inst));
    for(int i = 0; i < Inst_nResN(inst); i++) {
        sameResN[i] = VInt_create();
        sameResNCoef[i] = VInt_create();
    }


    for(int l = 0; l < ncols && cont > 0; l++) {
        conflicts[l] = VInt_create();
        ALLOCATE_VECTOR_INI(varInConfs[l], int, ncols);

        if (fabs(xf[l])>EPS ) {
            //if (fabs(xf[l])>1e-5 || maxRC != -1.0) {
            //  if(maxRC != -1.0)
            //    if(rdc[l] > maxRC) continue;

            lp_col_name( lp,l, name );

            if (tolower(name[0])=='x') {
                parseName( name, prefix, id );
                if (prefix[0]=='x') {

                    int jlp = id[0];
                    int mlp =id[1];
                    int tlp = id[2];

                    int ind = xIdx[jlp][mlp][tlp];
                    if(ind == -1) continue;
                    //  printf("ind %d xIdx[%d][%d][%d];\n",ind,jlp,mlp,tlp);
                    if(elements[ind]==0) continue;

                    VInt_pushBack(sameJob[jlp], ind);

                    const Job *job = Inst_job(inst,jlp);
                    const Mode *mode = Job_mode(job,mlp);
                    if(Mode_duration(mode)==0)continue;


                    for(int mmm = 0 ; mmm < Mode_nResN(mode); mmm++) {
                        int idxResN = Mode_idxResN(mode,mmm);
                        VInt_pushBack(sameResN[idxResN], ind);
                        VInt_pushBack(sameResNCoef[idxResN], Mode_useResN(mode,mmm));
                    }

                    int plp = Job_project(job);

                    vectorProjects[plp][nVProj[plp]].idx = ind;
                    vectorProjects[plp][nVProj[plp]].j = jlp;
                    vectorProjects[plp][nVProj[plp]].m = mlp;
                    vectorProjects[plp][nVProj[plp]].t = tlp;
                    vectorProjects[plp][nVProj[plp]].value = x[jlp][mlp][tlp];
                    nVProj[plp]++;

                }
            }
        }
    }


    for(int i = 0; i < Inst_nJobs(inst) ; i++) {
        int size = VInt_size(sameJob[i]);
        if(size<=0) continue;
        int nme=0;
        // printf("\nElem for job %d\n", i);
        for(int e = 0 ; e < size ; e++ ) {
            int elem = VInt_get(sameJob[i],e);
            //    printf("elem %d ", elem);
            VInt_pushBack(elemrow[elem], cont);
            VDbl_pushBack(coefelemrows[elem], 1);
            nme++;
        }
        if(nme>0) {
            double rhs = 1;
            VDbl_pushBack(rrhs,rhs);
            VInt_pushBack(rtype,RES_MODE);
            VInt_pushBack(rtimeslotmin,mininterval);
            VInt_pushBack(rtimeslotmax,maxinterval);
            char nameCC[256];
            sprintf(nameCC,"modeSel(%d)",i);
            // printf(" %s type %d \n",nameCC,  RES_MODE);
            VStr_pushBack(rname,nameCC);
            cont++;
        }
    }

    for(int r = 0; r < Inst_nResN(inst); r++) {
        if(cont >= lp_rows(lp)) break;
        int size = VInt_size(sameResN[r]);
        if(size<=0) continue;
        int nre =0 ;
        //printf("resource %d ", r);
        for(int elem = 0 ; elem < size ; elem++) {
            int l = VInt_get(sameResN[r],elem);
            double coef= (double) VInt_get(sameResNCoef[r],elem);
            //printf("%f*%d  ", coef, l);
            VInt_pushBack(elemrow[l], cont);
            VDbl_pushBack(coefelemrows[l], coef);
            nre++;
        }
        if(nre>0) {
            double rhs = (double) Inst_capResN(inst,r);
            //   printf("rhs %f\n", rhs);
            VDbl_pushBack(rrhs,rhs);
            VInt_pushBack(rtype,RES_NR);
            VInt_pushBack(rtimeslotmin,mininterval);
            VInt_pushBack(rtimeslotmax,maxinterval);

            char nameCC[256];
            sprintf(nameCC,"resN(%d)",r);
            //printf("%s type %d",nameCC, RES_NR );
            VStr_pushBack(rname,nameCC);
            cont++;
            nre=0;
        }
    }


    for(int p = 0 ; p < nProj && cont > 0; p++) {
        //  if(cont >= lp_rows(lp)) break;
        for(int elem = 0 ; elem < nVProj[p] ; elem++) {
            //if(cont >= lp_rows(lp)) break;
            int l = vectorProjects[p][elem].idx;
            int j = vectorProjects[p][elem].j, m = vectorProjects[p][elem].m, t =  vectorProjects[p][elem].t;
            //     printf("\n E (%d,%d,%d) %d ", j, m, t, l);
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
                    //    printf("Same Job (%d,%d,%d) %d ", j2, m2, t2, l2);
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
                        //printf("\ntimeEndPred %d > t2 %d || ((t2 %d - timeEndPred %d) < winTime %d)", timeEndPred,t2,t2,timeEndPred, winTime);
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
                        //         printf("tpd p(%d,%d,%d) %d ", j2, m2, t2, l2);
                        varInConfs[l][l2] = 1;

                    }
                }
            }
        }
    }


    int **conf;
    ALLOCATE_VECTOR(conf, int*, ncols);

    int *nConf;
    ALLOCATE_VECTOR_INI(nConf, int, ncols);
    int weight[lp_cols(lp)];

    //printf("find conflicts clique: end\n");fflush(stdout);
    for(int l = 0 ; l < ncols ; l++) {
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

                //  printf("confJ[o]%d", confJ[o]);
                conf[l][o]= confJ[o];
                //printf("conf[%d][%d] %d ", l,o, conf[l][o]);
                //getchar();
            }
            CutP_quick_sort_vec(conf[l],coef,sizeConflictsJ);
            if(xf[l]<=0.00001)
                weight[l] = 1;
            else
                weight[l] = xf[l]*1000.0;
            free(coef);
        } else
            weight[l] = 0.0;
    }

    _time = ( (double) timeLeft - (omp_get_wtime()-tinit) );
    CGraph *cgraph = build_cgraph_conflicts(conf, nConf, lp_cols(lp), _time);

    cgraph_set_weight( cgraph, weight );
    BronKerbosch *bk =  bk_create(cgraph);
    _time = ( (double) timeLeft - (omp_get_wtime()-tinit) );
  //  bk_set_min_weight(bk, weight);
   // bk_set_timelimit(bk, _time);
    bk_run(bk);

    const CliqueSet *clqSet = bk_get_clq_set(bk);
    int nCliqueSet =  clq_get_number_of_cliques(clqSet);

    for(int nc = 0 ; nc < nCliqueSet ; nc++) {
        if(cont >= lp_rows(lp)) break;
        const int* iset = clq_set_clique_elements(clqSet, nc);
        int size = clq_set_clique_size(clqSet, nc);
        if(size <= 1) continue;
        double coef[size];
        FILL(coef,0,size,1);

        //    printf(" Clique: ");
        double peso = 0 ;
        int ncs = 0;
        for(int riz = 0 ; riz <size ; riz++) {
            //VInt_pushBack(ridx[cont],iset[riz]);
            //VDbl_pushBack(rcoef[cont],coef[riz]);
            char namename[256];
            lp_col_name( lp, iset[riz], namename );
            parseName( namename, prefix, id );
            //            int j = id[0];
            //            int m =id[1];
            //            int t = id[2];
            if(xf[iset[riz]]<=0.00001) {
                //        printf(" %f: ",  xf[iset[riz]]);
                peso += 1;
            } else {
                //          printf(" %f: ",  xf[iset[riz]]);
                peso += xf[iset[riz]]*1000.0;
            }
            int ind = xIdx[id[0]][id[1]][id[2]];
            if(ind != -1) {
                if(elements[ind]==0 ) {
                    elements[ind] = 1;
                    idxelements[contElem] = ind;
                    VStr_pushBack(nameelements, namename);
                    contElem++;
                }

                VInt_pushBack(elemrow[iset[riz]], cont);
                VDbl_pushBack(coefelemrows[iset[riz]], coef[riz]);
                //printf("%f * %d (%d,%d,%d) ", coef[riz], iset[riz],j,m,t);
                ncs++;
            } //else
            //printf("ind -1-----------------------------\n");
            //CutP_quick_sort_vec(elemrow[iset[riz]],coefelemrows[iset[riz]]);
        }

        if(ncs>0) {
            double rhs = 1;
            VDbl_pushBack(rrhs,rhs);
            VInt_pushBack(rtype,RES_CONFCL);
            VInt_pushBack(rtimeslotmin,mininterval);
            VInt_pushBack(rtimeslotmax,maxinterval);

            //   printf(" <= 1 peso %f\n", peso);
            // char sense = 'L';
            //VChar_pushBack(rsense,sense);
            char nameCC[256];
            sprintf(nameCC,"conflictCL(%d)",cont);
            VStr_pushBack(rname,nameCC);
            cont++;
            // printf("rhs %f %s type %d \n", rhs, nameCC, RES_CONFCL );
            // getchar();
        }

    }


    CutCG *ccg  = CutCG_create(inst, cont, contElem, ncols);

    for(int r = 0 ; r < cont  ; r++) {
        double rhs = VDbl_get(rrhs,r);
        VDbl_pushBack(ccg->rowrhs,rhs);
        VDbl_pushBack(ccg->rrhs,rhs);
        int type = VInt_get(rtype,r);
        VInt_pushBack(ccg->rowtype,type);
        int timeslotmin = VInt_get(rtimeslotmin,r);
        int timeslotmax = VInt_get(rtimeslotmax,r);
        VInt_pushBack(ccg->rowtimeslotmin,timeslotmin);
        VInt_pushBack(ccg->rowtimeslotmax,timeslotmax);
        VStr_pushBack(ccg->rname, VStr_get(rname,r));
    }

    for(int e = 0 ; e < contElem ; e++) {
        int idxel = idxelements[e];
        ccg->idxelements[e] = idxelements[e];
        VDbl_pushBack(ccg->xfElemPP,xf[idxel]);
        VInt_pushBack(ccg->idxElemPP,e);
        VInt_pushBack(ccg->idxelementsvec,idxel);

        VStr_pushBack(ccg->nameelements,VStr_get(nameelements,e));
    }


    for(int el = 0 ; el < contElem ; el++) {
        int elem = idxelements[el];
        int elemPP = VInt_get(ccg->idxElemPP,el);
        int nRows = VInt_size(elemrow[elem]);
        for(int i  = 0 ; i < nRows ; i++ ) {
            int row = VInt_get(elemrow[elem], i);
            VInt_pushBack(ccg->elemrow[elem],row);
            double coef = VDbl_get(coefelemrows[elem], i);
            VDbl_pushBack(ccg->coefelemrows[elem], coef);
            VInt_pushBack(ccg->rowElem[row],elemPP);
            VDbl_pushBack(ccg->rowCoef[row],coef);
        }
    }

    bk_free(bk);
    cgraph_free(&cgraph);
    //clq_set_free(&clqSet);

    free(elements);
    free(idxelements);
    VStr_free(&rname);
    VDbl_free(&rrhs);
    VInt_free(&rtype);
    VInt_free(&rtimeslotmin);
    VInt_free(&rtimeslotmax);
    for(int i = 0; i < Inst_nResN(inst); i++) {
        VInt_free(&sameResN[i]);
        VInt_free(&sameResNCoef[i]);
    }
    for(int j = 0 ; j < Inst_nJobs(inst) ; j++)
        VInt_free(&sameJob[j]);

    for ( int c=0 ; (c<lp_cols(lp)) ; ++c ) {
        VInt_free(&elemrow[c]);
        VInt_free(&conflicts[c]);
        VDbl_free(&coefelemrows[c]);
        if(nConf[c]>0) free(conf[c]);
    }

    free(conflicts);
    free(conf);
    free(nConf);

    free(sameJob);
    free(sameResN);
    free(sameResNCoef);
    free(coefelemrows);
    free(elemrow);

    VStr_free(&nameelements);

    for(int p = 0 ; p < nProj ; p++)
        free(vectorProjects[p]);
    free(vectorProjects);
    free(nVProj);

    for(int l = 0 ; l < lp_cols(lp) ; l++)
        free(varInConfs[l]);
    free(varInConfs);


    free(xIdx[0][0]);
    free(xIdx[0]);
    free(xIdx);
    free(x[0][0]);
    free(x[0]);
    free(x);


    return ccg;

}

/*identify the set of important rows to compose the CG procedure with all constraints (RR, cutRR, NR, Job Selection,) and conflicts about precedence relationship and project delay*/
CutCG *CutCG_create_and_identify_rows_all_conflitcts_cgraph(  LinearProgram *lp, const CGraph *cgraph, const int *origCols, LinearProgram *origLP, const Instance *inst, double timeLeft, int mininstant, int maxinstant, int jump, int nround, int *horizon)
{

    double _time;
    double tinit = omp_get_wtime();

    assert(lp);
    assert(inst);
    assert(origLP);

    //    int nround = Res_getRound(res);

    const double *xf = lp_x(lp);
    //double maxRC = lp_get_max_reduced_cost(origLP);
    //const double *rdc = lp_reduced_cost(lp);
    //  const double *rdc2 = lp_reduced_cost(origLP);
    char name[256];

    /* maximum number of indexes in a var name */
#define MAX_IDX 16
    char prefix[STR_SIZE];
    int id[MAX_IDX];
    int maxT = 0;
    int maxHorizon=0;
    int ncols =lp_cols(lp);

    int nJobs = Inst_nJobs(inst);
    int nModes = Inst_nMaxModes(inst);
    // checking maximum time  to all variables //where there is some allocation ** if
    for ( int i=0 ; (i<lp_cols(lp)) ; ++i ) {
        if (fabs(xf[i])>EPS ) {
            // if (fabs(xf[i])>1e-5 || maxRC != -1.0) {
            //   if(maxRC != -1.0)
            //     if(rdc[i] > maxRC) continue;

            lp_col_name( origLP, origCols[i], name );
            if (tolower(name[0])=='x') {
                parseName( name, prefix, id );
                if (prefix[0]=='x') {
                    int j = id[0];
                    int m = id[1];
                    int t = id[2];
                    maxT = MAX( maxT, t );
                    const Job * job = Inst_job(inst,j);
                    const Mode *mode = Job_mode(job, m);
                    maxHorizon = MAX(maxHorizon,t+Mode_duration(mode)+1);
                    //      printf( "x(%d,%d,%d) xf: %g xfor %g, rdc: %f \n", j, m, t, xf[i], xf[origCols[i]], rdc[i]);
                }
            }
        }
    }
    int nTimes = maxT+1;
    *horizon = maxHorizon;
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

    // filling x from fractional solution
    for ( int i=0 ; (i<lp_cols(lp)) ; ++i ) {
        if (fabs(xf[i])>EPS ) {
            //if (fabs(xf[i])>1e-5 || maxRC != -1.0) {
            //  if(maxRC != -1.0)
            //    if(rdc[i] > maxRC) continue;

            lp_col_name( origLP, origCols[i], name );
            if (tolower(name[0]=='x')) {
                parseName( name, prefix, id );
                if (prefix[0]=='x')
                    x[id[0]][id[1]][id[2]] = xf[i];

            }
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
    for ( int i=0 ; (i<lp_cols(lp)) ; ++i ) {
        if (fabs(xf[i])>EPS ) {
            //if (fabs(xf[i])>1e-5 || maxRC != -1.0) {
            // if(maxRC != -1.0)
            //      if(rdc[i] > maxRC) continue;

            lp_col_name( origLP, origCols[i], name );
            if (tolower(name[0])=='x') {
                parseName( name, prefix, id );
                if (prefix[0]=='x') {
                    int jlp = id[0];
                    int mlp =id[1];
                    int tlp = id[2];
                    xIdx[jlp][mlp][tlp] = i;
                    //     printf("xIdx[%d][%d][%d] = %d ", jlp,mlp,tlp,i);
                }
            }
        }
    }

    //  getchar();
    //CHOSE CONSTRAINTS

    //    printf("iniciou alocacao memory\n");// getchar();
    //    int lpr = lp_rows(lp);

    VecStr *rname = VStr_create(256);
    VecStr *nameelements = VStr_create(256);
    VecDbl * rrhs = VDbl_create();
    VecInt * rtype = VInt_create();
    VecInt * rtimeslotmin = VInt_create();
    VecInt * rtimeslotmax = VInt_create();

    /*VecChar * rsense = VChar_create();
    VecInt ** ridx;
    ALLOCATE_VECTOR(ridx,VecInt*,lp_rows(lp));
    VecDbl ** rcoef;
    ALLOCATE_VECTOR(rcoef,VecDbl*,lp_rows(lp));

    for ( int r=0 ; (r<lp_rows(lp)) ; ++r ) {
        ridx[r] = VInt_create();
        rcoef[r] = VDbl_create();
    }
    */
//    CutP_create_job_set(inst,lp, maxT);

    //    int contUR=0, contURN = 0, contUM =0, contUCI = 0, contUCP = 0;
    int nJumps =0;
    int *elements;
    ALLOCATE_VECTOR_INI(elements,int,lp_cols(lp));
    int *idxelements;
    ALLOCATE_VECTOR_INI(idxelements,int,lp_cols(lp));
    VecInt **elemrow;
    ALLOCATE_VECTOR(elemrow,VecInt*,lp_cols(lp));
    VecDbl **coefelemrows;
    ALLOCATE_VECTOR(coefelemrows,VecDbl*,lp_cols(lp));
    for ( int c=0 ; (c<lp_cols(lp)) ; ++c ) {
        elemrow[c] = VInt_create();
        coefelemrows[c] = VDbl_create();
    }

    // printf("passou alocacao memory\n"); fflush(stdout); //getchar();

    int cont = 0, contElem = 0;
    int mininterval = mininstant, maxinterval = maxinstant;
    double maxvalue = 0;
    double sumvalue = 0;
    //if(maxinstant >= maxT) maxinstant = maxT;
    for ( int t=mininstant; (t<maxinstant && maxinstant <= maxT) ; ++t ) {
        for ( int j=0; (j<nJobs) ; ++j ) {
            for ( int m=0 ; (m<nModes) ; ++m ) {
                int ind = xIdx[j][m][t];
                if( ind == -1) continue;
                if (fabs(xf[ind])>EPS ) {
                    //if (fabs(xf[ind])>1e-5 || maxRC != -1.0) {
                    //  if(maxRC != -1.0)
                    //    if(rdc[ind] > maxRC) continue;

                    double intcloser = ROUND(xf[ind]);
                    double value = fabs(intcloser-xf[ind]);
                    //printf("int closer %f, xf[%d] %f, value %f\n", intcloser, ind, xf[ind], value);
                    sumvalue +=value;
                }
            }
        }
        if(t==maxinstant-1) {

            if(sumvalue >= (maxvalue)*1.25 ) {
                mininterval = mininstant;
                maxinterval = maxinstant;
                maxvalue = sumvalue;
                nJumps++;
                //   printf("sumvalue %f, maxvalue %f\n",  sumvalue, maxvalue );
            }
            //   printf("maxT %d, min interval %d, max interval %d, min instant %d, max instant %d, jump %d , sumvalue %f, maxvalue %f\n",maxT,mininterval, maxinterval, mininstant, maxinstant, jump, sumvalue, maxvalue );
            //if(VERBOSE==3)
            sumvalue = 0;
            mininstant += jump;
            maxinstant += jump;
            if(maxinstant >= maxT) maxinstant = maxT;
            if(jump !=0 )
                t = mininstant;
        }
    }

    for ( int r=0 ; (r<lp_rows(lp)) ; ++r ) {
        char nameR[256];
        lp_row_name(lp, r, nameR);
        // printf("nameR %s \n", nameR);
        // getchar();
        if(strncmp(nameR,"resR",4)== 0 || strncmp(nameR,"cutRR",5)== 0){// || strncmp(nameR,"cutE",4)== 0){
            parseName( nameR, prefix, id );
            if(id[1] < mininterval || id[1] > maxinterval ) continue;
            //printf("mininstant %d maxinstant %d %s(%d,%d)\n", mininterval, maxinterval,nameR, id[0], id[1]);
            int *idx;
            ALLOCATE_VECTOR_INI(idx,int,lp_cols(lp));
            double *coef;
            ALLOCATE_VECTOR_INI(coef,double,lp_cols(lp));
            int size = lp_row(lp,r,idx,coef);
            int aux = 0;
            if(size<=0) {
                free(idx);
                free(coef);
                continue;
            }
            for(int i  = 0 ; i < size ; i++) {
                char nameE[256];
                lp_col_name(lp,idx[i],nameE);
                parseName( nameE, prefix, id );
                int ind = xIdx[id[0]][id[1]][id[2]];
                if(ind != -1) {

                    if(elements[ind]==0 ) {
                        elements[ind] = 1;
                        idxelements[contElem] = ind;
                        VStr_pushBack(nameelements, nameE);
                        contElem++;
                    }
                    //VInt_pushBack(ridx[cont],idx[i]);
                    //VDbl_pushBack(rcoef[cont],coef[i]);
                    VInt_pushBack(elemrow[ind], cont);
                    VDbl_pushBack(coefelemrows[ind], coef[i]);

                    aux++;
                }
            }
            if(aux>0) {
                double rhs = lp_rhs(lp,r);
                VDbl_pushBack(rrhs,rhs);
                if(strncmp(nameR,"resR",4)== 0) VInt_pushBack(rtype,RES_RR);
                if(strncmp(nameR,"cutRR",5)== 0) VInt_pushBack(rtype,LPC_RR);
                //if(strncmp(nameR,"cutE",4)== 0) VInt_pushBack(rtype,LPC_ENERGETIC);

                VInt_pushBack(rtimeslotmin,mininterval);
                VInt_pushBack(rtimeslotmax,maxinterval);
                //char sense = lp_sense(lp,r);
                //VChar_pushBack(rsense,sense);
                VStr_pushBack(rname,nameR);
                cont++;
                //    printf("%s",nameR );
            }
            free(idx);
            free(coef);
        }
    }


    int **varInConfs;
    ALLOCATE_VECTOR_INI(varInConfs, int*, ncols);
    int nProj = Inst_nProjects(inst);
    int *nVProj;
    ALLOCATE_VECTOR_INI(nVProj,int,nProj);
    IntTriple **vectorProjects;
    ALLOCATE_VECTOR(vectorProjects,IntTriple*,nProj);
    for(int i = 0 ; i < nProj ; i++)
        ALLOCATE_VECTOR_INI(vectorProjects[i],IntTriple,ncols);


    VecInt **sameJob;
    ALLOCATE_VECTOR_INI(sameJob,VecInt*,Inst_nJobs(inst));
    for(int i = 0; i < Inst_nJobs(inst) ; i++)
        sameJob[i] = VInt_create();
    VecInt **sameResN,**sameResNCoef;
    ALLOCATE_VECTOR_INI(sameResN,VecInt*,Inst_nResN(inst));
    ALLOCATE_VECTOR_INI(sameResNCoef,VecInt*,Inst_nResN(inst));
    for(int i = 0; i < Inst_nResN(inst); i++) {
        sameResN[i] = VInt_create();
        sameResNCoef[i] = VInt_create();
    }


    double *xfcghaph;
    ALLOCATE_VECTOR_INI(xfcghaph, double, ncols);

    for(int l = 0; l < ncols && cont > 0; l++) {
        ALLOCATE_VECTOR_INI(varInConfs[l], int, ncols);

        if (fabs(xf[l])>EPS ) {
            //if (fabs(xf[l])>1e-5 || maxRC != -1.0) {
            //  if(maxRC != -1.0)
            //    if(rdc[l] > maxRC) continue;

            lp_col_name( lp,l, name );

            if (tolower(name[0])=='x') {
                parseName( name, prefix, id );
                if (prefix[0]=='x') {

                    int jlp = id[0];
                    int mlp =id[1];
                    int tlp = id[2];

                    int ind = xIdx[jlp][mlp][tlp];
                    if(ind == -1) continue;
                    //  printf("ind %d xIdx[%d][%d][%d];\n",ind,jlp,mlp,tlp);
                    if(elements[ind]==0) continue;

                    VInt_pushBack(sameJob[jlp], ind);

                    const Job *job = Inst_job(inst,jlp);
                    const Mode *mode = Job_mode(job,mlp);
                    if(Mode_duration(mode)==0)continue;


                    for(int mmm = 0 ; mmm < Mode_nResN(mode); mmm++) {
                        int idxResN = Mode_idxResN(mode,mmm);
                        VInt_pushBack(sameResN[idxResN], ind);
                        VInt_pushBack(sameResNCoef[idxResN], Mode_useResN(mode,mmm));
                    }

                    int plp = Job_project(job);

                    vectorProjects[plp][nVProj[plp]].idx = ind;
                    vectorProjects[plp][nVProj[plp]].j = jlp;
                    vectorProjects[plp][nVProj[plp]].m = mlp;
                    vectorProjects[plp][nVProj[plp]].t = tlp;
                    vectorProjects[plp][nVProj[plp]].value = x[jlp][mlp][tlp];
                    nVProj[plp]++;

                }
            }
        }
    }


    for(int i = 0; i < Inst_nJobs(inst) ; i++) {
        int size = VInt_size(sameJob[i]);
        if(size<=0) continue;
        int nme=0;
        // printf("\nElem for job %d\n", i);
        for(int e = 0 ; e < size ; e++ ) {
            int elem = VInt_get(sameJob[i],e);
            //    printf("elem %d ", elem);
            VInt_pushBack(elemrow[elem], cont);
            VDbl_pushBack(coefelemrows[elem], 1);
            nme++;
        }
        if(nme>0) {
            double rhs = 1;
            VDbl_pushBack(rrhs,rhs);
            VInt_pushBack(rtype,RES_MODE);
            VInt_pushBack(rtimeslotmin,mininterval);
            VInt_pushBack(rtimeslotmax,maxinterval);
            char nameCC[256];
            sprintf(nameCC,"modeSel(%d)",i);
            // printf(" %s type %d \n",nameCC,  RES_MODE);
            VStr_pushBack(rname,nameCC);
            cont++;
        }
    }

    for(int r = 0; r < Inst_nResN(inst); r++) {
        if(cont >= lp_rows(lp)) break;
        int size = VInt_size(sameResN[r]);
        if(size<=0) continue;
        int nre =0 ;
        //printf("resource %d ", r);
        for(int elem = 0 ; elem < size ; elem++) {
            int l = VInt_get(sameResN[r],elem);
            double coef= (double) VInt_get(sameResNCoef[r],elem);
            //printf("%f*%d  ", coef, l);
            VInt_pushBack(elemrow[l], cont);
            VDbl_pushBack(coefelemrows[l], coef);
            nre++;
        }
        if(nre>0) {
            double rhs = (double) Inst_capResN(inst,r);
            //   printf("rhs %f\n", rhs);
            VDbl_pushBack(rrhs,rhs);
            VInt_pushBack(rtype,RES_NR);
            VInt_pushBack(rtimeslotmin,mininterval);
            VInt_pushBack(rtimeslotmax,maxinterval);

            char nameCC[256];
            sprintf(nameCC,"resN(%d)",r);
            //printf("%s type %d",nameCC, RES_NR );
            VStr_pushBack(rname,nameCC);
            cont++;
            nre=0;
        }
    }


    for(int p = 0 ; p < nProj && cont > 0; p++) {
        //  if(cont >= lp_rows(lp)) break;
        for(int elem = 0 ; elem < nVProj[p] ; elem++) {
            //if(cont >= lp_rows(lp)) break;
            int l = vectorProjects[p][elem].idx;
            int j = vectorProjects[p][elem].j, m = vectorProjects[p][elem].m, t =  vectorProjects[p][elem].t;
            //     printf("\n E (%d,%d,%d) %d ", j, m, t, l);
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
                    xfcghaph[l] = vectorProjects[p][elem].value;
                    xfcghaph[l2] = vectorProjects[p][elem2].value;
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
                        //printf("\ntimeEndPred %d > t2 %d || ((t2 %d - timeEndPred %d) < winTime %d)", timeEndPred,t2,t2,timeEndPred, winTime);
                        int l2 = vectorProjects[p][elem2].idx;
                        if(varInConfs[l][l2]==1) continue;
                        //     printf("prec D(%d,%d,%d) %d ", j2, m2, t2, l2);
                        xfcghaph[l] = vectorProjects[p][elem].value;
                        xfcghaph[l2] = vectorProjects[p][elem2].value;
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
                        //VInt_pushBack(conflicts[l],l2);
                        xfcghaph[l] = vectorProjects[p][elem].value;
                        xfcghaph[l2] = vectorProjects[p2][elem2].value;
                        //         printf("tpd p(%d,%d,%d) %d ", j2, m2, t2, l2);
                        varInConfs[l][l2] = 1;

                    }
                }
            }
        }
    }


    _time = ( (double) timeLeft - (omp_get_wtime()-tinit) );
    CliqueSeparation *sep = clq_sep_create(cgraph);
    clq_sep_set_maxcliques(sep,  100);
    clq_sep_set_max_it_bk(sep, 1000);
    clq_sep_set_extend_method(sep,0);


    _time = ( (double) timeLeft - (omp_get_wtime()-tinit) );
    if(_time <= clq_sep_get_max_time_bk(sep))
        clq_sep_set_max_time_bk(sep,_time);

    clq_sep_separate( sep, xfcghaph);

    const CliqueSet *clqSet = clq_sep_get_cliques( sep );
    int nCliques = clq_get_number_of_cliques(clqSet);

     for(int nc = 0 ; nc < nCliques ; nc++) {
        if(cont >= lp_rows(lp)) break;
        const int* iset = clq_set_clique_elements(clqSet, nc);
        int size = clq_set_clique_size(clqSet, nc);
        if(size <= 1) continue;
        double coef[size];
        FILL(coef,0,size,1);

        //    printf(" Clique: ");
        int ncs = 0;
        for(int riz = 0 ; riz <size ; riz++) {
            char namename[256];
            lp_col_name( lp, iset[riz], namename );
            parseName( namename, prefix, id );
            int ind = xIdx[id[0]][id[1]][id[2]];
            if(ind != -1) {
                if(elements[ind]==0 ) {
                    elements[ind] = 1;
                    idxelements[contElem] = ind;
                    VStr_pushBack(nameelements, namename);
                    contElem++;
                }

                VInt_pushBack(elemrow[iset[riz]], cont);
                VDbl_pushBack(coefelemrows[iset[riz]], coef[riz]);
                //printf("%f * %d (%d,%d,%d) ", coef[riz], iset[riz],j,m,t);
                ncs++;
            } //else
            //printf("ind -1-----------------------------\n");
            //CutP_quick_sort_vec(elemrow[iset[riz]],coefelemrows[iset[riz]]);
        }

        if(ncs>0) {
            double rhs = 1;
            VDbl_pushBack(rrhs,rhs);
            VInt_pushBack(rtype,RES_CONFCL);
            VInt_pushBack(rtimeslotmin,mininterval);
            VInt_pushBack(rtimeslotmax,maxinterval);
            // char sense = 'L';
            //VChar_pushBack(rsense,sense);
            char nameCC[256];
            sprintf(nameCC,"conflictCL(%d)",cont);
            VStr_pushBack(rname,nameCC);
            cont++;
            // printf("rhs %f %s type %d \n", rhs, nameCC, RES_CONFCL );
            // getchar();
        }

    }


    CutCG *ccg  = CutCG_create(inst, cont, contElem, ncols);

    for(int r = 0 ; r < cont  ; r++) {
        double rhs = VDbl_get(rrhs,r);
        VDbl_pushBack(ccg->rowrhs,rhs);
        VDbl_pushBack(ccg->rrhs,rhs);
        int type = VInt_get(rtype,r);
        VInt_pushBack(ccg->rowtype,type);
        int timeslotmin = VInt_get(rtimeslotmin,r);
        int timeslotmax = VInt_get(rtimeslotmax,r);
        VInt_pushBack(ccg->rowtimeslotmin,timeslotmin);
        VInt_pushBack(ccg->rowtimeslotmax,timeslotmax);
        VStr_pushBack(ccg->rname, VStr_get(rname,r));
    }

    for(int e = 0 ; e < contElem ; e++) {
        int idxel = idxelements[e];
        ccg->idxelements[e] = idxelements[e];
        VDbl_pushBack(ccg->xfElemPP,xf[idxel]);
        VInt_pushBack(ccg->idxElemPP,e);
        VInt_pushBack(ccg->idxelementsvec,idxel);

        VStr_pushBack(ccg->nameelements,VStr_get(nameelements,e));
    }


    for(int el = 0 ; el < contElem ; el++) {
        int elem = idxelements[el];
        int elemPP = VInt_get(ccg->idxElemPP,el);
        int nRows = VInt_size(elemrow[elem]);
        for(int i  = 0 ; i < nRows ; i++ ) {
            int row = VInt_get(elemrow[elem], i);
            VInt_pushBack(ccg->elemrow[elem],row);
            double coef = VDbl_get(coefelemrows[elem], i);
            VDbl_pushBack(ccg->coefelemrows[elem], coef);
            VInt_pushBack(ccg->rowElem[row],elemPP);
            VDbl_pushBack(ccg->rowCoef[row],coef);
        }
    }

    clq_sep_free(&sep);
    //clq_set_free(&clqSet);

    free(xfcghaph);
    free(elements);
    free(idxelements);
    VStr_free(&rname);
    VDbl_free(&rrhs);
    VInt_free(&rtype);
    VInt_free(&rtimeslotmin);
    VInt_free(&rtimeslotmax);
    for(int i = 0; i < Inst_nResN(inst); i++) {
        VInt_free(&sameResN[i]);
        VInt_free(&sameResNCoef[i]);
    }
    for(int j = 0 ; j < Inst_nJobs(inst) ; j++)
        VInt_free(&sameJob[j]);

    for ( int c=0 ; (c<lp_cols(lp)) ; ++c ) {
        VInt_free(&elemrow[c]);
        VDbl_free(&coefelemrows[c]);
    }

    free(sameJob);
    free(sameResN);
    free(sameResNCoef);
    free(coefelemrows);
    free(elemrow);

    VStr_free(&nameelements);

    for(int p = 0 ; p < nProj ; p++)
        free(vectorProjects[p]);
    free(vectorProjects);
    free(nVProj);

    for(int l = 0 ; l < lp_cols(lp) ; l++)
        free(varInConfs[l]);
    free(varInConfs);


    free(xIdx[0][0]);
    free(xIdx[0]);
    free(xIdx);
    free(x[0][0]);
    free(x[0]);
    free(x);


    return ccg;

}


/* create a MIP to separate CG cuts and add it to the cutPool, also it applies the lift procedure over the generated cuts*/
void CutCG_add_cuts_model_CG_parallel( CutCG *ccg, const CGraph *cgraph, int cutType, LinearProgram *lp, const Instance *inst, double timeLeft,   double maxcuts, int nround, int horizon)
{

    double _time;
    double startT = omp_get_wtime();

    assert(lp);
    assert(inst);
    const double *xf = lp_x(lp);
    char prefix[STR_SIZE];
    int id[MAX_IDX];


    // int contUR=0, contURN = 0, contUM =0, contUCI = 0, contUCP = 0, nJumps =0;

    //CREATE MIP CHAVATAL
    //  printf("CREATE MIP CHAVATAL \n");
    _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
    //    printf("time %f", _time);fflush(stdout);

    LinearProgram *mipCGSep = lp_create_cgsep(xf, ccg->cont, ccg->contElem, ccg->idxelements, ccg->nameelements, ccg->elemrow, ccg->coefelemrows, ccg->rrhs, ccg->rname, _time);
    //  printf("Starting optimize \n");
    // lp_write_lp(mipCGSep,"cg.lp");



    lp_set_numeric_focus(mipCGSep, 3);
    _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
    //printf("time %f", _time);fflush(stdout);
    lp_set_max_seconds(mipCGSep,(int)_time);

    //lp_set_concurrentMIP(mipCGSep,1);
    //lp_set_method(mipCGSep,4);
    //lp_set_seed(mipCGSep,100000);
    int st = lp_optimize(mipCGSep);
    // printf("Status st %d \n", st);

    //  lp_write_sol(mipCGSep,"cg.sol");
    //  getchar();

    //CREAT CUT
    int nSol = lp_num_saved_sols(mipCGSep);
    if(VERBOSE==3)  printf("\nnSol %d \n", nSol);


     ALLOCATE_VECTOR(ccg->cutElem,VecInt*,nSol);
     ALLOCATE_VECTOR(ccg->cutCoef,VecDbl*,nSol);


   // REALLOCATE_VECTOR(ccg->cutElem,VecInt*,nSol);
   // REALLOCATE_VECTOR(ccg->cutCoef,VecDbl*,nSol);

    for(int s = 0 ; s <nSol && (st==1 || st ==3 || st ==5 ); s++) {
        //if(s>0) {
            ccg->cutElem[s] = VInt_create();
            ccg->cutCoef[s] = VDbl_create();
            ccg->nAllocations++;
        //}
        if(VERBOSE==3) printf("\nsol %d obj %f\n", s, lp_saved_sol_obj(mipCGSep,s));
        const double *af = lp_saved_sol_x(mipCGSep, s);
        double sum = 0.0;
        int nc = lp_cols(mipCGSep);
        //   int idx[nc];
        //  double coef[nc];
        double rhs = 0;
        int cc = 0;
        for ( int i=0 ; (i<nc) ; ++i ) {
            if(fabs(af[i]) <= EPS  ) continue;
            char nameC[STR_SIZE];
            lp_col_name( mipCGSep, i, nameC );
            // printf("E %d: %s %f\n",i, nameC, af[i]);

            if(strncmp(nameC,"aRHS",4)== 0)
                rhs = af[i];
            else if(strncmp(nameC,"a",1)== 0) {
                parseName( nameC, prefix, id );
                char nameS[STR_SIZE];
                sprintf(nameS,"x(%d,%d,%d)", id[0], id[1], id[2]);

                int idsep = lp_col_index(lp,nameS);

                const int elem = idsep;
                //idx[cc] = elem;
                //coef[cc] = af[i];
                VInt_pushBack(ccg->cutElem[s], elem);
                VDbl_pushBack(ccg->cutCoef[s], af[i]);
                cc++;
                sum+= af[i]*xf[elem];
            }
            if(strncmp(nameC,"uresR",5)== 0) {
                ccg->contUR++;
                if(VERBOSE==3)   printf("contUR %d\n", ccg->contUR);
            }
            if(strncmp(nameC,"uresN",5)== 0) {
                ccg->contURN++;
                if(VERBOSE==3)printf("contURN %d\n", ccg->contURN);
            }
            if(strncmp(nameC,"umodeSel",8)== 0) {
                ccg->contUM++;
                if(VERBOSE==3) printf("contUM %d\n", ccg->contUM);
            }
            if(strncmp(nameC,"uconflictWT",11)== 0) {
                ccg->contUCI++;
                if(VERBOSE==3) printf("contUCI %d\n", ccg->contUCI);
            }
            if(strncmp(nameC,"uconflictP",10)== 0) {
                ccg->contUCP++;
                if(VERBOSE==3) printf("contUCP %d\n", ccg->contUCP);
            }

        }
        if(VERBOSE==3) printf("\n\n");
        //getchar();
        if(cc==0) continue;

        char nameC[STR_SIZE];
        sprintf(nameC,"CG(%d)#%d",nround, lp_rows(lp)+s+1);

        int *cutIdx = VInt_getPtr(ccg->cutElem[s]);
        double *cutCoef = VDbl_getPtr(ccg->cutCoef[s]);
        VStr_pushBack(ccg->cutname, nameC);
        VDbl_pushBack(ccg->cutrhs, rhs);
        VInt_pushBack(ccg->cutsense, 0);
        VInt_pushBack(ccg->cutdominated,0);
        VInt_pushBack(ccg->cutnelem, cc);
        VDbl_pushBack(ccg->cutviolation, sum-rhs);

        //   printf("\nSENSE %d RHS %f VI %f\n", 2,0.0, VDbl_get(ccg->cutviolation,icc));fflush(stdout);
        //   getchar();
        /*for(int i = 0 ; i < cc ; i++){
            printf(" %f*%d ",VDbl_get(ccg->cutCoef[s],i),VInt_get(ccg->cutElem[s],i));
        }
        printf("\nEnd before \n");*/
        CutP_quick_sort_vec(cutIdx,cutCoef, cc);
        _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
        //double newrhs = CutP_model_lift_cgraph( inst,cgraph, VStr_get(ccg->cutname,s),VInt_get(ccg->cutnelem,s),VInt_getPtr(ccg->cutElem[s]), VDbl_getPtr(ccg->cutCoef[s]), VDbl_get(ccg->cutrhs,s), lp,_time, horizon);
        double newrhs = CutP_model_lift( inst, VStr_get(ccg->cutname,s),VInt_get(ccg->cutnelem,s),VInt_getPtr(ccg->cutElem[s]), VDbl_getPtr(ccg->cutCoef[s]), VDbl_get(ccg->cutrhs,s), lp,_time, horizon);
       // double newrhs = CutP_model_lift_basic( inst, VStr_get(ccg->cutname,s),VInt_get(ccg->cutnelem,s),VInt_getPtr(ccg->cutElem[s]), VDbl_getPtr(ccg->cutCoef[s]), VDbl_get(ccg->cutrhs,s), lp,_time, horizon);



        if(newrhs!=-1)
            VDbl_set(ccg->cutrhs,s, newrhs);

        /*for(int i = 0 ; i < cc ; i++){
            printf(" %f*%d ",VDbl_get(ccg->cutCoef[s],i),VInt_get(ccg->cutElem[s],i));
        }
        printf("\nEnd after sort \n");
        fflush(stdout);*/
        //    if( sum <= rhs+0.00001  || cc==0)
        //          continue;

        //        if(VERBOSE==3)   printf("\n sum %f <=rhs %f : %f \n", sum, rhs,  sum-rhs);

    }

    lp_free(&mipCGSep);
}

/*add cuts to cutPool and translate the elements id from the procedure of combinatorial CG CPU or GPU to the original elements,
also it applies the lift procedure and dominance procedure between cuts generated
*/
void CutCG_add_cuts_CG_parallel( CutCG *ccg, const CGraph *cgraph, int cutType,LinearProgram *lp, const Instance *inst, double timeLeft, int continuous,  double maxcuts, int nround, int horizon)
{

    assert(ccg);
    assert(lp);
    assert(inst);
    double init = omp_get_wtime();

    //   char prefix[STR_SIZE];ĵĵĵ
    //    int id[MAX_IDX];

    //  REALLOCATE_VECTOR(ccggpu->cutElem,VecInt*,VDbl_size(ccggpu->cutrhs));
    // REALLOCATE_VECTOR(ccggpu->cutCoef,VecDbl*,VDbl_size(ccggpu->cutrhs));

    int ncut2 = 0;
    //printf("num cuts %d\n", VDbl_size(ccg->cutrhs));
    for( int r = 0 ; r < VDbl_size(ccg->cutrhs) ; r++) {
        //   if(r>0){
        //     ccggpu->cutElem[r] = VInt_create();
        //    ccggpu->cutCoef[r] = VDbl_create();
        // }
        double value = 0;
        int sizeElem = VDbl_size(ccg->cutCoef[r]);
        // int idx[sizeElem];
        // double coef[sizeElem];
        int cc = 0;
        for( int el = 0 ; el < sizeElem ; el++) {
            int elemPP = VInt_get(ccg->cutElem[r],el);
            double coefv = VDbl_get(ccg->cutCoef[r],el);
            int  elem = ccg->idxelements[elemPP];
            VInt_set(ccg->cutElem[r],el,elem);
            VDbl_set(ccg->cutCoef[r],el,coefv);
            //  idx[el] = elem;
            // coef[el] = coefv;
            value += coefv * VDbl_get(ccg->xfElemPP,elemPP);
            //    printf(" %f * %d (x* = %f) ", VDbl_get(ccg->cutCoef[r],el),elemPP, VDbl_get(ccg->xfElemPP,elemPP));
            cc++;
        }
        if( cc==0)
            continue;

        //double rhs = VDbl_get(ccggpu->cutrhs,r);
        //  printf(" rhs %f \n",VDbl_get(ccg->cutrhs,r));
        // printf(" %f < %f\n", value, VDbl_get(ccg->cutrhs,r));
        char nameC[STR_SIZE];
        if(cutType==LPC_CGCPU)
            sprintf(nameC,"cutCGCPU(%d)#%d", nround, lp_rows(lp)+r+1);
       if(cutType==LPC_CGGPUR2)
            sprintf(nameC,"CGGPUR2(%d)#%d", nround, lp_rows(lp)+r+1);

        int *cutIdx = VInt_getPtr(ccg->cutElem[r]);
        double *cutCoef = VDbl_getPtr(ccg->cutCoef[r]);
        VStr_pushBack(ccg->cutname, nameC);
        /* já está calculado dentro do método da GPU
         VInt_set(ccggpu->cutsense,r, 0);
         VInt_set(ccggpu->cutdominated,r, 0);
         VInt_set(ccggpu->cutnelem,r, cc);
         VDbl_set(ccggpu->cutviolation,r,value-rhs);
        */
        ncut2++;
        CutP_quick_sort_vec(cutIdx,cutCoef, cc);
        double _time = timeLeft -(omp_get_wtime()-init);
        double newrhs = CutP_model_lift( inst, VStr_get(ccg->cutname,r),VInt_get(ccg->cutnelem,r),VInt_getPtr(ccg->cutElem[r]), VDbl_getPtr(ccg->cutCoef[r]), VDbl_get(ccg->cutrhs,r), lp,_time, horizon);

        if(newrhs!=-1)
            VDbl_set(ccg->cutrhs,r, newrhs);

        //   if( value <= rhs+0.00001  || cc==0)
        //      continue;
        // if(VERBOSE==3) printf("\n value %f <= rhs %f : %f \n", value, rhs,  value-rhs);

        //   printf("\nSENSE %d RHS %f VI %f\n", 2,0.0, VDbl_get(ccggpu->cutviolation,icc));fflush(stdout);
        //   getchar();
    }

    //    double timeinitdomresource = omp_get_wtime();
    for(int ncdA = 0; ncdA < ncut2 ; ncdA++) {
        for(int ncdB = ncdA+1; ncdB < ncut2 ; ncdB++) {
            if(VInt_get(ccg->cutdominated,ncdA) == 1) break;
            if(VInt_get(ccg->cutdominated,ncdB) == 1) continue;
            CutP_dominanceBetweenTwo( lp, ncdA, VInt_getPtr(ccg->cutElem[ncdA]), VDbl_getPtr(ccg->cutCoef[ncdA]), VDbl_get(ccg->cutrhs,ncdA), VInt_get(ccg->cutsense,ncdA), VInt_get(ccg->cutnelem,ncdA), ncdB, VInt_getPtr(ccg->cutElem[ncdB]), VDbl_getPtr(ccg->cutCoef[ncdB]), VDbl_get(ccg->cutrhs,ncdB), VInt_get(ccg->cutsense,ncdB), VInt_get(ccg->cutnelem,ncdB),  ccg->cutdominated);
        }
    }
    //    printf("time compute dominance cg %f\n",omp_get_wtime()-timeinitdomresource);

}


/*get and print functions*/

VecDbl *CutCG_getrowrhs(CutCG *ccg)
{
    return ccg->rowrhs;
}

void CutCG_print(CutCG *ccg)
{
    printf("\nRows: \n");
    for( int r = 0 ; r < VDbl_size(ccg->rowrhs) ; r++) {
        //if(VDbl_size(ccg->rowCoef[r])<=0) continue;
        for( int el = 0 ; el < VDbl_size(ccg->rowCoef[r]) ; el++) {
            int elem = VInt_get(ccg->rowElem[r],el);
            printf(" %f * %d (x* = %f) ", VDbl_get(ccg->rowCoef[r],el),elem, VDbl_get(ccg->xfElemPP,elem));
        }
        printf(" rhs %f \n",VDbl_get(ccg->rowrhs,r));
        getchar();
    }
}

void CutCG_printCut(CutCG *ccg)
{


    printf("num cuts %d\n", VDbl_size(ccg->cutrhs));
    for( int r = 0 ; r < VDbl_size(ccg->cutrhs) ; r++) {
        double value = 0;
        for( int el = 0 ; el < VDbl_size(ccg->cutCoef[r]) ; el++) {
            int elem = VInt_get(ccg->cutElem[r],el);
            value += VDbl_get(ccg->cutCoef[r],el) * VDbl_get(ccg->xfElemPP,elem);
            printf(" %f * %d (x* = %f) ", VDbl_get(ccg->cutCoef[r],el),elem, VDbl_get(ccg->xfElemPP,elem));
        }
        printf(" rhs %f \n",VDbl_get(ccg->cutrhs,r));
        printf(" %f < %f \n", value, VDbl_get(ccg->cutrhs,r));
        //            getchar();
    }

}

VecDbl *CutCG_getxfElemPP(CutCG *ccg)
{
    return ccg->xfElemPP;
}

/*free memory*/
void CutCG_free( CutCG **_cutCG )
{

    CutCG *cutCG = *_cutCG;
    free(cutCG->elements);
    free(cutCG->idxelements);
    VInt_free(&cutCG->idxelementsvec);
    VDbl_free(&cutCG->xfElemPP);

    VStr_free(&cutCG->rname);
    //    VChar_free(&rsense);
    VDbl_free(&cutCG->rrhs);
    VInt_free(&cutCG->idxElemPP);

    for ( int c=0 ; (c<cutCG->nCols) ; ++c ) {
        VInt_free(&cutCG->elemrow[c]);
        VDbl_free(&cutCG->coefelemrows[c]);

    }
    free(cutCG->coefelemrows);
    free(cutCG->elemrow);

    for(int nr = 0 ; nr < VInt_size(cutCG->rowtype); nr++) {
        VInt_free(&cutCG->rowElem[nr]);
        VDbl_free(&cutCG->rowCoef[nr]);
    }

    free(cutCG->rowElem);
    free(cutCG->rowCoef);
    VDbl_free(&cutCG->rowrhs);
    VInt_free(&cutCG->rowtype);
    VInt_free(&cutCG->rowtimeslotmin);
    VInt_free(&cutCG->rowtimeslotmax);
    VStr_free(&cutCG->nameelements);

    //for(int nc = 0; nc < VInt_size(cutCG->cutsense); nc++) {
    for(int nc = 0; nc < cutCG->nAllocations; nc++) {
        VInt_free(&cutCG->cutElem[nc]);
        VDbl_free(&cutCG->cutCoef[nc]);
    }

    cutCG->nAllocations = 0;

    VDbl_free(&cutCG->cutrhs);
    VDbl_free(&cutCG->cutviolation);
    VInt_free(&cutCG->cutnelem);
    VInt_free(&cutCG->cutsense);
    VInt_free(&cutCG->cutdominated);
    VStr_free(&cutCG->cutname);


    free(cutCG->cutElem);
    free(cutCG->cutCoef);
    free( cutCG );
    *_cutCG = NULL;

}

