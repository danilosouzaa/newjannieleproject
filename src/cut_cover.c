
/* GOAL PSP Solver II
 * Resource-Constrained Project Scheduling Problems (RCPSP)
 *
 * Develop as part of the D.Sc. thesis of Araujo, Janniele A. S., with collaboration
 *                                   of Santos, H.G.
 */

#include "cut_cover.h"
#define VERBOSE  1

/*create and initialize the structure of cutC*/
CutC *CutC_create(const Instance *inst)
{

    CutC *ccr;
    ALLOCATE_INI(ccr,CutC);

    ALLOCATE_VECTOR(ccr->cutrhs,VecDbl*, Inst_nResR(inst));
    ALLOCATE_VECTOR(ccr->cutdominated,VecInt*, Inst_nResR(inst));
    ALLOCATE_VECTOR(ccr->cutviolation,VecDbl*, Inst_nResR(inst));
    ALLOCATE_VECTOR(ccr->cutnelem,VecInt*, Inst_nResR(inst));
    ALLOCATE_VECTOR(ccr->cutsense,VecInt*, Inst_nResR(inst));
    ALLOCATE_VECTOR(ccr->cutname,VecStr*, Inst_nResR(inst));
    ALLOCATE_VECTOR(ccr->nAllocations,int, Inst_nResR(inst));


    for(int r = 0 ; r < Inst_nResR(inst); r++) {
        ccr->nAllocations[r] = 0;
        ccr->cutrhs[r] = VDbl_create();
        ccr->cutdominated[r] = VInt_create();
        ccr->cutviolation[r] = VDbl_create();
        ccr->cutnelem[r] = VInt_create();
        ccr->cutsense[r] = VInt_create();
        ccr->cutname[r] = VStr_create(256);
    }

    ccr->inst = inst;

    return ccr;
}

/*Lifted MIP to separate cover cuts,  also verify dominance between generated cuts*/
void CutC_add_cuts_cover_model_parallel( CutC *ccr, LinearProgram *lp, const int *origCols,  LinearProgram *origLP, const Instance *inst, double timeLeft, int continuous, int lifting,  int nround )
{

    assert(lp);
    assert(inst);
    assert(origLP);

    double init = omp_get_wtime();
    int nCut = 0;
    double maxRC = lp_get_max_reduced_cost(origLP);
    const double *rdc = lp_reduced_cost(lp);

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
    int firstTJ[nJobs];
    FILL( firstTJ, 0, nJobs, INT_MAX );
    int lastTJ[nJobs];
    FILL( lastTJ, 0, nJobs, 0 );

    // checking maximum time where there is some allocation
    for ( i=0 ; (i<lp_cols(lp)) ; ++i ) {
        //if (fabs(xf[i])>1e-5) {
        if(maxRC != -1.0)
            if(rdc[i] > maxRC) continue;

        lp_col_name( origLP, origCols[i], name );
        parseName( name, prefix, idx );

        if (prefix[0]=='x') {
            const int j = idx[0];
            const int t = idx[2];

            maxT = MAX( maxT, t );
            firstTJ[j] = MIN(firstTJ[j], t);
            lastTJ[j] = MAX(lastTJ[j], t);

        }
    }

    int nTimes = maxT+1;
    double ****sumX;//[nTimes][Inst_nResR(inst)][Inst_nJobs(inst)][Inst_nMaxModes(inst)];
    double ****sumRC;//[nTimes][Inst_nResR(inst)][Inst_nJobs(inst)][Inst_nMaxModes(inst)];


    //double sumX[nTimes][Inst_nResR(inst)][Inst_nJobs(inst)][Inst_nMaxModes(inst)];
    //double sumRC[nTimes][Inst_nResR(inst)][Inst_nJobs(inst)][Inst_nMaxModes(inst)];

    ALLOCATE_VECTOR( sumX, double ***, nTimes );
    ALLOCATE_VECTOR( sumRC, double ***, nTimes );
    for ( int t=0; (t<nTimes) ; ++t ){
        ALLOCATE_VECTOR( sumX[t], double **, Inst_nResR(inst) );
        ALLOCATE_VECTOR( sumRC[t], double **, Inst_nResR(inst) );
        for ( int r=0 ; (r<Inst_nResR(inst)) ; ++r ){
            ALLOCATE_VECTOR_INI( sumX[t][r], double *, Inst_nJobs(inst) );
            ALLOCATE_VECTOR_INI( sumRC[t][r], double * , Inst_nJobs(inst));
            for ( int j=0 ; (j<Inst_nJobs(inst)) ; ++j ){
                ALLOCATE_VECTOR_INI( sumX[t][r][j], double , Inst_nMaxModes(inst) );
                ALLOCATE_VECTOR_INI( sumRC[t][r][j], double , Inst_nMaxModes(inst));
            }
        }
    }


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
    for ( i=0 ; (i<lp_cols(lp)) ; ++i ) {
        if(maxRC != -1.0)
            if(rdc[i] > maxRC) continue;

        lp_col_name( origLP, origCols[i], name );
        if (tolower(name[0]=='x')) {
            parseName( name, prefix, idx );
            //  if (fabs(xf[i])>1e-5) {
            if (prefix[0]=='x')
                x[idx[0]][idx[1]][idx[2]] = xf[i];
            // }
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
    for ( i=0 ; (i<lp_cols(lp)) ; ++i ) {

        if(maxRC != -1.0) {
            //  printf("\nmaxRC %f , rdc[%d] %f\n", maxRC, i, rdc[i]);
            if(rdc[i] > maxRC) continue;
        }

        lp_col_name( origLP, origCols[i], name );
        if (tolower(name[0])=='x') {
            parseName( name, prefix, idx );
            // if (fabs(xf[i])>1e-5) {
            if (prefix[0]=='x')
                xIdx[idx[0]][idx[1]][idx[2]] = i;
            // }
        }
    }


    /* z indexes new problem*/
    int ****addedTRJM;
    ALLOCATE_VECTOR(addedTRJM, int***,nTimes);

    VecStr ***namesbin;
    ALLOCATE_VECTOR(namesbin, VecStr**, nTimes);
    int ****zIdxm;
    ALLOCATE_VECTOR( zIdxm, int ***, nTimes );
    int ****wIdx;
    ALLOCATE_VECTOR( wIdx, int ***, nTimes );
    for(int t  = 0 ; t < nTimes ; t++) {
        ALLOCATE_VECTOR(namesbin[t], VecStr*,  Inst_nResR(inst));
        ALLOCATE_VECTOR( addedTRJM[t], int **, Inst_nResR(inst) );

        ALLOCATE_VECTOR( zIdxm[t], int **, Inst_nResR(inst) );
        ALLOCATE_VECTOR( wIdx[t], int **, Inst_nResR(inst) );
        for ( int r=0 ; (r<Inst_nResR(inst)) ; ++r ) {
            namesbin[t][r] = VStr_create( STR_SIZE );
            ALLOCATE_VECTOR( addedTRJM[t][r], int *, nJobs );
            ALLOCATE_VECTOR( zIdxm[t][r], int *, nJobs );
            ALLOCATE_VECTOR( wIdx[t][r], int *, nJobs );
            for(int j = 0 ; j < nJobs ; j++) {
                ALLOCATE_VECTOR( addedTRJM[t][r][j], int, nModes );
                ALLOCATE_VECTOR( zIdxm[t][r][j], int, nModes );
                ALLOCATE_VECTOR( wIdx[t][r][j], int, nModes );
                for(int m = 0 ; m <nModes ; m++) {
                    addedTRJM[t][r][j][m] = 0;
                    zIdxm[t][r][j][m] = -1;
                    wIdx[t][r][j][m] = -1;
                    sumX[t][r][j][m] = 0;
                    sumRC[t][r][j][m] = 0;

                }
            }
        }
    }


    VecInt ***idxJ;
    ALLOCATE_VECTOR( idxJ, VecInt**, nTimes );
    int ***yIdx;
    ALLOCATE_VECTOR( yIdx, int **, nTimes );
    for(int t  = 0 ; t < nTimes ; t++) {
        ALLOCATE_VECTOR( idxJ[t], VecInt*, Inst_nResR(inst) );
        ALLOCATE_VECTOR( yIdx[t], int *, Inst_nResR(inst) );
        for ( int r=0 ; (r<Inst_nResR(inst)) ; ++r ) {
            idxJ[t][r] = VInt_create();
            ALLOCATE_VECTOR( yIdx[t][r], int, nJobs );
            for(int j = 0 ; j < nJobs ; j++)
                yIdx[t][r][j] = -1;
        }
    }



    VecInt ***idxResRZ;
    VecInt ***idxResRW;
    ALLOCATE_VECTOR( idxResRZ, VecInt**, Inst_nResR(inst) );
    ALLOCATE_VECTOR_INI( idxResRZ[0], VecInt*, Inst_nResR(inst)*(nTimes) );
    ALLOCATE_VECTOR( idxResRW, VecInt**, Inst_nResR(inst) );
    ALLOCATE_VECTOR_INI( idxResRW[0], VecInt*, Inst_nResR(inst)*(nTimes) );
    for ( int i=1 ; (i<Inst_nResR(inst)) ; ++i ) {
        idxResRZ[i] = idxResRZ[i-1] + nTimes;
        idxResRW[i] = idxResRW[i-1] + nTimes;
    }
    for ( int r=0 ; (r<Inst_nResR(inst)) ; ++r )
        for ( int t=0 ; (t<(maxT+1)) ; ++t ) {
            idxResRZ[r][t] = VInt_create();
            idxResRW[r][t] = VInt_create();
        }

    VecDbl ***coefResRZ;
    VecDbl ***coefResRW;
    ALLOCATE_VECTOR( coefResRZ, VecDbl**, Inst_nResR(inst) );
    ALLOCATE_VECTOR_INI( coefResRZ[0], VecDbl*, Inst_nResR(inst)*(maxT+1) );
    ALLOCATE_VECTOR( coefResRW, VecDbl**, Inst_nResR(inst) );
    ALLOCATE_VECTOR_INI( coefResRW[0], VecDbl*, Inst_nResR(inst)*(maxT+1) );
    for ( int i=1 ; (i<Inst_nResR(inst)) ; ++i ) {
        coefResRZ[i] = coefResRZ[i-1] + nTimes;
        coefResRW[i] = coefResRW[i-1] + nTimes;
    }

    for ( int r=0 ; (r<Inst_nResR(inst)) ; ++r )
        for ( int t=0 ; (t<nTimes) ; ++t ) {
            coefResRZ[r][t] = VDbl_create();
            coefResRW[r][t] = VDbl_create();
        }



    VecInt ****vecTSameJZ;
    VecInt ****vecTSameJW;
    // VecDbl *****vecTSameJZcoef;
    double ****vecTSameJWcoef;
    // VecInt *****vecJMWR;
    ALLOCATE_VECTOR( vecTSameJZ, VecInt***, nTimes);
    ALLOCATE_VECTOR( vecTSameJW, VecInt***, nTimes);
    //ALLOCATE_VECTOR( vecTSameJZcoef, VecDbl****, nTimes);
    ALLOCATE_VECTOR( vecTSameJWcoef, double***, nTimes);
    //ALLOCATE_VECTOR( vecJMWR, VecInt****, nTimes);
    for ( int i=0 ; (i<nTimes) ; ++i ) {
        ALLOCATE_VECTOR( vecTSameJZ[i], VecInt**, nJobs);
        ALLOCATE_VECTOR( vecTSameJW[i], VecInt**, nJobs);
        // ALLOCATE_VECTOR( vecTSameJZcoef[i], VecDbl***, nJobs);
        ALLOCATE_VECTOR( vecTSameJWcoef[i], double**, nJobs);
        //ALLOCATE_VECTOR( vecJMWR[i], VecInt***, nJobs);
        for ( int j=0 ; (j<nJobs) ; ++j ) {
            ALLOCATE_VECTOR( vecTSameJZ[i][j], VecInt*, Inst_nResR(inst));
            ALLOCATE_VECTOR( vecTSameJW[i][j], VecInt*, Inst_nResR(inst));
            //ALLOCATE_VECTOR( vecTSameJZcoef[i][j], VecDbl**, Inst_nResR(inst));
            ALLOCATE_VECTOR( vecTSameJWcoef[i][j], double*, Inst_nResR(inst));
            //ALLOCATE_VECTOR( vecJMWR[i][j], VecInt**, Inst_nResR(inst));
            for ( int r=0 ; (r<Inst_nResR(inst)) ; ++r ) {
                vecTSameJZ[i][j][r] = VInt_create();
                vecTSameJW[i][j][r] = VInt_create();
                //ALLOCATE_VECTOR( vecTSameJZcoef[i][j][r], VecDbl*, Inst_nMaxModes(inst));
                //ALLOCATE_VECTOR( vecTSameJWcoef[i][j][r], VecDbl*, Inst_nMaxModes(inst));
                ALLOCATE_VECTOR( vecTSameJWcoef[i][j][r], double, Inst_nMaxModes(inst));
                // ALLOCATE_VECTOR( vecJMWR[i][j][r], VecInt*, Inst_nMaxModes(inst));
                for ( int m=0 ; (m<Inst_nMaxModes(inst)) ; ++m ) {
                    //  vecTSameJZcoef[i][j][r][m] = VDbl_create();
                    //  vecTSameJWcoef[i][j][r][m] = VDbl_create();
                    vecTSameJWcoef[i][j][r][m] = 0;
                }

            }
        }
    }


    int **nColB;
    ALLOCATE_VECTOR_INI(nColB,int*,nTimes);
    double ***objbin;
    ALLOCATE_VECTOR(objbin,double**,nTimes);
    for(int o = 0 ; o < nTimes; o++) {
        ALLOCATE_VECTOR_INI(objbin[o],double*,Inst_nResR(inst));
        ALLOCATE_VECTOR_INI(nColB[o],int,Inst_nResR(inst));
        for ( int r=0 ; (r<Inst_nResR(inst)) ; ++r )
            ALLOCATE_VECTOR_INI(objbin[o][r],double,nJobs*Inst_nMaxModes(inst)*nTimes);
    }


    for ( int t=0 ; t< nTimes; ++t ) {
        for ( int j= 0 ; j< nJobs ; ++j ) {
            const Job *job = Inst_job( inst, j );
            const int nModes = Job_nModes( job );
            for ( int m=0 ; (m<nModes) ; ++m ) {
                const Mode *mode = Job_mode( job, m );
                if(Mode_duration(mode)==0) continue;
                // if(xIdx[j][m][t] == -1 || x[j][m][t] >= 0.9995) continue;
                //if(xIdx[j][m][t] == -1 || fabs(x[j][m][t]) == 1.0) continue;

                if (!Mode_isFeasible(inst,mode))
                    continue;


                /* adding resource constraints */
                for ( int ir=0 ; (ir<Mode_nResR(mode)) ; ++ir ) {
                    int idzR = Mode_idxResR( mode, ir );
                    double use = (double) Mode_useResR( mode, ir );

                    for ( int tj=t-Mode_duration(mode)+1 ; (tj<=t && tj<nTimes) ; ++tj )
                        //for ( int tj=t ; (tj<t+Mode_duration(mode) && tj<nTimes) ; ++tj )
                    {
                        if(tj<0) {
                            tj = 0;
                            // continue;
                        }

                        if(xIdx[j][m][tj] == -1 || fabs(x[j][m][tj]) == 1.0) continue;
//

                        // if( wIdx[tj][idzR][j][m] == -1 && zIdxm[tj][idzR][j][m] == -1 ) continue;
                        sumX[t][idzR][j][m] += x[j][m][tj];
                        sumRC[t][idzR][j][m] += rdc[xIdx[j][m][tj]];




                        if(addedTRJM[t][idzR][j][m] == 0) {

                            char vnameZ[STR_SIZE];
                            sprintf( vnameZ, "z(%d,%d)", j, m);
                            zIdxm[t][idzR][j][m] =  nColB[t][idzR] ;
                            nColB[t][idzR]++;
                            VStr_pushBack( namesbin[t][idzR], vnameZ );
                            if(yIdx[t][idzR][j]==-1) {
                                char vnameY[STR_SIZE];
                                sprintf( vnameY, "y(%d)", j);
                                objbin[t][idzR][nColB[t][idzR]] = 0.0;
                                yIdx[t][idzR][j] = nColB[t][idzR];
                                VStr_pushBack( namesbin[t][idzR], vnameY );
                                VInt_pushBack( idxJ[t][idzR], j );
                                nColB[t][idzR]++;
                            }

                            char vnameW[STR_SIZE];
                            sprintf( vnameW, "w(%d,%d)", j, m);
                            wIdx[t][idzR][j][m] =  nColB[t][idzR];
                            nColB[t][idzR]++;
                            VStr_pushBack( namesbin[t][idzR], vnameW );

                            //  printf("T%d,R%d,J%d,M%d: (maxRC %f - sumRC %f)* sumX %f: objint %f\n", tj, idzR,j,m,maxRC,sumRC[tj][idzR][j][m],sumX[tj][idzR][j][m], objbin[tj][idzR][zIdx[t][idzR][j][m]]);
                            // getchar();
                            // double obj = ((maxRC-sumRC[tj][idzR][j][m])*sumX[tj][idzR][j][m]);
                            // if(obj==0) obj= 0.1;
                            objbin[t][idzR][zIdxm[t][idzR][j][m]] = - 0.1;


                            //objbin[tj][idzR][zIdxm[tj][idzR][j][m]] = 0.1;
                            objbin[t][idzR][wIdx[t][idzR][j][m]] = 0.0;


                            VInt_pushBack( vecTSameJZ[t][j][idzR], zIdxm[t][idzR][j][m]);
                            VInt_pushBack( vecTSameJW[t][j][idzR], wIdx[t][idzR][j][m] );

                            //VDbl_pushBack( vecTSameJZcoef[tj][j][idzR][m], use );
                            //VDbl_pushBack( vecTSameJWcoef[tj][j][idzR][m], use );
                            vecTSameJWcoef[t][j][idzR][m] = use;
                            VInt_pushBack( idxResRZ[idzR][t], zIdxm[t][idzR][j][m] );
                            VDbl_pushBack( coefResRZ[idzR][t], use );

                            VInt_pushBack( idxResRW[idzR][t], wIdx[t][idzR][j][m] );
                            VDbl_pushBack( coefResRW[idzR][t], use );

                            addedTRJM[t][idzR][j][m] = 1;
                        }
                    }


                }
            }
        }
    }


    LinearProgram *mipCutm;
    ALLOCATE_VECTOR(ccr->cutElem,VecInt**,Inst_nResR(inst));
    ALLOCATE_VECTOR(ccr->cutCoef,VecDbl**,Inst_nResR(inst));
    int *nCutR;
    ALLOCATE_VECTOR_INI(nCutR,int,Inst_nResR(inst));


    for(int r = 0 ; r < Inst_nResR(inst); r++ ) {
        ALLOCATE_VECTOR(ccr->cutElem[r],VecInt*,nTimes);
        ALLOCATE_VECTOR(ccr->cutCoef[r],VecDbl*,nTimes);
        for(int nc = 0 ; nc < nTimes; nc++ ) {
            ccr->cutElem[r][nc] = VInt_create();
            ccr->cutCoef[r][nc] = VDbl_create();
            ccr->nAllocations[r]++;
        }
    }



    for ( int t=0 ; t<nTimes; ++t ) {
        for ( int r=0 ; (r<Inst_nResR(inst)) ; ++r ) {

            if(nColB[t][r]==0) continue;
            mipCutm = lp_create();
            double _time = (double) timeLeft - (omp_get_wtime()-init);
            lp_set_print_messages(mipCutm,0);
            lp_set_max_seconds(mipCutm,_time);
            lp_add_bin_cols( mipCutm, nColB[t][r], objbin[t][r], VStr_ptr(namesbin[t][r]) );

            VecStr *namesint = VStr_create(STR_SIZE);
            char vnamee[256];
            sprintf( vnamee, "e");
            VStr_pushBack( namesint, vnamee );
            char vnamev[256];
            sprintf( vnamev, "v");
            VStr_pushBack( namesint, vnamev );

            double *objint;
            ALLOCATE_VECTOR_INI(objint,double,2);

            char *integer;
            ALLOCATE_VECTOR(integer,char,2);
            integer[0] = True;
            integer[1] = False;

            double *lb;
            ALLOCATE_VECTOR_INI(lb,double,2);
            double *up;
            ALLOCATE_VECTOR_INI(up,double,2);
            FILL(lb,0,2,0.0);
            FILL(up,0,2,INT_MAX);
            objint[0] = 0.0;
            objint[1] = -100000;

            lp_add_cols(mipCutm,2,objint,lb,up,integer,VStr_ptr( (VecStr *)namesint));


            free(integer);
            free(lb);
            free(up);
            free(objint);
            VStr_free(&namesint);
            int nColsLP = lp_cols(lp);

            const VecInt *vidzm = idxResRZ[r][t];
            const VecInt *vidwm = idxResRW[r][t];

            const int nwt = VInt_size(vidwm);
            const int nzt = VInt_size(vidzm);

            const VecInt *vj = idxJ[t][r];
            const int nj = VInt_size( vj );
            int *idj =  VInt_getPtr( (VecInt *) vj );


            if(nj && nwt && nzt ) {

                int idwr[nj+nj*nzt+2];
                double coefwr[nj+nj*nzt+2];
                int ijm = 0;
                int idzyv[nj+nj*nzt+2];
                double coefzyv[nj+nj*nzt+2];
                int ijz = 0;

                for(int iy = 0 ; iy <nj ; iy++) {
                    if( yIdx[t][r][idj[iy]]==-1)continue;
                    // int idzyv2[nj+nj*nzt+2];
                    //double coefzyv2[nj+nj*nzt+2];
                    //int ijz2 = 0;

                    int iii = 0;
                    const Job* jobw = Inst_job(inst,idj[iy]);
                    int idSameJW[Job_nModes(jobw)+1];
                    double coefSameJW[Job_nModes(jobw)+1];
                    FILL( coefSameJW, 0, Job_nModes(jobw)+1, -1.0 );

                    int contjm=0;
                    for(int jmw = 0 ; jmw < Job_nModes(jobw) ; jmw++) {

                        int idxmoderes[Job_nModes(jobw)+1];
                        double coefmoderes[Job_nModes(jobw)+1];
                        int imr = 0;
                        if( wIdx[t][r][idj[iy]][jmw] == -1 || zIdxm[t][r][idj[iy]][jmw] ==-1 ) continue;
                        idSameJW[iii] = wIdx[t][r][idj[iy]][jmw];
                        iii++;
                        idwr[ijm] = wIdx[t][r][idj[iy]][jmw];
                        coefwr[ijm] = vecTSameJWcoef[t][idj[iy]][r][jmw];//vecTSameJWcoef[t][idj[iy]][r][jmw];
                        ijm++;


                        if(sumX[t][r][idj[iy]][jmw]>EPS) {
                            idzyv[ijz] = zIdxm[t][r][idj[iy]][jmw];
                            coefzyv[ijz] = sumX[t][r][idj[iy]][jmw];
                            ijz++;
                        }

                        // idzyv2[ijz2] = zIdxm[t][r][idj[iy]][jmw];
                        // coefzyv2[ijz2] = 1.0;
                        // ijz2++;

                        // double sumXobj = ((maxRC-sumRC[t][r][idj[iy]][jmw])*sumX[t][r][idj[iy]][jmw]);
                        //if(sumXobj==0) sumXobj= 0.1;
                        //coefzyv[ijz] = sumXobj;

                        //idwz[0] = wIdx[t][r][idj[iy]][jmw];
                        //idwz[1] = zIdxm[t][r][idj[iy]][jmw];
                        //coefwz[0] = 1.0;
                        //coefwz[1] = -1.0;
                        // printf("CoefSumX %f\n", coefzyv[ijz]); fflush(stdout);
                        contjm++;
                        for(int jmw2 = 0 ; jmw2 < Job_nModes(jobw); jmw2++) {
                            if( wIdx[t][r][idj[iy]][jmw2] == -1 ) continue;
                            //if(jmw2==jmw) continue;
                            if(vecTSameJWcoef[t][idj[iy]][r][jmw2]<=vecTSameJWcoef[t][idj[iy]][r][jmw]) {
                                // char var[256];
                                //sprintf(var,"w(%d,%d)",idj[iy],jmw2);
                                idxmoderes[imr] = wIdx[t][r][idj[iy]][jmw2];// lp_col_index(mipCut,var);
                                coefmoderes[imr] = -1.0;
                                imr++;
                            }
                        }
                        idxmoderes[imr] = zIdxm[t][r][idj[iy]][jmw];
                        coefmoderes[imr] = 1.0;
                        imr++;


                        char namewz[STR_SIZE];
                        sprintf( namewz, "cons(3).(%d,%d)#(%d,%d)",t,r,idj[iy],jmw);
                        // printf("cons(3).(%d,%d)#(%d,%d)",t,r,idj[iy],jmw);
                        lp_add_row( mipCutm, imr, idxmoderes, coefmoderes, namewz, 'L', 0.0);

                        //char namewz2[STR_SIZE];
                        // sprintf( namewz2, "cons(3.1).(%d,%d)#(%d,%d)",t,r,idj[iy],jmw);
                        // printf("cons(3).(%d,%d)#(%d,%d)",t,r,idj[iy],jmw);
                        // lp_add_row( mipCutm, 2, idwz, coefwz, namewz2, 'E', 0.0);

                    }



                    int coly = yIdx[t][r][idj[iy]];
                    idSameJW[iii] = coly;
                    coefSameJW[iii] = 1.0;
                    char namewy[STR_SIZE];
                    iii++;
                    sprintf( namewy, "cons(1).(%d,%d)(%d)",t,r,idj[iy]);
                    lp_add_row( mipCutm, iii, idSameJW, coefSameJW, namewy, 'E',0.0);

                    idzyv[ijz] = coly;
                    coefzyv[ijz] =-1.0;//contjw;
                    ijz++;
                    // idzyv2[ijz2] = coly;
                    //coefzyv2[ijz2] = -contjm;
                    //ijz2++;


                    // char name3[STR_SIZE];
                    // sprintf( name3, "cons(7).(%d,%d)",t,r);
                    //lp_add_row( mipCutm, ijz2, idzyv2,coefzyv2, name3, 'L', 0 );

                }

                idwr[ijm] = lp_col_index(mipCutm,"e");
                coefwr[ijm] = -1.0;
                ijm++;

                char name[STR_SIZE];
                sprintf( name, "cons(2).(%d,%d)",t,r);
                lp_add_row( mipCutm, ijm, idwr, coefwr, name, 'E', Inst_capResR(inst, r) );


                //  char name1[STR_SIZE];
                // sprintf( name1, "cons(2.2).(%d,%d)",t,r);
                // lp_add_row( mipCutm, ijm-1, idwr, coefwr, name1, 'G', Inst_capResR(inst, r)+1 );

                char nameel[STR_SIZE];
                sprintf( nameel, "cons(4).(%d,%d)",t,r);
                int idel[1];
                double coefel[1];
                idel[0] = lp_col_index(mipCutm,"e");
                coefel[0] = 1.0;
                lp_add_row( mipCutm, 1, idel, coefel, nameel, 'G', 1 );

                idzyv[ijz] = lp_col_index(mipCutm,"v");
                fflush(stdout);
                coefzyv[ijz] = -1.0;
                ijz++;

                char namezyv[STR_SIZE];
                sprintf( namezyv, "cons(5).(%d,%d)",t,r);
                // printf("ijz %d, size(idzyv) %d size(coefzyv) %d\n", ijz, sizeof(idzyv)/sizeof(idzyv[0]), sizeof(coefzyv)/sizeof(coefzyv[0])); fflush(stdout);
                lp_add_row( mipCutm, ijz, idzyv, coefzyv, namezyv, 'E', -1 );



                char namev[STR_SIZE];
                sprintf( namev, "cons(6).(%d,%d)",t,r);
                int idev[1];
                double coefv[1];
                idev[0] = lp_col_index(mipCutm,"v");
                coefv[0] = 1;
                lp_add_row( mipCutm, 1, idev, coefv,namev, 'G', 0.005 );

                //lp_set_concurrentMIP(mipCutm,1);
                //lp_set_method(mipCutm,4);
                //lp_set_seed(mipCutm,100000);
                int st = lp_optimize(mipCutm);


                /* if(st==0)
                 {
                     printf("creating mipCutm RR: end\n");
                     fflush(stdout);
                     lp_write_sol(mipCutm,"mipCutm.sol");
                     zf2 = lp_x(mipCutm);
                     getchar();

                 lp_write_lp(mipCutm,"mipCutm.lp");
                 lp_conflicts_model(mipCutm);
                 }
                 */



                if(st==0) {
                    double valueZ =1;
                    char vvalue[256];
                    sprintf(vvalue,"v");
                    int vidxvalue = lp_col_index(mipCutm,vvalue);
                    const double *zf2 = lp_x(mipCutm);
                    valueZ = zf2[vidxvalue];
                    //valueZ = lp_obj_value(mipCutm);
                    //  printf("lp_obj_value(mipCutm) %f, zf2[vidxvalue] %f\n",lp_obj_value(mipCutm), zf2[vidxvalue] );
                    if(valueZ >= 0.0005) {


                        // printf("find real cuts mipCut RR to resource %d on time %d : start\n", r,t);fflush(stdout);
                        //    lp_write_lp(mipCut,"lp_cut.lp");
                        //    printf(" writing mipCut st %d end \n", st);
                        //    printf("valueZ %g t %d : ", valueZ,t);
                        //   getchar();

                        IntDblPair *cutIdxCo;
                        ALLOCATE_VECTOR(cutIdxCo, IntDblPair, nColsLP);


                        // filling x from fractional solution
                        int c=0, nCutReal =0, cs =0; //minDuration =INT_MAX, cMinDuration=0 ;
                        IntTriple *varCut;
                        ALLOCATE_VECTOR_INI(varCut,IntTriple, nJobs*Inst_nMaxModes(inst)*nTimes)
                        IntTriple *varCutOrig;
                        ALLOCATE_VECTOR_INI(varCutOrig,IntTriple, nJobs*Inst_nMaxModes(inst)*nTimes)

                        // int cZidx[nJobs*nModes*nTimes];

                        int nc = lp_cols(mipCutm);


                        double slack = 0;

                       // printf("\n new t %d: ", t);
                        int ny = 0;
                        /* VecInt **jobsincluded;
                         ALLOCATE_VECTOR_INI(jobsincluded,int*,nJobs);
                         FILL(jobsincluded,0,nJobs,0);
                         for(int iikj = 0 ; iikj < nJobs ; iikj++)
                             jobsincluded[iikj] = VInt_create();
                         */
                        for ( i=0 ; (i<nc) ; ++i ) {
                            char nameC[STR_SIZE];

                            lp_col_name( mipCutm, i, nameC );

                            if (tolower(nameC[0]=='y')) {

                                if (zf2[i]>EPS) {
                                    parseName( nameC, prefix, auxidZ );
                                    ny++;
                                }
                            }
                            if (tolower(nameC[0]=='z')) {
                                parseName( nameC, prefix, auxidZ );

                                if (zf2[i]>EPS) {
                                    // printf(" NEW zf[%d] %f > EPS %f ", i, zf2[i], EPS); fflush(stdout);
                                    if (prefix[0]=='z') {


                                        int j = auxidZ[0];

                                        int m = auxidZ[1];
                                        // VInt_pushBack(jobsincluded[j],m);
                                        const Job *jobOnCut = Inst_job(inst,j);
                                        const Mode* modeOnCut = Job_mode(jobOnCut,m);
                                        if(Mode_duration(modeOnCut)==0) continue;
                                        for(int tt = t-Mode_duration(modeOnCut)+1 ; tt<=t ; tt++) {

                                            assert(j<nJobs);
                                            assert(m<nModes);

                                            assert(c<nJobs*nModes*nTimes);
                                            if(xIdx[j][m][tt] == -1) continue;

                                            // cZidx[c] = xIdx[j][m][tt];
                                            varCut[c].j = j;
                                            varCut[c].m = m;
                                            varCut[c].t = tt;
                                            varCut[c].value = zf2[i];

                                            varCutOrig[cs].j = j;
                                            varCutOrig[cs].m = m;
                                            varCutOrig[cs].t = tt;
                                            varCutOrig[cs].value = zf2[i];

                                            cutIdxCo[c].a =  xIdx[j][m][tt];//cZidx[c];
                                            cutIdxCo[c].b = 1.0;
                                            slack += xf[cutIdxCo[c].a]*cutIdxCo[c].b;


                                            nCutReal++;
                                            c++;
                                          //  printf("O (%d,%d,%d), ", j,m,tt);
                                            fflush(stdout);
                                            /*  if( tt+Mode_duration(modeOnCut) <= minDuration)
                                              {
                                                  minDuration = tt+Mode_duration(modeOnCut);
                                                  cMinDuration = cs;
                                              }*/
                                            cs++;
                                        }
                                    }
                                }
                            }
                        }


                        // printf("slack-(ny-1) %f\n", slack-(ny-1));
                        //if((slack-(ny-1)) >= 0.0002)
                        // {

                        char nameCut[STR_SIZE];
                        //printf( "cutRR(%d,%d)#%d\n", r, t, nCut );
                        //   if(r==1 && t==5 && nCut == 1274) getchar();
                        //   if(r==1 && t==6 && nCut == 1588) getchar();
                        sprintf( nameCut, "cutRR(%d,%d)#%d#%d", r, t, nround, nCut );
                        VStr_pushBack(ccr->cutname[r],nameCut);
                        // double coe[c];
//                        getchar();
                        CutP_quick_sort(cutIdxCo, c);
                        for(int sel  = 0 ; sel <c; sel++) {
                            int elem =cutIdxCo[sel].a;
                            //  cZidx[sel] = elem;
                            // coe[sel] = cutIdxCo[sel].b;

                            VInt_pushBack(ccr->cutElem[r][nCutR[r]],elem);
                            VDbl_pushBack(ccr->cutCoef[r][nCutR[r]],cutIdxCo[sel].b);
                        }
                        VDbl_pushBack(ccr->cutrhs[r],ny-1);
                        VInt_pushBack(ccr->cutsense[r],0);
                        VInt_pushBack(ccr->cutdominated[r],0);
                        VInt_pushBack(ccr->cutnelem[r],c);
                        VDbl_pushBack(ccr->cutviolation[r],slack-(ny-1));

                        //double newrhs = CutP_model_lift( inst, VStr_get(ccr->cutname[r],nCutR[r]), VInt_get(ccr->cutnelem[r],nCutR[r]),VInt_getPtr(ccr->cutElem[r][nCutR[r]]), VDbl_getPtr(ccr->cutCoef[r][nCutR[r]]), VDbl_get(ccr->cutrhs[r],nCutR[r]), lp, _time);
                        //ir(newrhs!=-1)
                          //  VDbl_set(ccr->cutrhs,r,newrhs);
                        //   int *cutIdx = VInt_getPtr(ccr->cutElem[nCut]);
                        // double *cutCoef = VDbl_getPtr(ccr->cutCoef[nCut]);
                        // CutP_quick_sort_vec(cutIdx,cutCoef, c);
                        nCutR[r]++;

                        nCut++;
                        //}
                        free(cutIdxCo);

                        free(varCut);
                        free(varCutOrig);
                    }
                }
            }
            lp_free(&mipCutm);

        }
    }
    free(nCutR);
    for(int r = 0 ; r < Inst_nResR(inst); r++ ) {
        // printf("VInt_size(ccr->cutrhs[r]) %d: %d\n",VInt_size(ccr->cutrhs[r]) ,  VInt_size(ccr->cutrhs[r]) *VInt_size(ccr->cutrhs[r]) );
        for(int ncdA = 0; ncdA < VDbl_size(ccr->cutrhs[r]) ; ncdA++) {
            for(int ncdB = ncdA+1; ncdB < VDbl_size(ccr->cutrhs[r]) ; ncdB++) {
                if(VInt_get(ccr->cutdominated[r],ncdA) == 1) break;
                if(VInt_get(ccr->cutdominated[r],ncdB) == 1) continue;
                CutP_dominanceBetweenTwo( lp, ncdA, VInt_getPtr(ccr->cutElem[r][ncdA]), VDbl_getPtr(ccr->cutCoef[r][ncdA]), VDbl_get(ccr->cutrhs[r],ncdA), VInt_get(ccr->cutsense[r],ncdA), VInt_get(ccr->cutnelem[r],ncdA), ncdB, VInt_getPtr(ccr->cutElem[r][ncdB]), VDbl_getPtr(ccr->cutCoef[r][ncdB]), VDbl_get(ccr->cutrhs[r],ncdB), VInt_get(ccr->cutsense[r],ncdB), VInt_get(ccr->cutnelem[r],ncdB),  ccr->cutdominated[r]);
                // getchar();
            }
        }
    }
    for ( int i=0 ; (i<nTimes) ; ++i ) {
        for ( int j=0 ; (j<nJobs) ; ++j ) {
            for ( int r=0 ; (r<Inst_nResR(inst)) ; ++r ) {
                VInt_free( &vecTSameJZ[i][j][r]);
                VInt_free( &vecTSameJW[i][j][r]);
                //for ( int m=0 ; (m<Inst_nMaxModes(inst)) ; ++m )
                //{
                //  free(vecTSameJWcoef[i][j][r][m]);
                //}
                free( vecTSameJWcoef[i][j][r]);
            }
            free( vecTSameJZ[i][j]);
            free( vecTSameJW[i][j]);
            free( vecTSameJWcoef[i][j]);
        }
        free( vecTSameJZ[i]);
        free( vecTSameJW[i]);
        free( vecTSameJWcoef[i]);
    }
    free( vecTSameJZ);
    free( vecTSameJW);
    free( vecTSameJWcoef);

    for ( int i=0 ; (i<Inst_nResR( inst )*(nTimes)) ; ++i ) {
        VInt_free( &idxResRZ[0][i] );
        VDbl_free( &coefResRZ[0][i] );
        VInt_free( &idxResRW[0][i] );
        VDbl_free( &coefResRW[0][i] );
    }
    free(idxResRZ[0]);
    free(idxResRW[0]);
    free(coefResRZ[0]);
    free(coefResRW[0]);
    free( idxResRZ );
    free( coefResRZ );
    free( idxResRW );
    free( coefResRW );
    for(int t  = 0 ; t < nTimes ; t++) {
        free(nColB[t]);
        for(int r  = 0 ; r < Inst_nResR(inst) ; r++) {
            free(yIdx[t][r]);

            VStr_free(&namesbin[t][r]);
            free(objbin[t][r]);
            VInt_free(&idxJ[t][r]);
            for(int j = 0 ; j < nJobs ; j++) {
                free(zIdxm[t][r][j]);
                free(wIdx[t][r][j]);
                free(addedTRJM[t][r][j]);
                free(sumX[t][r][j]);
                free(sumRC[t][r][j]);
            }
            free(sumRC[t][r]);
            free(sumX[t][r]);
            free(idxJ[t][r]);
            free(addedTRJM[t][r]);
            free(zIdxm[t][r]);
            free(wIdx[t][r]);
        }
        free(objbin[t]);
        free(namesbin[t]);
        free(addedTRJM[t]);
        free(zIdxm[t]);
        free(wIdx[t]);
        free(yIdx[t]);
        free(idxJ[t]);
        free(sumX[t]);
        free(sumRC[t]);
    }
    free(sumRC);
    free(sumX);
    free(nColB);
    free(objbin);
    free(namesbin);
    free(addedTRJM);
    free(zIdxm);
    free(wIdx);
    free(xIdx[0][0]);
    free(xIdx[0]);
    free(xIdx);
    free(yIdx);
    free(idxJ);
    free(x[0][0]);
    free(x[0]);
    free(x);
#undef MAX_IDX
}

/*Cover procedure based on Zhu 2006, also verify dominance between generated cuts*/
void CutC_add_cuts_cover_parallel( CutC *ccr, LinearProgram *lp, const int *origCols,  LinearProgram *origLP, const Instance *inst, double timeLeft, int continuous, int lifting,  int nround )
{

    assert(lp);
    assert(inst);
    assert(origLP);

//    double init = omp_get_wtime();
    int nCut = 0;
    double maxRC = lp_get_max_reduced_cost(origLP);
    const double *rdc = lp_reduced_cost(lp);

    const double *xf = lp_x(lp);
    char name[256];
    char name2[256];
    int i;

    /* maximum number of indexes in a var name */
#define MAX_IDX 16
    char prefix[STR_SIZE];
    int idx[MAX_IDX];
    int auxidZ[MAX_IDX];
    int maxT = 0;
    sprintf(prefix,"");
    const int nJobs = Inst_nJobs(inst);
    const int nModes = Inst_nMaxModes(inst);

    // first and last time job has some allocation
    int firstTJ[nJobs];
    FILL( firstTJ, 0, nJobs, INT_MAX );
    int lastTJ[nJobs];
    FILL( lastTJ, 0, nJobs, 0 );

    // checking maximum time where there is some allocation
    for ( i=0 ; (i<lp_cols(lp)) ; ++i ) {
        //if (fabs(xf[i])>1e-5) {
        if(maxRC != -1.0)
            if(rdc[i] > maxRC) continue;

        lp_col_name( origLP, origCols[i], name );
        parseName( name, prefix, idx );

        if (prefix[0]=='x') {
            const int j = idx[0];
            const int t = idx[2];

            maxT = MAX( maxT, t );
            firstTJ[j] = MIN(firstTJ[j], t);
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
    for ( i=0 ; (i<lp_cols(lp)) ; ++i ) {
        if(maxRC != -1.0)
            if(rdc[i] > maxRC) continue;

        lp_col_name( origLP, origCols[i], name );
        if (tolower(name[0]=='x')) {
            parseName( name, prefix, idx );
            //  if (fabs(xf[i])>1e-5) {
            if (prefix[0]=='x')
                x[idx[0]][idx[1]][idx[2]] = xf[i];
            // }
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
    for ( i=0 ; (i<lp_cols(lp)) ; ++i ) {

        if(maxRC != -1.0) {
            //  printf("\nmaxRC %f , rdc[%d] %f\n", maxRC, i, rdc[i]);
            if(rdc[i] > maxRC) continue;
        }

        lp_col_name( origLP, origCols[i], name );
        if (tolower(name[0])=='x') {
            parseName( name, prefix, idx );
            // if (fabs(xf[i])>1e-5) {
            if (prefix[0]=='x')
                xIdx[idx[0]][idx[1]][idx[2]] = i;
            // }
        }
    }

    /* z indexes new problem*/
    int *****addedTRJM;
    ALLOCATE_VECTOR(addedTRJM, int****,nTimes);

    VecStr ***names;
    ALLOCATE_VECTOR(names, VecStr**, nTimes);
    int *****zIdx;
    ALLOCATE_VECTOR( zIdx, int ****, nTimes );
    for(int t  = 0 ; t < nTimes ; t++) {
        ALLOCATE_VECTOR(names[t], VecStr*,  Inst_nResR(inst));
        ALLOCATE_VECTOR( addedTRJM[t], int ***, Inst_nResR(inst) );

        ALLOCATE_VECTOR( zIdx[t], int ***, Inst_nResR(inst) );
        for ( int r=0 ; (r<Inst_nResR(inst)) ; ++r ) {
            names[t][r] = VStr_create( STR_SIZE );
            ALLOCATE_VECTOR( addedTRJM[t][r], int **, nJobs );
            ALLOCATE_VECTOR( zIdx[t][r], int **, nJobs );
            for(int j = 0 ; j < nJobs ; j++) {
                ALLOCATE_VECTOR( addedTRJM[t][r][j], int*, nModes );
                ALLOCATE_VECTOR( zIdx[t][r][j], int*, nModes );
                 for(int m = 0 ; m <nModes ; m++) {
                    //addedTRJM[t][r][j][m] = 0;
                    //zIdx[t][r][j][m] = -1;
                    ALLOCATE_VECTOR( addedTRJM[t][r][j][m], int, nTimes );
                    ALLOCATE_VECTOR( zIdx[t][r][j][m], int, nTimes );
                    for(int tt = 0 ; tt <nTimes ; tt++){
                        addedTRJM[t][r][j][m][tt] = 0;
                        zIdx[t][r][j][m][tt] = -1;
                    }
                }
            }
        }
    }

    VecInt ***idxResR;
    VecDbl ***coefResR;
    ALLOCATE_VECTOR( idxResR, VecInt**, Inst_nResR(inst) );
    ALLOCATE_VECTOR_INI( idxResR[0], VecInt*, Inst_nResR(inst)*(nTimes) );
    ALLOCATE_VECTOR( coefResR, VecDbl**, Inst_nResR(inst) );
    ALLOCATE_VECTOR_INI( coefResR[0], VecDbl*, Inst_nResR(inst)*(maxT+1) );

    for ( int i=1 ; (i<Inst_nResR(inst)) ; ++i ){
        idxResR[i] = idxResR[i-1] + nTimes;
        coefResR[i] = coefResR[i-1] + nTimes;
    }
    for ( int r=0 ; (r<Inst_nResR(inst)) ; ++r ){
        for ( int t=0 ; (t<nTimes) ; ++t ){
            idxResR[r][t] = VInt_create();
            coefResR[r][t] = VDbl_create();
        }
    }


    int **nColC;
    ALLOCATE_VECTOR_INI(nColC,int*,nTimes);
    double ***obj;
    ALLOCATE_VECTOR(obj,double**,nTimes);
    int ***zx;
    ALLOCATE_VECTOR(zx,int**,nTimes);

    for(int o = 0 ; o < nTimes; o++) {
        ALLOCATE_VECTOR_INI(obj[o],double*,Inst_nResR(inst));
        ALLOCATE_VECTOR_INI(nColC[o],int,Inst_nResR(inst));
        ALLOCATE_VECTOR_INI(zx[o],int*,Inst_nResR(inst));
       for ( int r=0 ; (r<Inst_nResR(inst)) ; ++r ) {
            ALLOCATE_VECTOR_INI(obj[o][r],double,nJobs*Inst_nMaxModes(inst)*nTimes);
            ALLOCATE_VECTOR_INI(zx[o][r],int,nJobs*Inst_nMaxModes(inst)*nTimes);
        }
    }

  //  printf("find cuts RR: start\n");fflush(stdout);


    for ( int t=0 ; t< nTimes; ++t ) {
        for ( int j= 0 ; j< nJobs ; ++j ) {
            const Job *job = Inst_job( inst, j );
            const int nModes = Job_nModes( job );
            for ( int m=0 ; (m<nModes) ; ++m ) {
                const Mode *mode = Job_mode( job, m );
                if(Mode_duration(mode)==0) continue;
                if (!Mode_isFeasible(inst,mode))
                    continue;
                if(xIdx[j][m][t] == -1 || fabs(x[j][m][t]) == 1.0) continue;
                /* adding resource constraints */
                for ( int ir=0 ; (ir<Mode_nResR(mode)) ; ++ir ) {
                    int idzR = Mode_idxResR( mode, ir );
                    double use = (double) Mode_useResR( mode, ir );
                    for ( int tj=t ; (tj<t+Mode_duration(mode) && tj<nTimes) ; ++tj ) {
                      //      if(j==3 && m==2 && (t == 5 || t==1) && idzR==1){
                       // printf("Trad(t%d,r%d): xIdx[%d][%d][%d] %d fabs(x[j][m][tj])  %f = 1, x[j][m][tj] %f \n ", tj, idzR, j, m, t, xIdx[j][m][t],  (double) fabs(x[j][m][t]),  (double)x[j][m][t] );
                      //  getchar();

                         //   }
                            char vnamein[STR_SIZE];
                            sprintf( vnamein, "z(%d,%d,%d)", j, m, t );
                            obj[tj][idzR][nColC[tj][idzR]] = fabs(1.0-x[j][m][t])+0.0001;
                          //  zIdx[tj][idzR][j][m][t] =  nColC[tj][idzR];
                            zx[tj][idzR][nColC[tj][idzR]] = nColC[tj][idzR] ;
                           // coefzx[tj][idzR][nColC[tj][idzR]] = fabs(x[j][m][t]) ;

                            VStr_pushBack( names[tj][idzR], vnamein );
                            VInt_pushBack( idxResR[idzR][tj],nColC[tj][idzR] );
                            VDbl_pushBack( coefResR[idzR][tj], use );
                            nColC[tj][idzR]++;
                    }
                }
            }
        }
    }

    LinearProgram *mipCut;
    ALLOCATE_VECTOR(ccr->cutElem,VecInt**,Inst_nResR(inst));
    ALLOCATE_VECTOR(ccr->cutCoef,VecDbl**,Inst_nResR(inst));
    int *nCutR;
    ALLOCATE_VECTOR_INI(nCutR,int,Inst_nResR(inst));


    for(int r = 0 ; r < Inst_nResR(inst); r++ ) {
        ALLOCATE_VECTOR(ccr->cutElem[r],VecInt*,nTimes);
        ALLOCATE_VECTOR(ccr->cutCoef[r],VecDbl*,nTimes);
        for(int nc = 0 ; nc < nTimes; nc++ ) {
            ccr->cutElem[r][nc] = VInt_create();
            ccr->cutCoef[r][nc] = VDbl_create();
            ccr->nAllocations[r]++;
        }
    }

    for ( int t=0 ; t<nTimes; ++t ) {
        for ( int r=0 ; (r<Inst_nResR(inst)) ; ++r ) {

            if(nColC[t][r]==0) continue;
            //  if(continuous) {
            //  _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
            // if(_time < 1) {
            //    printf( "Mip cut RR time is over %f \n", _time);
            //    goto FREE_MEMORY;
            // }
            // }
           // printf("creating mipCut RR: start\n");fflush(stdout);
            mipCut = lp_create();
            lp_set_print_messages(mipCut,0);
            //lp_set_max_nodes(mipCut, 500);
            //lp_set_max_solutions(mipCut,10);
            lp_add_bin_cols( mipCut, nColC[t][r], obj[t][r], VStr_ptr(names[t][r]));

            sprintf( name, "resR(%d,%d)", r, t );
            sprintf( name2, "ineqValid(%d,%d)", r, t );

            const VecInt *vidz = idxResR[r][t];
            const VecDbl *vcoef = coefResR[r][t];

            const int nz = VInt_size( vidz );
            int *idz =  VInt_getPtr( (VecInt *) vidz );
            double *coef =  VDbl_getPtr( (VecDbl *) vcoef );

            int nColsLP = lp_cols(lp);

            if (nz) {


                lp_add_row( mipCut, nz, idz, coef, name, 'G', Inst_capResR(inst, r)+1 );

                lp_add_row( mipCut, nColC[t][r], zx[t][r], obj[t][r], name2, 'L', 0.9995 );
                //lp_write_lp(mipCut,"lp_cut.lp");

                //lp_set_concurrentMIP(mipCut,1);
                //lp_set_method(mipCut,4);
                //lp_set_seed(mipCut,100000);
                int st = lp_optimize(mipCut);



                double valueZ = lp_obj_value(mipCut);
                //printf("creating mipCut RR: end\n");fflush(stdout);
                if(fabs(valueZ) < 0.9995 && st==0) {


                //   printf("find real cuts mipCut RR to resource %d on time %d : start\n", r,t);fflush(stdout);
                  //lp_write_lp(mipCut,"lp_cut.lp");
                    //printf(" writing mipCut st %d end \n", st);
                   // printf("valueZ %g t %d : ", valueZ,t);
                    //getchar();

                    IntDblPair *cutIdxCo;
                    ALLOCATE_VECTOR(cutIdxCo, IntDblPair, nColsLP);

                    const double *zf = lp_x(mipCut);

                    // filling x from fractional solution
                    int c=0, nCutReal =0, cs =0 ;
                    IntTriple *varCut;
                    ALLOCATE_VECTOR_INI(varCut,IntTriple, nJobs*Inst_nMaxModes(inst)*nTimes)
                    IntTriple *varCutOrig;
                    ALLOCATE_VECTOR_INI(varCutOrig,IntTriple, nJobs*Inst_nMaxModes(inst)*nTimes)

                    // int cZidx[nJobs*nModes*nTimes];

                    int nc = lp_cols(mipCut);

                    //   for(int ia = 0 ; ia <nz ; ia++){
                    //       printf("%f * %d (z*[%d] %f) + ", coef[ia], idz[ia], idz[ia], zf[idz[ia]]);
                    //   }
                    //  printf(" < %d \n", Inst_capResR(inst, r)+1);

                    double slack = 0;

                  //  printf("\n t %d: ", t);
                    for ( i=0 ; (i<nc) ; ++i ) {
                        char nameC[STR_SIZE];

                        lp_col_name( mipCut, i, nameC );
                        if (tolower(nameC[0]=='z')) {
                            parseName( nameC, prefix, auxidZ );

                            if (zf[i]>EPS) {
                                //    printf(" zf[%d] %f > EPS %f ", i, zf[i], EPS);
                                if (prefix[0]=='z') {

                                    int j = auxidZ[0];
                                    int m = auxidZ[1];
                                    int tt = auxidZ[2];
                                    const Job *jobOnCut = Inst_job(inst,j);
                                    const Mode* modeOnCut = Job_mode(jobOnCut,m);
                                    if(Mode_duration(modeOnCut)==0) continue;

                                    assert(j<nJobs);
                                    assert(m<nModes);
                                    assert(tt<nTimes);
                                    assert(c<nJobs*nModes*nTimes);
                                    assert( xIdx[j][m][tt] != -1);

                                    // cZidx[c] = xIdx[j][m][tt];
                                    varCut[c].j = j;
                                    varCut[c].m = m;
                                    varCut[c].t = tt;
                                    varCut[c].value = zf[i];

                                    varCutOrig[cs].j = j;
                                    varCutOrig[cs].m = m;
                                    varCutOrig[cs].t = tt;
                                    varCutOrig[cs].value = zf[i];

                                    cutIdxCo[c].a =  xIdx[j][m][tt];//cZidx[c];
                                    cutIdxCo[c].b = 1.0;
                                    slack += xf[cutIdxCo[c].a]*cutIdxCo[c].b;


                                    nCutReal++;
                                    c++;
                                    cs++;
                                //    printf("O (%d,%d,%d), ", j,m,tt);

                                }
                            }
                        }
                    }



                      // printf("slack-(cs-1) %f\n", slack-(cs-1));

                        char nameCut[STR_SIZE];
                        //printf( "cutRR(%d,%d)#%d\n", r, t, nCut );
                        //   if(r==1 && t==5 && nCut == 1274) getchar();
                        //   if(r==1 && t==6 && nCut == 1588) getchar();
                        sprintf( nameCut, "cutRR(%d,%d)#%d#%d", r, t, nround, nCut );
                        VStr_pushBack(ccr->cutname[r],nameCut);
                        // double coe[c];
//                        getchar();
                        CutP_quick_sort(cutIdxCo, c);
                        for(int sel  = 0 ; sel <c; sel++) {
                            int elem =cutIdxCo[sel].a;
                            //  cZidx[sel] = elem;
                            // coe[sel] = cutIdxCo[sel].b;

                            VInt_pushBack(ccr->cutElem[r][nCutR[r]],elem);
                            VDbl_pushBack(ccr->cutCoef[r][nCutR[r]],cutIdxCo[sel].b);
                        }
                        VDbl_pushBack(ccr->cutrhs[r],nCutReal-1);
                        VInt_pushBack(ccr->cutsense[r],0);
                        VInt_pushBack(ccr->cutdominated[r],0);
                        VInt_pushBack(ccr->cutnelem[r],c);
                        VDbl_pushBack(ccr->cutviolation[r],slack-(nCutReal-1));

                        //double newrhs = CutP_model_lift( inst, VStr_get(ccr->cutname[r],nCutR[r]), VInt_get(ccr->cutnelem[r],nCutR[r]),VInt_getPtr(ccr->cutElem[r][nCutR[r]]), VDbl_getPtr(ccr->cutCoef[r][nCutR[r]]), VDbl_get(ccr->cutrhs[r],nCutR[r]), lp, _time);
                        //ir(newrhs!=-1)
                          //  VDbl_set(ccr->cutrhs,r,newrhs);
                        //   int *cutIdx = VInt_getPtr(ccr->cutElem[nCut]);
                        // double *cutCoef = VDbl_getPtr(ccr->cutCoef[nCut]);
                        // CutP_quick_sort_vec(cutIdx,cutCoef, c);
                        nCutR[r]++;

                        nCut++;

                        free(cutIdxCo);

                        free(varCut);
                        free(varCutOrig);
                }
            }
            lp_free(&mipCut);
        }
    }


    free(nCutR);
    //    double timeinitdomresource = omp_get_wtime();

    for(int r = 0 ; r < Inst_nResR(inst); r++ ) {
        // printf("VInt_size(ccr->cutrhs[r]) %d: %d\n",VInt_size(ccr->cutrhs[r]) ,  VInt_size(ccr->cutrhs[r]) *VInt_size(ccr->cutrhs[r]) );
        for(int ncdA = 0; ncdA < VDbl_size(ccr->cutrhs[r]) ; ncdA++) {
            for(int ncdB = ncdA+1; ncdB < VDbl_size(ccr->cutrhs[r]) ; ncdB++) {
                if(VInt_get(ccr->cutdominated[r],ncdA) == 1) break;
                if(VInt_get(ccr->cutdominated[r],ncdB) == 1) continue;
                CutP_dominanceBetweenTwo( lp, ncdA, VInt_getPtr(ccr->cutElem[r][ncdA]), VDbl_getPtr(ccr->cutCoef[r][ncdA]), VDbl_get(ccr->cutrhs[r],ncdA), VInt_get(ccr->cutsense[r],ncdA), VInt_get(ccr->cutnelem[r],ncdA), ncdB, VInt_getPtr(ccr->cutElem[r][ncdB]), VDbl_getPtr(ccr->cutCoef[r][ncdB]), VDbl_get(ccr->cutrhs[r],ncdB), VInt_get(ccr->cutsense[r],ncdB), VInt_get(ccr->cutnelem[r],ncdB),  ccr->cutdominated[r]);
                // getchar();
            }
        }
    }
    //    printf("time compute dominance resource %f\n",omp_get_wtime()-timeinitdomresource);

    for ( int i=0 ; (i<Inst_nResR( inst )*(nTimes)) ; ++i ) {
        VInt_free( &idxResR[0][i] );
        VDbl_free( &coefResR[0][i] );
    }
    free( idxResR[0] );
    free( coefResR[0] );
    free( idxResR );
    free( coefResR );


    for(int t  = 0 ; t < nTimes ; t++) {

        for(int r  = 0 ; r < Inst_nResR(inst) ; r++) {
            VStr_free(&names[t][r]);
            for(int j = 0 ; j < nJobs ; j++) {
               for(int m = 0 ; m <nModes ; m++){
                   free(zIdx[t][r][j][m]);
                   free(addedTRJM[t][r][j][m]);
                }
                free(zIdx[t][r][j]);
                free(addedTRJM[t][r][j]);
            }
            free(zIdx[t][r]);
            free(addedTRJM[t][r]);
            free(obj[t][r]);
            free(zx[t][r]);
        }
        free(zIdx[t]);
        free(addedTRJM[t]);
        free(names[t]);
        free(nColC[t]);
        free(obj[t]);
        free(zx[t]);
    }
    free( nColC);
    free(obj);
    free(zx);
    free(names);
    free(zIdx);
    free(addedTRJM);
    free(xIdx[0][0]);
    free(xIdx[0]);
    free(xIdx);
    free(x[0][0]);
    free(x[0]);
    free(x);

    //    ccr->timeseparation = (double) omp_get_wtime()-tinit;

#undef MAX_IDX
}

/*free memory*/
void CutC_free( CutC **_cutC )
{

    CutC *cutC = *_cutC;
    for(int r = 0 ; r < Inst_nResR(cutC->inst); r++) {
        //      printf(" cutC->nAllocations[r] %d size %d",  cutC->nAllocations[r], VDbl_size(cutC->cutrhs[r]));fflush(stdout);
        for(int nr = 0 ; nr < cutC->nAllocations[r]; nr++) {
            VInt_free(&cutC->cutElem[r][nr]);
            VDbl_free(&cutC->cutCoef[r][nr]);
        }
        VStr_free(&cutC->cutname[r]);
        VDbl_free(&cutC->cutrhs[r]);
        VDbl_free(&cutC->cutviolation[r]);
        VInt_free(&cutC->cutnelem[r]);
        VInt_free(&cutC->cutdominated[r]);
        VInt_free(&cutC->cutsense[r]);
        if( cutC->nAllocations[r]!=0) {
            free(cutC->cutElem[r]);
            free(cutC->cutCoef[r]);
        }
    }
    free(cutC->nAllocations);
    free(cutC->cutElem);
    free(cutC->cutCoef);

    free(cutC->cutname);
    free(cutC->cutrhs);
    free(cutC->cutviolation);
    free(cutC->cutnelem);
    free(cutC->cutdominated);
    free(cutC->cutsense);

    free( cutC );
    *_cutC = NULL;

}


