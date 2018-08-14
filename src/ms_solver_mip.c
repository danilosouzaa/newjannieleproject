
/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * ModeSet MIP Solver
 * formulation variables:
 *     y_j   starting time of mode j
 *     d_j   duration of job j
 *     s_jm  if mode m is selected for job j
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include <omp.h>
#include "ms_solver_mip.h"
#include "lp.h"
#include "macros.h"
#include "vec_str.h"
#include "mode_set.h"
#define STR_SIZE        256

const double W_END_PROJECT = 500.0;
const double W_DURATION    =   1.0;

struct _MSM_Solver {
    const Instance *inst;
    LinearProgram *mip;

    ModeSet *modeSet;

    /* indexes of variables  */
    int *yIdx; /* index of variable of job i */
    int *dIdx;
    int **sIdx;

    int firstJob;
    int lastJob;
};

static void MSM_getSol( MSM_Solver *msm );

MSM_Solver *MSM_create( const Instance *inst, int firstJob, int lastJob, double timeLeft)
{
    double _time = 0;
    double startT = omp_get_wtime();

    assert( inst && firstJob >= 0 && lastJob >= firstJob && lastJob < Inst_nJobs(inst) );
    const int maxModes = 3;//Inst_nMaxModes( inst );
    const int nVars = Inst_nJobs(inst)*2+Inst_nJobs(inst)*maxModes;
    double lb[nVars];
    double ub[nVars];
    double obj[nVars];
    char integer[nVars];
    char name[STR_SIZE];
    //ALLOCATE_VECTOR(name, char, STR_SIZE);

    VecStr *names = VStr_create( STR_SIZE );

    MSM_Solver *mss;
    ALLOCATE( mss, MSM_Solver );

    mss->inst = inst;
    mss->firstJob = firstJob;
    mss->lastJob = lastJob;

    ALLOCATE_VECTOR(  mss->yIdx, int, Inst_nJobs(inst) );
    ALLOCATE_VECTOR(  mss->dIdx, int, Inst_nJobs(inst) );
    ALLOCATE_VECTOR(  mss->sIdx, int*, Inst_nJobs(inst) );
    ALLOCATE_VECTOR(  mss->sIdx[0], int, Inst_nJobs(inst)*maxModes );
    for ( int i=1 ; (i<Inst_nJobs(inst)) ; ++i )
        mss->sIdx[i] = mss->sIdx[i-1] + maxModes;

    //ALLOCATE_VECTOR(  mss->xIdx[0], int, Inst_nJobs(inst)*maxModes );

    _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
    if(_time <= 0) {
        printf( "\nMSM Create: Time is over %f Time Left %f \n", _time, timeLeft);
        exit(0);
    }

    LinearProgram *mip = lp_create();
    int j=0;
    for ( int i=firstJob ; (i<=lastJob) ; ++i,++j ) {
        sprintf( name, "y(%d)", i );
        VStr_pushBack( names, name );
        lb[j] = 0.0;
        ub[j] = DBL_MAX;
        const Job *job = Inst_job(inst,i);
        char endProjDummy = ((Job_maxDuration(job)==0) && (Job_nPred(job)>=1));
        obj[j] = endProjDummy ? W_END_PROJECT : 0.0;
        integer[j] = True;
        mss->yIdx[i] = j;
        _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
        if(_time <= 0) {
            printf( "\nMSM Create: Time is over %f Time Left %f \n", _time, timeLeft);
            exit(0);
        }

    }
    for ( int i=firstJob ; (i<=lastJob) ; ++i,++j ) {
        sprintf( name, "d(%d)", i );
        VStr_pushBack( names, name );
        const Job *job = Inst_job( inst, i );
        lb[j] = Job_minDuration( job );
        ub[j] = Job_maxDuration( job );
        obj[j] = W_DURATION;
        integer[j] = True;
        mss->dIdx[i] = j;
        _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
        if(_time <= 0) {
            printf( "\nMSM Create: Time is over %f Time Left %f \n", _time, timeLeft);
            exit(0);
        }
    }

    for ( int i=firstJob ; (i<=lastJob) ; ++i ) {
        // for ( int i=firstJob+1 ; (i<lastJob) ; ++i ) {            //janniele
        const Job *job = Inst_job( inst, i );
        _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
        if(_time <= 0) {
            printf( "\nMSM Create: Time is over %f Time Left %f \n", _time, timeLeft);
            exit(0);
        }
        for ( int m=0 ; (m<Job_nModes(job)) ; ++m,++j ) {
            const Mode *mode = Job_mode( job, m );
            lb[j] = 0.0;
            ub[j] = Mode_isFeasible( inst, mode ) ? 1.0 : 0.0;
            obj[j] = 0.0;
            integer[j] = True;
            sprintf( name, "s(%d,%d)", i, m );
            VStr_pushBack( names, name );
            mss->sIdx[i][m] = j;
        }
        _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
        if(_time <= 0) {
            printf( "\nMSM Create: Time is over %f Time Left %f \n", _time, timeLeft);
            exit(0);
        }

    }

    lp_add_cols( mip, j, obj, lb, ub, integer, VStr_ptr(names) );

    {
        int idx[maxModes+1];
        double coef[maxModes+1];
        // sel mode constraints
        for ( int i=firstJob ; (i<=lastJob) ; ++i ) {
            const Job *job = Inst_job( inst, i );
            for ( int m=0 ; (m<Job_nModes(job)) ; ++m ) {
                idx[m] = mss->sIdx[i][m];
                coef[m] = 1.0;
            }
            sprintf( name, "selMode(%d)", i );
            _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
            if(_time <= 0) {
                printf( "\nMSM Create: Time is over %f Time Left %f \n", _time, timeLeft);
                exit(0);
            }
            lp_add_row( mip, Job_nModes(job), idx, coef, name, 'E', 1.0 );
        }



        // duration constraints
        for ( int i=firstJob ; (i<lastJob) ; ++i ) {
            const Job *job = Inst_job( inst, i );
            if(Job_nModes(job)==1&&Job_minDuration(job)==0) continue; //janniele
            for ( int m=0 ; (m<Job_nModes(job)) ; ++m ) {
                const Mode *mode = Job_mode( job, m );
                idx[m] = mss->sIdx[i][m];
                coef[m] = Mode_duration(mode);
            }
            idx[Job_nModes(job)] = mss->dIdx[i];
            coef[Job_nModes(job)] = -1.0;
            sprintf( name, "duration(%d)", i );
            _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
            if(_time <= 0) {
                printf( "\nMSM Create: Time is over %f Time Left %f \n", _time, timeLeft);
                exit(0);
            }
            lp_add_row( mip, Job_nModes(job)+1, idx, coef, name, 'E', 0.0 );
        }

        // precedence constraints
        for ( int i=firstJob ; (i<=lastJob) ; ++i ) {
            const Job *job = Inst_job( inst, i );
            idx[0] = mss->yIdx[i];
            coef[0] = 1.0;
            for ( int p=0 ; (p<Job_nPred(job)) ; ++p ) {
                const int ip = Job_pred( job, p );
                idx[1] = mss->yIdx[ip];
                coef[1] = -1.0;
                idx[2] = mss->dIdx[ip];
                coef[2] = -1.0;

                sprintf( name, "pred(%d,%d)", i, ip );
                _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
                if(_time <= 0) {
                    printf( "\nMSM Create: Time is over %f Time Left %f \n", _time, timeLeft);
                    exit(0);
                }
                lp_add_row( mip, 3, idx, coef, name, 'G', 0.0 );
            }
        }

    }

    /* non-renewable resources constraint */
    {
        int  maxNZ = Inst_nJobs(inst)*maxModes;
        int idx[maxNZ];
        double coef[maxNZ];

        for ( int k=0 ; (k<Inst_nResN(inst)) ; ++k ) {
            int nz = 0;
            for ( int i=firstJob ; (i<=lastJob) ; ++i ) {
                const Job *job = Inst_job( inst, i );
                for ( int m=0 ; (m<Job_nModes(job)) ; ++m ) {
                    const Mode *mode = Job_mode( job, m );
                    for ( int l=0 ; l<Mode_nResN(mode) ; ++l ) {
                        if (Mode_idxResN(mode,l)==k) {
                            idx[nz] = mss->sIdx[i][m];
                            coef[nz] = Mode_useResN(mode,l);
                            ++nz;
                        }
                        _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
                        if(_time <= 0) {
                            printf( "\nMSM Create: Time is over %f Time Left %f \n", _time, timeLeft);
                            exit(0);
                        }
                    }
                }
            }
            if (nz) {
                sprintf(name,"useN(%d)",k);
                lp_add_row( mip, nz, idx, coef, name, 'L', Inst_capResN(inst,k) );
            }
        } // all resources
    }

    mss->mip = mip;
    mss->modeSet = Modes_create( inst );

    VStr_free( &names );

    //lp_write_lp(mip,"modes");

    return mss;
}

int MSM_solve( MSM_Solver *msm, double timeLeft)
{
    double startT = omp_get_wtime();

    double _time;

    LinearProgram *mip = msm->mip;
    lp_set_print_messages( mip, 0 );

    _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
    if(_time < 1) {
        printf( "Mip ModeMin Solve time is over %f \n", _time);
        exit(0);
    }

    lp_set_max_seconds(mip,_time);

    //lp_set_concurrentMIP(mip,1);
    //lp_set_method(mip,4);
    //lp_set_seed(mip,100000);
    int status =  lp_optimize( mip );
    // _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
    if(_time < 1) {
        printf( "Mip ModeMin Solve time is over %f \n", _time);
        exit(0);
    }

#ifdef DEBUG
    printf("status: %d (LP_OPTIMAL: %d, LP_FEASIBLE: %d, LP_INFEASIBLE: %d)  \n", status, LP_OPTIMAL, LP_FEASIBLE, LP_INFEASIBLE);
#endif // DEBUG
    //assert( status==LP_OPTIMAL || status==LP_FEASIBLE );
    return status;
    //MSM_getSol( msm );
}

LinearProgram *MSM_lp(const MSM_Solver *msm)
{
    return msm->mip;
}

void MSM_solve_as_continuous( MSM_Solver *msm )
{
    LinearProgram *mip = msm->mip;
    lp_set_print_messages( mip, 1 );

    //lp_set_concurrentMIP(mip,1);
    //lp_set_method(mip,4);
    //lp_set_seed(mip,100000);
    int status = lp_optimize_as_continuous(mip);
#ifdef DEBUG
    printf("status: %d (LP_OPTIMAL: %d, LP_FEASIBLE: %d, LP_INFEASIBLE: %d)  \n", status, LP_OPTIMAL, LP_FEASIBLE, LP_INFEASIBLE);
#endif // DEBUG
    assert( status==LP_OPTIMAL || status==LP_FEASIBLE );
    MSM_getSol( msm );
}

const ModeSet *MSM_modes( const MSM_Solver *msm )
{
    return msm->modeSet;
}

void MSM_free( MSM_Solver **_msm )
{
    MSM_Solver *msm = *_msm;
    free( msm->yIdx );
    free( msm->dIdx );

    free( msm->sIdx[0] );
    free( msm->sIdx );
    Modes_free( &msm->modeSet );
    lp_free( &msm->mip );

    free( msm );
    *_msm = NULL;
}

char MSM_changeModes( MSM_Solver *msm, const ModeSet *current, int minChanges,  int maxChanges, const int **residency )
{
    const Instance *inst = msm->inst;
    int minR = INT_MAX;
    int maxR = 0;
    int maxD = 0;
    int maxObj = 0;
    for ( int i=msm->firstJob ; (i<=msm->lastJob) ; ++i ) {
        const Job *job = Inst_job( inst, i );
        for ( int m=0 ; (m<Job_nModes(job)) ; ++m ) {
            minR = MIN( minR, residency[i][m] );
            maxR = MAX( maxR, residency[i][m] );
            const Mode *mode = Job_mode( job, m );
            maxD = MAX( maxD, Mode_duration(mode) );
            maxObj++;
        }
    }

    int idx[maxObj];
    double obj[maxObj];
    double coef[maxObj];
    int nCh = 0;
    double resRange = maxR-minR;

    /* changing weights in modes variables */
    for ( int i=msm->firstJob ; (i<=msm->lastJob) ; ++i ) {
        const Job *job = Inst_job( inst, i );
        for ( int m=0 ; (m<Job_nModes(job)) ; ++m ) {
            idx[nCh] = msm->sIdx[i][m];
            if (fabs(resRange)<1e-10)
                obj[nCh] = 0.0;
            else
                obj[nCh] = ((double)maxD*2.0)*(((double)(residency[i][m]-minR))/resRange);
            ++nCh;
        }
    }
    lp_chg_obj( msm->mip, nCh, idx, obj );

    /* removing old rows for change */
    for ( int i=lp_rows(msm->mip)-1 ; i>=0 ; --i ) {
        char rName[STR_SIZE];
        lp_row_name( msm->mip, i, rName );
        if (strstr( rName, "Change" ))
            lp_remove_row( msm->mip,i );
        else
            break;
    }

    /* row with minimum mode changes */
    int nz = 0;
    for ( int i=msm->firstJob ; (i<=msm->lastJob) ; ++i ) {
        const Job *job = Inst_job( inst, i );
        for ( int m=0 ; (m<Job_nModes(job)) ; ++m ) {
            if ( Modes_job( msm->modeSet, i ) != m ) {
                idx[nz] = msm->sIdx[i][m];
                coef[nz] = 1.0;
                ++nz;
            }
        }
    }
    char name[STR_SIZE];
    sprintf( name, "minChange");
    lp_add_row( msm->mip, nz, idx, coef, name, 'G', minChanges );


    //lp_write_lp(msm->mip,"modesCH");

    //lp_set_concurrentMIP(msm->mip,1);
    //lp_set_method(msm->mip,4);
    //lp_set_seed( msm->mip,100000);
    int status = lp_optimize( msm->mip );

    if ((status==LP_OPTIMAL || status==LP_FEASIBLE)) {
        MSM_getSol( msm );
        return True;
    }



    return False;
}


int MSM_getSIdx( const MSM_Solver *msm, int j, int m)
{

    assert(msm != NULL);
    assert(j>= 0 && j < Inst_nJobs(msm->inst));
    //sprintf( "s(%d,%d) %d", j, m, msm->sIdx[j][m] );

    return msm->sIdx[j][m];

}

void MSM_getSol( MSM_Solver *msm )
{
    const Instance *inst = msm->inst;
    LinearProgram *mip = msm->mip;

    const double *x = lp_x( mip );
    int j=0,c=0, aux = 0;

    for ( int i=msm->firstJob ; (i<=msm->lastJob) ; ++i ) {

        const Job *job = Inst_job(inst,i);
        for ( int m=0 ; (m<Job_nModes(job)) ; ++m ) {
            if ( x[ msm->sIdx[i][m] ] >= 0.98 ) {
                Modes_modify( msm->modeSet, i, m);
                j++;
                aux = 1;
            }
        }
        if(aux==0) {
            if ( x[ msm->yIdx[i] ] >= 0.98 )
                c++;
        }
        aux = 0;
    }
}

