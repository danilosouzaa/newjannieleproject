/* GOAL PSP Solver II
 * Resource-Constrained Project Scheduling Problem (RCPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G.
 */

#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "instance.h"
#include "macros.h"
#include "vec_int.h"
#include "vec_str.h"
#include "tokenizer.h"
#include "str_utils.h"
#include "stack.h"

typedef struct {
    int a;
    int b;
} IntPair;

int cmp_int_pair_b( const void *v1, const void *v2 )
{
    const IntPair *ip1 = (const IntPair *) v1;
    const IntPair *ip2 = (const IntPair *) v2;

    if ( ip1->b!=ip2->b )
        return ip1->b-ip2->b;

    return ip1->a-ip2->a;
}

struct _Mode {
    int index;

    int duration;

    /* renewable resources */
    int nResR;
    int *idxResR;
    int *useResR;

    /* non-renewable resources */
    int nResN;
    int *idxResN;
    int *useResN;
};


struct _Job {
    int index;
    int idxOnProject;
    int idxProject;

    int nModes;
    int nInfeasModes;
    int infeaseWithAll;
    Mode *modes;

    int nSucc;
    int *succ;

    int nPred;
    int *pred;

    IntPair *idxModesSort;

    /*improving search of pred and succs direct and indirect.
    Janniele*/

    int *hasPred;
    int *hasSucc;

    //int *indPred;
    //int *indSucc;

    int *hasIndPred;
    int *hasIndSucc;

    /*end*/

    int minDuration;
    int idxModeMinDuration;
    int maxDuration;

    int est;

    /* if it is an instance job, says to which
       project job it relates */
    const struct _Job *origJob;
};

struct _Project {
    int index;
    int releaseDate;
    int criticalPath;

    int dueDate;
    int tardCost;
    int mpmTime;

    /* idx of the first project job in the instance */
    int idxFirstJob;

    int nJobs;
    Job *jobs;

};

struct _Instance {
    int nProjects;

    Project *projects;
    /* resource capacity information */
    int nResR;
    int *capR;

    int nResRGlobal;

    int nResN;
    int *capN;

    int nJobs;
    Job *jobs;

    int *idxIniResNProj;
    int sumTPD;
    int sumTMS;

    int **matMaxD;

    int ***matMaxDJM;

    int **maxDurationPath;
    int ****maxDurationPathInterMode;

   //  int ****unitsResourceES;

    int **jobmodeinfeasible;

    VecInt **paths;

    int nMaxModes;
};


static void Inst_readProject( const char *dir, const char *projectFile, int projIdx, Instance *inst, int nres, const int globalRCap[], VecInt *capR, VecInt *capNR );
static void Inst_readProject_rcpsp( const char *dir, const char *projectFile, int projIdx, Instance *inst, VecInt *capR, VecInt *capNR );
static void Inst_readProject_mrcpsp( const char *dir, const char *projectFile, int projIdx, Instance *inst, VecInt *capR, VecInt *capNR );
static void Inst_computeCompPaths( const Instance *inst, VecInt **paths);


static void Inst_freeJobVectorContents(Instance * inst, Job *jStart, Job *jEnd );


/* copies mode contents from one mode to another */
static void Mode_cpy( Mode *mTarget, const Mode *mSource);

/* copies job contents */
static void Job_cpy( const Instance *inst, Job *jTarget, const Job* jSource, int* nRemovedModes );

static void Inst_computeEST( Instance *inst );

static void Inst_computeESTJob( Instance *inst, int jIdx );

/* reads a series of integers from a string and returns the number of elements readed
   aborts if number of elements > maxElements */
int readIntVectorStr( const char *str, Tokenizer *tok, int minElements, int maxElements, int el[] );


void Project_rcpsp_print(Instance *inst)
{
    printf("\n************************************************************************");
    printf("\nPRECEDENCE RELATIONS:");
    printf("\njobnr.\t#modes\t#successors\tsuccessors");

    for(int i=0; i<inst->projects[0].nJobs; ++i) {
        printf("\n%d\t%d\t%d\t\t", i+1, 1, inst->projects[0].jobs[i].nSucc);
        for(int j=0; j<inst->projects[0].jobs[i].nSucc; ++j)
            printf("%d\t", inst->projects[0].jobs[i].succ[j]+1 );
    }

    printf("\n************************************************************************");
    printf("\nREQUESTS/DURATIONS:");
    printf("\njobnr.\tmode\tduration\tR 1\tR 2\tR 3\tR 4");

    for(int i=0; i<inst->projects[0].nJobs; ++i) {
        printf("\n%d\t%d\t%d\t\t", i+1, 1, inst->projects[0].jobs[i].modes[0].duration);
        //for(int j=0; j<inst->projects[0].jobs[i].modes[i].nResR; ++j)
        //printf("%d\t", inst->projects[0].jobs[i].modes[0].useResR[0] );


    }


}

void Project_mrcpsp_print(Instance *inst)
{
    printf("\n************************************************************************");
    printf("\nPRECEDENCE RELATIONS:");
    printf("\njobnr.\t#modes\t#successors\tsuccessors");

    for(int i=0; i<inst->projects[0].nJobs; ++i) {
        printf("\n%d\t%d\t%d\t\t", i+1, 1, inst->projects[0].jobs[i].nSucc);
        for(int j=0; j<inst->projects[0].jobs[i].nSucc; ++j)
            printf("%d\t", inst->projects[0].jobs[i].succ[j]+1 );
    }

    printf("\n************************************************************************");
    printf("\nREQUESTS/DURATIONS:");
    printf("\njobnr.\tmode\tduration\tR 1\tR 2\tR 3\tR 4");

    for(int i=0; i<inst->projects[0].nJobs; ++i) {
        printf("\n%d", i+1);
        for(int j=0; j<inst->projects[0].jobs[i].nModes; ++j)
            printf("\t%d\t%d\t\t\n", j+1, inst->projects[0].jobs[i].modes[j].duration);
    }

}

/* loads instance from a file RCPSP */
Instance *Inst_create_rcpsp( const char *dir, const char *instance )
{
    char fname[FILE_NAME_SIZE];
    sprintf( fname, "%s/%s", dir, instance );
    int nRemovedModes = 0;

    Instance *inst;
    ALLOCATE_INI( inst, Instance );

    inst->nProjects = 1;

    ALLOCATE_VECTOR( inst->projects, Project, inst->nProjects );
    ALLOCATE_VECTOR_INI( inst->idxIniResNProj, int, inst->nProjects );



    int totalJobs = 0; //0
    inst->projects[0].index = 0;
    /* current index for renewableResources */

    VecInt *tcapR = VInt_create();
    VecInt *tcapNR = VInt_create();

    inst->projects[0].index = 0;

    Inst_readProject_rcpsp( dir, instance, 0, inst, tcapR, tcapNR );
    totalJobs = inst->projects[0].nJobs;

    inst->nResR = VInt_size( tcapR );
    inst->nResN = VInt_size( tcapNR );

    ALLOCATE_VECTOR_INI( inst->capR, int, VInt_size( tcapR ) );
    memcpy( inst->capR, VInt_getPtr(tcapR), sizeof(int)*VInt_size( tcapR ) );

    ALLOCATE_VECTOR_INI( inst->capN, int, VInt_size( tcapNR) );
    memcpy( inst->capN, VInt_getPtr(tcapNR), sizeof(int)*VInt_size( tcapNR ) );

    VInt_free( &tcapR );
    VInt_free( &tcapNR );

    /* adding global jobs */
    inst->nJobs = totalJobs;
    ALLOCATE_VECTOR_INI( inst->jobs, Job, totalJobs );

    int jIdx = 0;
    Job *job = inst->jobs;

    /* temporarily used as counters */
    /* adding project jobs */

    const int firstJobP = jIdx;
    Project *p = &inst->projects[0];
    p->idxFirstJob = jIdx;

    for ( int j=0; j<p->nJobs; ++j,jIdx++,job++ ) {

        job->index = jIdx;

        const Job *projJob = &p->jobs[j];
        Job_cpy( inst, job, projJob, &nRemovedModes );

        job->index = jIdx;
        job->idxOnProject = projJob->index;
        job->idxProject = p->index;
        job->origJob = projJob;

#ifdef DEBUG
        // projects must have dummy jobs too
        if ( j == 0 || j==p->nJobs-1  ) {
            assert( projJob->nModes==1 );
            assert( projJob->modes[0].duration == 0 );
        }
#endif

        ALLOCATE_VECTOR( job->hasPred, int, inst->nJobs );
        ALLOCATE_VECTOR( job->hasSucc, int, inst->nJobs );
        ALLOCATE_VECTOR_INI( job->hasIndPred, int, inst->nJobs );
        ALLOCATE_VECTOR_INI( job->hasIndSucc, int, inst->nJobs );

        // updating indexes
        for ( int k=0; k<job->nSucc; ++k ) {
            job->succ[k] += firstJobP;
            job->hasSucc[job->succ[k]] = 1; // 1 to indicate if a job is successor Janniele
        }

        for ( int k=0; k<job->nPred; ++k ) {
            job->pred[k] += firstJobP;
            job->hasPred[job->pred[k]] = 1; // 1 to indicate if a job is pred Janniele
        }

    } // all jobs

    /* updating indirect preds */

    Stack *stackJobs = Stk_create(inst->nJobs);

    for(int jj = 0; jj < inst->nJobs; jj++) {

        const Job* job = Inst_job(inst, jj);

        int idxJob = Job_index(job);
        Stk_push(stackJobs, idxJob);

        while (!Stk_isEmpty(stackJobs)) {
            idxJob = Stk_pop(stackJobs);

            const Job *jobInd = Inst_job(inst,idxJob);

            for (int pred = 0; pred < Job_nPred(jobInd); pred++) {
                if (!job->hasIndPred[Job_pred(jobInd,pred)]) {
                    job->hasIndPred[Job_pred(jobInd,pred)] = 1;
                    Stk_push(stackJobs,Job_pred(jobInd,pred));
                    //                    printf("%d: %d \n", Job_pred(jobInd,pred), job->hasIndPred[Job_pred(jobInd,pred)] );
                }
            }
        }

        /* updating indirect preds*/
        idxJob = Job_index(job);
        Stk_push(stackJobs, idxJob);

        while (!Stk_isEmpty(stackJobs)) {
            idxJob = Stk_pop(stackJobs);

            const Job *jobInd = Inst_job(inst,idxJob);
            //printf("Succs:\n");
            for (int succ = 0; succ < Job_nSucc(jobInd); succ++) {
                if (!job->hasIndSucc[Job_succ(jobInd, succ)]) {
                    job->hasIndSucc[Job_succ(jobInd, succ)] = 1;
                    Stk_push(stackJobs,Job_succ(jobInd, succ));
                    //      printf("%d: %d \n", Job_succ(jobInd,succ), job->hasIndSucc[Job_succ(jobInd,succ)] );
                }
            }
        }


    }

    Stk_free(&stackJobs);

    /* maximum and minimum duration for a job */
    ALLOCATE_VECTOR_INI(inst->maxDurationPath, int*, inst->nJobs);
    ALLOCATE_VECTOR(inst->maxDurationPathInterMode, int***, inst->nJobs);

    for ( int i=0; i<inst->nJobs; ++i ) {
        Job *job = inst->jobs + i;
        job->minDuration = INT_MAX;
        job->idxModeMinDuration = 0;
        job->maxDuration = 0;
        ALLOCATE_VECTOR(job->idxModesSort, IntPair, job->nModes);
        ALLOCATE_VECTOR_INI(inst->maxDurationPath[i], int, job->nModes);
        ALLOCATE_VECTOR(inst->maxDurationPathInterMode[i], int**, job->nModes);
        for ( int m=0; m<job->nModes; ++m ) {
            ALLOCATE_VECTOR(inst->maxDurationPathInterMode[i][m], int*, inst->nJobs);
            job->idxModesSort[m].a = m;
            job->idxModesSort[m].b = job->modes[m].duration;
            if( job->modes[m].duration <  job->minDuration)
                job->idxModeMinDuration = m;
            job->minDuration = MIN( job->minDuration, job->modes[m].duration );
            job->maxDuration = MAX( job->maxDuration, job->modes[m].duration );
            for(int j = 0 ; j < inst->nJobs ; j++) {
                const Job *job2 = inst->jobs + j;
                ALLOCATE_VECTOR_INI(inst->maxDurationPathInterMode[i][m][j], int, job2->nModes);
            }
        }
        qsort( job->idxModesSort, job->nModes, sizeof(IntPair), cmp_int_pair_b);
    }




    if ( nRemovedModes )
        printf("%d invalid modes were removed.\n", nRemovedModes );

    inst->projects[0].criticalPath = Project_calc_criticalPath(inst);
    /*
        FILE *file = fopen("CPD_j120.txt", "a");
        fprintf(file, "\n%d", inst->projects[0].criticalPath);
        fclose(file);
        exit(0);
    */
    Inst_computeEST( inst );

    for(int jm = 0; jm < inst->nJobs; jm++) {
        const Job* job = Inst_job(inst, jm);
        for(int mm = 0; mm < job->nModes; mm++) {
            const Mode *mode = Job_mode(job,mm);
            inst->maxDurationPath[jm][mm] = Mode_duration(mode)+somaMinDurationPathsByJob(inst, jm, job->nSucc);
            for(int jms = 0; jms < inst->nJobs; jms++) {
                const Job* job2 = Inst_job(inst, jms);
                for(int mms = 0; mms < job2->nModes; mms++) {
                    const Mode *mode2 = Job_mode(job2,mms);
                    inst->maxDurationPathInterMode[jm][mm][jms][mms] = Mode_duration(mode)+somaMinDurationPathsByJobAndInter(inst, jm, job->nSucc, jms, Mode_duration(mode2));
                }
            }
        }
    }

    ALLOCATE_VECTOR(inst->matMaxD, int*,inst->nJobs);
    for(int i  = 0 ; i < inst->nJobs ; i ++)
        ALLOCATE_VECTOR_INI(inst->matMaxD[i],int,inst->nJobs);

    floydWarshallMax(inst, inst->matMaxD, inst->nJobs);

    ALLOCATE_VECTOR(inst->matMaxDJM, int**,inst->nJobs);
    for(int i  = 0 ; i < inst->nJobs ; i ++) {
        const Job *job = Inst_job(inst, i);
        int nMode = Job_nModes(job);
        ALLOCATE_VECTOR(inst->matMaxDJM[i],int*,nMode);
        for(int m  = 0 ; m < nMode ; m ++)
            ALLOCATE_VECTOR_INI(inst->matMaxDJM[i][m],int,inst->nJobs);
    }

    maxDistanceByModes(inst, inst->matMaxD, inst->matMaxDJM);

    VecInt **paths;
    ALLOCATE_VECTOR(paths, VecInt*, Inst_nJobs(inst));
    ALLOCATE_VECTOR( inst->jobmodeinfeasible, int*, inst->nJobs);
    for(int i = 0 ; i < Inst_nJobs(inst) ; i++) {
        paths[i] = VInt_create();
        ALLOCATE_VECTOR_INI(inst->jobmodeinfeasible[i], int, inst->nMaxModes);
    }

    inst->paths = paths;

//    ALLOCATE_VECTOR(inst->unitsResourceES,int****,inst->nJobs);
    Inst_computeCompPaths( inst, inst->paths);

    //    Inst_jobsInfeasible(inst);

    return inst;
}

/* loads instance from a file MRCPSP */
Instance *Inst_create_mrcpsp( const char *dir, const char *instance )
{
    char fname[FILE_NAME_SIZE];
    sprintf( fname, "%s/%s", dir, instance );
    int nRemovedModes = 0;

    Instance *inst;
    ALLOCATE_INI( inst, Instance );

    inst->nProjects = 1;

    ALLOCATE_VECTOR( inst->projects, Project, inst->nProjects );
    ALLOCATE_VECTOR_INI( inst->idxIniResNProj, int, inst->nProjects );

    int totalJobs = 0;
    /* current index for renewableResources */

    VecInt *tcapR = VInt_create();
    VecInt *tcapNR = VInt_create();

    inst->projects[0].index = 0;

    /* processing subprojects */
    Inst_readProject_mrcpsp( dir, instance, 0, inst, tcapR, tcapNR );
    totalJobs = inst->projects[0].nJobs;

    inst->nResR = VInt_size( tcapR );
    inst->nResN = VInt_size( tcapNR );

    ALLOCATE_VECTOR_INI( inst->capR, int, VInt_size( tcapR ) );
    memcpy( inst->capR, VInt_getPtr(tcapR), sizeof(int)*VInt_size( tcapR ) );

    ALLOCATE_VECTOR_INI( inst->capN, int, VInt_size( tcapNR) );
    memcpy( inst->capN, VInt_getPtr(tcapNR), sizeof(int)*VInt_size( tcapNR ) );

    VInt_free( &tcapR );
    VInt_free( &tcapNR );

    /* adding global jobs */
    inst->nJobs = totalJobs;
    ALLOCATE_VECTOR_INI( inst->jobs, Job, totalJobs );

    int jIdx = 0;
    Job *job = inst->jobs;

    /* temporarily used as counters */
    /* adding project jobs */

    const int firstJobP = jIdx;
    Project *p = &inst->projects[0];
    p->idxFirstJob = jIdx;

    for ( int j=0; j<p->nJobs; ++j,jIdx++,job++ ) {
        job->index = jIdx;
        const Job *projJob = &p->jobs[j];
        Job_cpy( inst, job, projJob, &nRemovedModes );

        job->index = jIdx;
        job->idxOnProject = projJob->index;
        job->idxProject = p->index;
        job->origJob = projJob;

#ifdef DEBUG
        /* projects must have dummy jobs too */
        if ( j == 0 || j==p->nJobs-1  ) {
            assert( projJob->nModes==1 );
            assert( projJob->modes[0].duration == 0 );
        }
#endif

        ALLOCATE_VECTOR( job->hasPred, int, inst->nJobs );
        ALLOCATE_VECTOR( job->hasSucc, int, inst->nJobs );
        ALLOCATE_VECTOR_INI( job->hasIndPred, int, inst->nJobs );
        ALLOCATE_VECTOR_INI( job->hasIndSucc, int, inst->nJobs );

        /* updating indexes */
        for ( int k=0; k<job->nSucc; ++k ) {
            job->succ[k] += firstJobP;
            job->hasSucc[job->succ[k]] = 1; /* 1 to indicate if a job is successor
                                                Janniele*/
        }

        for ( int k=0; k<job->nPred; ++k ) {
            job->pred[k] += firstJobP;
            job->hasPred[job->pred[k]] = 1; /* 1 to indicate if a job is pred
                                Janniele*/
        }

    } /* all jobs */

    /* updating indirect preds */

    Stack *stackJobs = Stk_create(inst->nJobs);

    for(int jj = 0; jj < inst->nJobs; jj++) {

        const Job* job = Inst_job(inst, jj);

        int idxJob = Job_index(job);
        Stk_push(stackJobs, idxJob);

        while (!Stk_isEmpty(stackJobs)) {
            idxJob = Stk_pop(stackJobs);
            //            printf("\nJob: %d\n", idxJob);
            const Job *jobInd = Inst_job(inst,idxJob);
            //          printf("Preds:\n");
            for (int pred = 0; pred < Job_nPred(jobInd); pred++) {
                if (!job->hasIndPred[Job_pred(jobInd,pred)]) {
                    job->hasIndPred[Job_pred(jobInd,pred)] = 1;
                    Stk_push(stackJobs,Job_pred(jobInd,pred));
                    //                    printf("%d: %d \n", Job_pred(jobInd,pred), job->hasIndPred[Job_pred(jobInd,pred)] );
                }
            }
        }

        /* updating indirect preds*/
        idxJob = Job_index(job);
        Stk_push(stackJobs, idxJob);

        while (!Stk_isEmpty(stackJobs)) {
            idxJob = Stk_pop(stackJobs);

            const Job *jobInd = Inst_job(inst,idxJob);
            //printf("Succs:\n");
            for (int succ = 0; succ < Job_nSucc(jobInd); succ++) {
                if (!job->hasIndSucc[Job_succ(jobInd, succ)]) {
                    job->hasIndSucc[Job_succ(jobInd, succ)] = 1;
                    Stk_push(stackJobs,Job_succ(jobInd, succ));
                    //      printf("%d: %d \n", Job_succ(jobInd,succ), job->hasIndSucc[Job_succ(jobInd,succ)] );
                }
            }
        }


    }

    Stk_free(&stackJobs);

    /* maximum and minimum duration for a job */
    ALLOCATE_VECTOR_INI(inst->maxDurationPath, int*, inst->nJobs);
    ALLOCATE_VECTOR(inst->maxDurationPathInterMode, int***, inst->nJobs);

    for ( int i=0; i<inst->nJobs; ++i ) {
        Job *job = inst->jobs + i;
        job->minDuration = INT_MAX;
        job->maxDuration = 0;
        ALLOCATE_VECTOR(job->idxModesSort, IntPair, job->nModes);
        ALLOCATE_VECTOR_INI(inst->maxDurationPath[i], int, job->nModes);
        ALLOCATE_VECTOR(inst->maxDurationPathInterMode[i], int**, job->nModes);
        for ( int m=0; m<job->nModes; ++m ) {
            ALLOCATE_VECTOR(inst->maxDurationPathInterMode[i][m], int*, inst->nJobs);
            job->idxModesSort[m].a = m;
            job->idxModesSort[m].b = job->modes[m].duration;
            job->minDuration = MIN( job->minDuration, job->modes[m].duration );
            job->maxDuration = MAX( job->maxDuration, job->modes[m].duration );
            for(int j = 0 ; j < inst->nJobs ; j++) {
                const Job *job2 = inst->jobs + j;
                ALLOCATE_VECTOR_INI(inst->maxDurationPathInterMode[i][m][j], int, job2->nModes);
            }
        }
        qsort( job->idxModesSort, job->nModes, sizeof(IntPair), cmp_int_pair_b);
    }



    if ( nRemovedModes )
        printf("%d invalid moves were removed.\n", nRemovedModes );

    inst->projects[0].criticalPath = Project_calc_criticalPath(inst);

    Inst_computeEST( inst );

    for(int jm = 0; jm < inst->nJobs; jm++) {
        const Job* job = Inst_job(inst, jm);
        for(int mm = 0; mm < job->nModes; mm++) {
            const Mode *mode = Job_mode(job,mm);
            inst->maxDurationPath[jm][mm] = Mode_duration(mode)+somaMinDurationPathsByJob(inst, jm, job->nSucc);
            for(int jms = 0; jms < inst->nJobs; jms++) {
                const Job* job2 = Inst_job(inst, jms);
                for(int mms = 0; mms < job2->nModes; mms++) {
                    const Mode *mode2 = Job_mode(job2,mms);
                    inst->maxDurationPathInterMode[jm][mm][jms][mms] = Mode_duration(mode)+somaMinDurationPathsByJobAndInter(inst, jm, job->nSucc, jms, Mode_duration(mode2));
                }
            }
        }
    }

    ALLOCATE_VECTOR(inst->matMaxD, int*,inst->nJobs);
    for(int i  = 0 ; i < inst->nJobs ; i ++)
        ALLOCATE_VECTOR_INI(inst->matMaxD[i],int,inst->nJobs);

    floydWarshallMax(inst, inst->matMaxD, inst->nJobs);

    ALLOCATE_VECTOR(inst->matMaxDJM, int**,inst->nJobs);
    for(int i  = 0 ; i < inst->nJobs ; i ++) {
        const Job *job = Inst_job(inst, i);
        int nMode = Job_nModes(job);
        ALLOCATE_VECTOR(inst->matMaxDJM[i],int*,nMode);
        for(int m  = 0 ; m < nMode ; m ++)
            ALLOCATE_VECTOR_INI(inst->matMaxDJM[i][m],int,inst->nJobs);
    }

    maxDistanceByModes(inst, inst->matMaxD, inst->matMaxDJM);
//    ALLOCATE_VECTOR(inst->unitsResourceES,int****,inst->nJobs);
    VecInt **paths;
    ALLOCATE_VECTOR(paths, VecInt*, Inst_nJobs(inst) );
    for(int i = 0 ; i < Inst_nJobs(inst); i++)
        paths[i] = VInt_create();

    inst->paths = paths;

    Inst_computeCompPaths( inst, inst->paths);

    /* for(int i = 0 ; i < Inst_nJobs(inst); i++) {
         int ss = VInt_size(paths[i]);
         printf("%d : ", i);
         for(int j = 0 ; j < ss ; j++)
             printf(" %d ", VInt_get(paths[i], j));
         printf("\n");
     }*/


    return inst;
}

Instance *Inst_create( const char *dir, const char *instance )
{
    char fname[FILE_NAME_SIZE];
    sprintf( fname, "%s/%s", dir, instance );
    int nRemovedModes = 0;

    VecStr *linesInst = VStr_create( STR_SIZE );
    VStr_readFrom( linesInst, fname, True );

    Instance *inst;
    ALLOCATE_INI( inst, Instance );
    int r;

    r = sscanf( VStr_get( linesInst, 0 ), "%d", &(inst->nProjects) );
    assert( r==1 );

    ALLOCATE_VECTOR( inst->projects, Project, inst->nProjects );
    ALLOCATE_VECTOR_INI( inst->idxIniResNProj, int, inst->nProjects );

    VecStr *projectFiles = VStr_create( FILE_NAME_SIZE );

    int line = 1;
    for ( int i=0; (i<inst->nProjects); ++i ) {
        Project *p = &(inst->projects[ i ]);
        p->index = i;
        r = sscanf( VStr_get( linesInst, line++ ), "%d", &(p->releaseDate) );
        assert( r==1 );
        r = sscanf( VStr_get( linesInst, line++ ), "%d", &(p->criticalPath) );
        assert( r==1 );
        // subproject
        char fnamesp[FILE_NAME_SIZE];
        r = sscanf( VStr_get( linesInst, line++ ), "%s", fnamesp );
        assert( r==1 );

        VStr_pushBack( projectFiles, fnamesp );

        //printf("> %d %d %d %s\n", p->index, p->releaseDate, p->criticalPath, fnamesp );
    }

    int nres;
    r = sscanf( VStr_get( linesInst, line++ ), "%d", &nres );
    assert( r==1 );

    //printf("\nT: %d\n", nres); getchar();

    int *globalRCap;
    ALLOCATE_VECTOR( globalRCap, int, nres );

    Tokenizer *tok = Tok_create();

    char cleanLine[LINE_SIZE];
    strcpy( cleanLine, VStr_get( linesInst, line++ ) );
    strRemoveDblSpaces( cleanLine );

    Tok_parse( tok, cleanLine, ' ' );

    {
        int k = 0;
        for ( int i=0; i<Tok_nTokens(tok); ++i ) {
            char cell[STR_SIZE];
            strcpy( cell, Tok_token( tok, i ) );
            strRemoveSpsEol( cell );
            if ( !strlen( cell ) )
                continue;

            globalRCap[k++] = atoi( cell );
        }

        int nGlobalResR = 0;
        for ( int i=0; (i<k); ++i )
            if ( globalRCap[i] != -1 )
                nGlobalResR++;

        assert( k==nres );
        inst->nResRGlobal = nGlobalResR;
    }


    Tok_free( &tok );
    VStr_free( &linesInst );


    int totalJobs = r-r; //0
    /* current index for renewableResources */

    VecInt *tcapR = VInt_create();
    VecInt *tcapNR = VInt_create();

    /* global resources will be the first renewable resources added */
    for ( int i=0; (i<nres); ++i )
        if ( globalRCap[i] != -1 )
            VInt_pushBack( tcapR, globalRCap[i] );

    /* processing subprojects */
    for ( int i=0; (i<VStr_size(projectFiles)); ++i ) {
        Inst_readProject( dir, VStr_get( projectFiles, i), i, inst, nres, globalRCap, tcapR, tcapNR );
        totalJobs += inst->projects[i].nJobs;
    }
    inst->nResR = VInt_size( tcapR );
    inst->nResN = VInt_size( tcapNR );

    ALLOCATE_VECTOR_INI( inst->capR, int, VInt_size( tcapR ) );
    memcpy( inst->capR, VInt_getPtr(tcapR), sizeof(int)*VInt_size( tcapR ) );

    ALLOCATE_VECTOR_INI( inst->capN, int, VInt_size( tcapNR) );
    memcpy( inst->capN, VInt_getPtr(tcapNR), sizeof(int)*VInt_size( tcapNR ) );

    VInt_free( &tcapR );
    VInt_free( &tcapNR );

    VStr_free( &projectFiles );
    free( globalRCap );

    /* adding global jobs */
    inst->nJobs = totalJobs;
    ALLOCATE_VECTOR_INI( inst->jobs, Job, totalJobs );

    int jIdx = 0;
    Job *job = inst->jobs;

    /* temporarily used as counters */
    /* adding project jobs */
    for ( int i=0; (i<inst->nProjects); ++i ) {
        const int firstJobP = jIdx;
        Project *p = &inst->projects[i];
        p->idxFirstJob = jIdx;

        for ( int j=0; j<p->nJobs; ++j,jIdx++,job++ ) {
            job->index = jIdx;
            const Job *projJob = &p->jobs[j];
            Job_cpy( inst, job, projJob, &nRemovedModes );

            job->index = jIdx;
            job->idxOnProject = projJob->index;
            job->idxProject = p->index;
            job->origJob = projJob;

#ifdef DEBUG
            /* projects must have dummy jobs too */
            if ( j == 0 || j==p->nJobs-1  ) {
                assert( projJob->nModes==1 );
                assert( projJob->modes[0].duration == 0 );
            }
#endif

            ALLOCATE_VECTOR( job->hasPred, int, inst->nJobs );
            ALLOCATE_VECTOR( job->hasSucc, int, inst->nJobs );
            ALLOCATE_VECTOR_INI( job->hasIndPred, int, inst->nJobs );
            ALLOCATE_VECTOR_INI( job->hasIndSucc, int, inst->nJobs );


            /* updating indexes */
            for ( int k=0; k<job->nSucc; ++k ) {
                job->succ[k] += firstJobP;
                job->hasSucc[job->succ[k]] = 1; /* 1 to indicate if a job is successor
                                                    Janniele*/
            }

            for ( int k=0; k<job->nPred; ++k ) {
                job->pred[k] += firstJobP;
                job->hasPred[job->pred[k]] = 1; /* 1 to indicate if a job is pred
                                    Janniele*/
            }

        } /* all jobs */


    } /* all projects */


    /* updating indirect preds */

    Stack *stackJobs = Stk_create(inst->nJobs);

    for(int jj = 0; jj < inst->nJobs; jj++) {

        const Job* job = Inst_job(inst, jj);

        int idxJob = Job_index(job);
        Stk_push(stackJobs, idxJob);

        while (!Stk_isEmpty(stackJobs)) {
            idxJob = Stk_pop(stackJobs);
            //            printf("\nJob: %d\n", idxJob);
            const Job *jobInd = Inst_job(inst,idxJob);
            //          printf("Preds:\n");
            for (int pred = 0; pred < Job_nPred(jobInd); pred++) {
                if (!job->hasIndPred[Job_pred(jobInd,pred)]) {
                    job->hasIndPred[Job_pred(jobInd,pred)] = 1;
                    Stk_push(stackJobs,Job_pred(jobInd,pred));
                    //                    printf("%d: %d \n", Job_pred(jobInd,pred), job->hasIndPred[Job_pred(jobInd,pred)] );
                }
            }
        }

        /* updating indirect preds*/
        idxJob = Job_index(job);
        Stk_push(stackJobs, idxJob);

        while (!Stk_isEmpty(stackJobs)) {
            idxJob = Stk_pop(stackJobs);
            //printf("\nJob: %d\n", idxJob);
            const Job *jobInd = Inst_job(inst,idxJob);
            //printf("Succs:\n");
            for (int succ = 0; succ < Job_nSucc(jobInd); succ++) {
                if (!job->hasIndSucc[Job_succ(jobInd, succ)]) {
                    job->hasIndSucc[Job_succ(jobInd, succ)] = 1;
                    Stk_push(stackJobs,Job_succ(jobInd, succ));
                    //      printf("%d: %d \n", Job_succ(jobInd,succ), job->hasIndSucc[Job_succ(jobInd,succ)] );
                }
            }
        }


    }

    Stk_free(&stackJobs);

    /* maximum and minimum duration for a job */
    ALLOCATE_VECTOR_INI(inst->maxDurationPath, int*, inst->nJobs);
    ALLOCATE_VECTOR(inst->maxDurationPathInterMode, int***, inst->nJobs);
    for ( int i=0; i<inst->nJobs; ++i ) {
        Job *job = inst->jobs + i;
        job->minDuration = INT_MAX;
        job->maxDuration = 0;
        ALLOCATE_VECTOR(job->idxModesSort, IntPair, job->nModes);
        ALLOCATE_VECTOR_INI(inst->maxDurationPath[i], int, job->nModes);
        ALLOCATE_VECTOR(inst->maxDurationPathInterMode[i], int**, job->nModes);
        for ( int m=0; m<job->nModes; ++m ) {
            ALLOCATE_VECTOR(inst->maxDurationPathInterMode[i][m], int*, inst->nJobs);
            job->idxModesSort[m].a = m;
            job->idxModesSort[m].b = job->modes[m].duration;
            job->minDuration = MIN( job->minDuration, job->modes[m].duration );
            job->maxDuration = MAX( job->maxDuration, job->modes[m].duration );
            for(int j = 0 ; j < inst->nJobs ; j++) {
                const Job *job2 = inst->jobs + j;
                ALLOCATE_VECTOR_INI(inst->maxDurationPathInterMode[i][m][j], int, job2->nModes);
            }
        }
        qsort( job->idxModesSort, job->nModes, sizeof(IntPair), cmp_int_pair_b);
    }



    if ( nRemovedModes )
        printf("%d invalid moves were removed.\n", nRemovedModes );

    Inst_computeEST( inst );

    for(int jm = 0; jm < inst->nJobs; jm++) {
        const Job* job = Inst_job(inst, jm);
        for(int mm = 0; mm < job->nModes; mm++) {
            const Mode *mode = Job_mode(job,mm);
            inst->maxDurationPath[jm][mm] = Mode_duration(mode)+somaMinDurationPathsByJob(inst, jm, job->nSucc);
            for(int jms = 0; jms < inst->nJobs; jms++) {
                const Job* job2 = Inst_job(inst, jms);
                for(int mms = 0; mms < job2->nModes; mms++) {
                    const Mode *mode2 = Job_mode(job2,mms);
                    inst->maxDurationPathInterMode[jm][mm][jms][mms] = Mode_duration(mode)+somaMinDurationPathsByJobAndInter(inst, jm, job->nSucc, jms, Mode_duration(mode2));
                }
            }
        }
    }

    ALLOCATE_VECTOR(inst->matMaxD, int*,inst->nJobs);
    for(int i  = 0 ; i < inst->nJobs ; i ++)
        ALLOCATE_VECTOR_INI(inst->matMaxD[i],int,inst->nJobs);

    floydWarshallMax(inst, inst->matMaxD, inst->nJobs);

    ALLOCATE_VECTOR(inst->matMaxDJM, int**,inst->nJobs);
    for(int i  = 0 ; i < inst->nJobs ; i ++) {
        const Job *job = Inst_job(inst, i);
        int nMode = Job_nModes(job);
        ALLOCATE_VECTOR(inst->matMaxDJM[i],int*,nMode);
        for(int m  = 0 ; m < nMode ; m ++)
            ALLOCATE_VECTOR_INI(inst->matMaxDJM[i][m],int,inst->nJobs);
    }

    maxDistanceByModes(inst, inst->matMaxD, inst->matMaxDJM);
//    ALLOCATE_VECTOR(inst->unitsResourceES,int****,inst->nJobs);
    VecInt **paths;
    ALLOCATE_VECTOR(paths, VecInt*, Inst_nJobs(inst));
    for(int i = 0 ; i < Inst_nJobs(inst) ; i++)
        paths[i] = VInt_create();

    inst->paths = paths;

    Inst_computeCompPaths( inst, inst->paths);



    return inst;
}


/*int Job_durationModeSort( const Job *job, int m)
{
    return job->idxModesSort[m].b;
}
*/

/*int Job_idxModeSort( const Job *job, int m)
{
    return job->idxModesSort[m].a;
}*/
/* when reading instance resource usage is informed in format R R N N,
 * where fist the renewable resources are used and latter the non renewable ones */
static void Inst_fillModeResUsage( const Instance *inst, Mode *m,
                                   int nresR, int nresN, const int *usage,
                                   int idxStartR, int idxStartNR, const int globalRCap[] )
{


    /* usage of renewable resources for this mode */
    for ( int i=0; (i<nresR); ++i ) // counting
        m->nResR += (int) (usage[i]>0);

    /* adding non zero entries of renewable resources */
    int nAddedGlobal = 0, nAddedLocal = 0;
    if ( m->nResR ) {
        ALLOCATE_VECTOR( m->idxResR, int, m->nResR );
        ALLOCATE_VECTOR( m->useResR, int, m->nResR );
        m->nResR = 0;
        for ( int i=0; (i<nresR); ++i )
            if (usage[i] && globalRCap[i] == -1) {
                m->useResR[ m->nResR ] = usage[i];
                m->idxResR[ m->nResR ] = idxStartR + nAddedLocal;   // index refers to local resource
                ++nAddedLocal;
                m->nResR++;
            } else {
                if(usage[i]) {
                    m->useResR[ m->nResR ] = usage[i];
                    m->idxResR[ m->nResR ] = nAddedGlobal;              // index refers to global resource
                    m->nResR++;
                }
                if(globalRCap[i] >0)
                    ++nAddedGlobal;
            }

    } else {
        m->idxResR = NULL;
        m->useResR = NULL;
    }

    /* local non-renewable resources */
    m->nResN = 0;
    for ( int k=0; (k<nresN); ++k )
        m->nResN += (int) (usage[nresR+k]>0);

    if ( m->nResN ) {
        ALLOCATE_VECTOR( m->idxResN, int, m->nResN );
        ALLOCATE_VECTOR( m->useResN, int, m->nResN );
        m->nResN = 0;
        for ( int k=0; (k<nresN); ++k )
            if ( usage[nresR+k]>0 ) {
                m->idxResN[m->nResN] = idxStartNR + k;
                m->useResN[m->nResN] = usage[nresR+k];
                m->nResN++;
            }
    } else {
        m->idxResN = NULL;
        m->useResN = NULL;
    }
}

Instance *Inst_read( char ** argv, int argc )
{
    Instance *inst;

    for(int n = 0; n< argc ; n++) {
        if (strcmp(argv[n],"-rcpsp") == 0) {
            inst = Inst_create_rcpsp( argv[1], argv[2] );
            return inst;
        }
        if (strcmp(argv[n],"-mrcpsp") == 0) {
            inst = Inst_create_mrcpsp( argv[1], argv[2] );
            return inst;
        }
        if (strcmp(argv[n],"-mmrcmpsp") == 0) {
            inst = Inst_create( argv[1], argv[2] );
            return inst;
        }
    }

    printf("\nTipo de problema nao especificado. [ -rcpsp, -mrcpsp ou -mmrcmpsp ]\n");
    exit(1);

}

static void Inst_readProject( const char *dir, const char *projectFile, int projIdx,
                              Instance *inst, int nres, const int globalRCap[],
                              VecInt *capR, VecInt *capNR )
{
    char fname[FILE_NAME_SIZE];
    sprintf( fname, "%s/%s", dir, projectFile );
    char line[LINE_SIZE];
    char *s;

    FILE *f = fopen( fname, "r" );
    assert( f );

    Project *p = &(inst->projects[projIdx]);


    /* contents in the format   field(substr)   :  value  */
#define N_FIELDS 5
    const char field[N_FIELDS][STR_SIZE] = { "projects", "supersource", " renewable", "nonrenewable", "doubly" };
    int nprojs = 0, jobs = 0, nrenew = 0, nnonrenew = 0, ndoubly = 0;
    int *fieldPtr[ N_FIELDS ];
    fieldPtr[0] = &nprojs;
    fieldPtr[1] = &jobs;
    fieldPtr[2] = &nrenew;
    fieldPtr[3] = &nnonrenew;
    fieldPtr[4] = &ndoubly;
    int currField = 0;
    p->nJobs = 0;

    Tokenizer *tok = Tok_create();

    // int totalModes = 0;

    while ( (s=fgets( line, LINE_SIZE, f)) ) {
        /* first fields */
        if (strstr( s, "rel.date" )) {
            int head[6];    /* pronr.  #jobs rel.date duedate tardcost  MPM-Time */

            for(int i=0; i<6; ++i)
                fscanf(f, "%d", &head[i]);

            //printf("\nRD: %d", head[2]); getchar();
            //p->releaseDate = head[2];
            p->dueDate = head[3];
            p->tardCost = head[4];
            p->mpmTime = head[5];
        }

        if ( currField<N_FIELDS ) {
            if (strstr( line, field[currField] )) {
                char *sp = strstr( s, ":" );
                assert( sp );
                ++sp;
                assert( *sp != '\0' );
                int r = sscanf( sp, "%d", fieldPtr[currField] );
                assert( r==1 );
                //currField++;
                currField= currField +1 + (r-r);
                continue;
            }
        }

        if (strstr( s, "#successors" )) {
            assert( p->nJobs==0 );
            assert( jobs >= 2 );
            p->nJobs = jobs;
            ALLOCATE_VECTOR_INI( p->jobs, Job, p->nJobs  );

            int succs[jobs+3];

            for ( int i=0; i<jobs; ++i ) {
                s=fgets( line, LINE_SIZE, f);
                assert( s );
                strRemoveDblSpaces( s );
                readIntVectorStr( s, tok, 3, jobs+3, succs );

#ifdef DEBUG
                const int jIdx   = succs[0];
                assert( jIdx==i+1 );
#endif

                const int jModes = succs[1];
                const int jnSucc = succs[2];


                assert( jModes>=1 );
                assert( jnSucc>=0 && jnSucc<=jobs-1 );

                Job *job = &(p->jobs[i]);
                job->index = i;
                job->idxOnProject = i;
                job->nModes = jModes;
                job->nSucc  = jnSucc;
                job->succ   = NULL;
                if ( job->nSucc )
                    ALLOCATE_VECTOR( job->succ, int, job->nSucc );

                ALLOCATE_VECTOR_INI( job->modes, Mode, job->nModes );
                for ( int j=0; (j<job->nModes); ++j )
                    job->modes[j].index = j;

                // totalModes += job->nModes;
                if (job->nSucc)
                    memcpy( job->succ, &succs[3], sizeof(int)*job->nSucc );
                for ( int j=0; (j<job->nSucc); ++j ) {
                    job->succ[j]--;
                    assert( job->succ[j]>=1 && job->succ[j]<jobs );
                    p->jobs[job->succ[j]].nPred++;
                }
            }

            /* allocating pred */
            for ( int i=0; (i<p->nJobs); ++i ) {
                p->jobs[i].pred = NULL;
                if ( p->jobs[i].nPred )
                    ALLOCATE_VECTOR( p->jobs[i].pred, int, p->jobs[i].nPred );
                p->jobs[i].nPred = 0;
            }

            /* filling preds */
            for ( int i=0; (i<p->nJobs); ++i ) {
                for ( int j=0; (j<p->jobs[i].nSucc); ++j ) {
                    Job *jobSucc = &(p->jobs[p->jobs[i].succ[j]]);
                    jobSucc->pred[jobSucc->nPred]= i;
                    jobSucc->nPred++;
                }
            }

            continue;
        }

        if (strstr( s, "mode duration" )) {
            assert( countChar( s, 'R') == nrenew );
            assert( countChar( s, 'N') == nnonrenew );
            assert( firstOccurrence( s, 'N') > lastOccurrence( s, 'R' ) );

            int durLine[ nres+3 ];

            const int idxStartRLocalR  = VInt_size( capR );
            const int idxStartRLocalNR = VInt_size( capNR );
            inst->idxIniResNProj[p->index] = idxStartRLocalNR;

            s=fgets( line, LINE_SIZE, f);
            assert( s );
            assert( strstr(s, "--") );

            for ( int i=0; i<p->nJobs; ++i ) {
                Job *job = &p->jobs[i];
                job->minDuration = INT_MAX;
                job->maxDuration = 0;
                s=fgets( line, LINE_SIZE, f);
                assert( s );

                readIntVectorStr( s, tok, nres+3, nres+3, durLine );
                job->modes[0].duration = durLine[2];
                assert( nrenew>0 && nnonrenew>0 );
                Inst_fillModeResUsage( inst, job->modes, nrenew, nnonrenew, durLine+3, idxStartRLocalR, idxStartRLocalNR, globalRCap );

                for ( int j=1; j<job->nModes; ++j ) {
                    s=fgets( line, LINE_SIZE, f);
                    assert( s );
                    readIntVectorStr( s, tok, nres+2, nres+2, durLine );
                    job->modes[j].duration = durLine[1];
                    job->minDuration = MIN( job->minDuration, job->modes[j].duration );
                    job->maxDuration = MAX( job->maxDuration, job->modes[j].duration );

                    Mode *m = job->modes + j;

                    Inst_fillModeResUsage( inst, m, nrenew, nnonrenew, durLine+2, idxStartRLocalR, idxStartRLocalNR, globalRCap );
                }

                inst->nMaxModes = MAX( inst->nMaxModes, job->nModes );
            }


            continue;
        }

        if (strstr( s, "RESOURCEAVAILABILITIES" )) {
            s=fgets( line, LINE_SIZE, f);
            assert( s );
            assert( nres == nrenew+nnonrenew );

#ifdef DEBUG
            /* checking if it is in the format R R (renewable)  followed by non renewable */
            char *pt = NULL;
            char *sp = s;
            for ( int i=0; i<nrenew; ++i ) {
                assert( (pt=strstr(sp,"R")) );
                sp = pt+1;
                pt = NULL;
            }
#endif

            s=fgets( line, LINE_SIZE, f);
            assert( s );

            //int resAV[ nrenew ];
            int *resAV = (int*)malloc( nres * sizeof(int));

            readIntVectorStr( s, tok, nres, nres, resAV );

            for ( int i=0; i<nrenew; ++i )
                if ( globalRCap[i] == -1 )
                    VInt_pushBack( capR, resAV[i] );

            {
                for ( int i=nrenew; i<nres; ++i)
                    VInt_pushBack( capNR, resAV[i] );
            }

            free(resAV);

            continue;
        }

    }

    Tok_free( &tok );

    /* validating data readed before */
    assert( nres == nrenew+nnonrenew );
    assert( currField == N_FIELDS );
    assert( nprojs == 1 );  /* should be only one project here */
    assert( ndoubly == 0 );

    fclose( f );
}

static void Inst_readProject_mrcpsp( const char *dir, const char *projectFile, int projIdx,
                                     Instance *inst, VecInt *capR, VecInt *capNR )
{
    char fname[FILE_NAME_SIZE];
    sprintf( fname, "%s/%s", dir, projectFile );
    char line[LINE_SIZE];
    char *s;
    int nres = 0;

    FILE *f = fopen( fname, "r" );
    assert( f );

    Project *p = &(inst->projects[projIdx]);


    /* contents in the format   field(substr)   :  value  */
#define N_FIELDS 5
    const char field[N_FIELDS][STR_SIZE] = { "projects", "supersource", " renewable", "nonrenewable", "doubly" };
    int nprojs = 0, jobs = 0, nrenew = 0, nnonrenew = 0, ndoubly = 0;
    int *fieldPtr[ N_FIELDS ];
    fieldPtr[0] = &nprojs;
    fieldPtr[1] = &jobs;
    fieldPtr[2] = &nrenew;
    fieldPtr[3] = &nnonrenew;
    fieldPtr[4] = &ndoubly;
    int currField = 0;
    p->nJobs = 0;

    Tokenizer *tok = Tok_create();

    //    int totalModes = 0;
    int *globalRCap;

    while ( (s=fgets( line, LINE_SIZE, f)) ) {

        if (strstr( s, "RESOURCES" )) {
            int Rnum, Nnum;
            char vetCharRN[100];

            for(int h=0; h<3; ++h)
                fscanf(f, "%s", (char*)&vetCharRN);

            fscanf(f, "%d", &Rnum);

            nres += Rnum;

            for(int h=0; h<4; ++h)
                fscanf(f, "%s", (char*)&vetCharRN);

            fscanf(f, "%d", &Nnum);

            nres += Nnum;
        }
    }

    fseek(f, 0, SEEK_SET);

    ALLOCATE_VECTOR( globalRCap, int, nres );

    /* first: obtain project resources (required by: mrcpsp and rcpsp) */
    while ( (s=fgets( line, LINE_SIZE, f)) ) {
        if (strstr( s, "RESOURCEAVAILABILITIES" )) {

            int locRes[nres];
            int num;
            char charRN;
            char vetCharRN[nres];

            int z=0;
            while(z < nres) {
                fscanf(f, "%c", &charRN);
                fscanf(f, "%d", &num);

                if(charRN == ' ')
                    continue;

                ++z;

                if(charRN == 'R')
                    vetCharRN[num-1] = charRN;
                else
                    vetCharRN[num-1+(nres/2)] = charRN;
            }

            for(int i=0; i<nres; ++i)
                fscanf(f, "%d", &locRes[i]);

            for(int i=0; i<nres; ++i)
                if(vetCharRN[i] == 'R')
                    globalRCap[i] = locRes[i];
                else
                    globalRCap[i] = -1;
        }
    }

    //printf("\n %d %d %d %d", globalRCap[0], globalRCap[1], globalRCap[2], globalRCap[3] ); getchar();

    fseek(f, 0, SEEK_SET);

    /* global resources will be the first renewable resources added */
    for ( int i=0; (i<nres); ++i )
        if ( globalRCap[i] != -1 )
            VInt_pushBack( capR, globalRCap[i] );


    int nGlobalResR = 0;
    for ( int i=0; (i<nres); ++i )
        if ( globalRCap[i] != -1 )
            nGlobalResR++;

    inst->nResRGlobal = nGlobalResR;

    //printf("\n %d ", nGlobalResR); getchar();

    while ( (s=fgets( line, LINE_SIZE, f)) ) {
        /* first fields */
        if (strstr( s, "rel.date" )) {
            int head[6];    /* pronr.  #jobs rel.date duedate tardcost  MPM-Time */

            for(int i=0; i<6; ++i)
                fscanf(f, "%d", &head[i]);

            p->releaseDate = head[2];
            p->dueDate = head[3];
            p->tardCost = head[4];
            p->mpmTime = head[5];
        }

        if ( currField<N_FIELDS ) {
            if (strstr( line, field[currField] )) {
                char *sp = strstr( s, ":" );
                assert( sp );
                ++sp;
                assert( *sp != '\0' );
                int r = sscanf( sp, "%d", fieldPtr[currField] );
                assert( r==1 );
                //currField++;
                currField= currField +1 + (r-r);

                continue;
            }
        }

        if (strstr( s, "#successors" )) {
            assert( p->nJobs==0 );
            assert( jobs >= 2 );
            p->nJobs = jobs;
            ALLOCATE_VECTOR_INI( p->jobs, Job, p->nJobs  );

            int succs[jobs+3];

            for ( int i=0; i<jobs; ++i ) {
                s=fgets( line, LINE_SIZE, f);
                assert( s );
                strRemoveDblSpaces( s );
                readIntVectorStr( s, tok, 3, jobs+3, succs );

#ifdef DEBUG
                const int jIdx   = succs[0];
                assert( jIdx==i+1 );
#endif

                const int jModes = succs[1];
                const int jnSucc = succs[2];

                assert( jModes>=1 );
                assert( jnSucc>=0 && jnSucc<=jobs-1 );

                Job *job = &(p->jobs[i]);
                job->index = i;
                job->idxOnProject = i;
                job->nModes = jModes;
                job->nSucc  = jnSucc;
                job->succ   = NULL;
                if ( job->nSucc )
                    ALLOCATE_VECTOR( job->succ, int, job->nSucc );

                ALLOCATE_VECTOR_INI( job->modes, Mode, job->nModes );
                for ( int j=0; (j<job->nModes); ++j )
                    job->modes[j].index = j;

                //                totalModes += job->nModes;
                if (job->nSucc)
                    memcpy( job->succ, &succs[3], sizeof(int)*job->nSucc );
                for ( int j=0; (j<job->nSucc); ++j ) {
                    job->succ[j]--;
                    assert( job->succ[j]>=1 && job->succ[j]<jobs );
                    p->jobs[job->succ[j]].nPred++;
                }
            }

            /* allocating pred */
            for ( int i=0; (i<p->nJobs); ++i ) {
                p->jobs[i].pred = NULL;
                if ( p->jobs[i].nPred )
                    ALLOCATE_VECTOR( p->jobs[i].pred, int, p->jobs[i].nPred );
                p->jobs[i].nPred = 0;
            }

            /* filling preds */
            for ( int i=0; (i<p->nJobs); ++i ) {
                for ( int j=0; (j<p->jobs[i].nSucc); ++j ) {
                    Job *jobSucc = &(p->jobs[p->jobs[i].succ[j]]);
                    jobSucc->pred[jobSucc->nPred]= i;
                    jobSucc->nPred++;
                }
            }

            continue;
        }

        if (strstr( s, "mode duration" )) {
            assert( countChar( s, 'R') == nrenew );
            assert( countChar( s, 'N') == nnonrenew );
            assert( firstOccurrence( s, 'N') > lastOccurrence( s, 'R' ) );

            int durLine[ nres+3 ];

            const int idxStartRLocalR  = VInt_size( capR );
            const int idxStartRLocalNR = VInt_size( capNR );
            inst->idxIniResNProj[p->index] = idxStartRLocalNR;

            s=fgets( line, LINE_SIZE, f);
            assert( s );
            assert( strstr(s, "--") );

            for ( int i=0; i<p->nJobs; ++i ) {
                Job *job = &p->jobs[i];
                job->minDuration = INT_MAX;
                job->maxDuration = 0;
                s=fgets( line, LINE_SIZE, f);
                assert( s );

                readIntVectorStr( s, tok, nres+3, nres+3, durLine );
                job->modes[0].duration = durLine[2];
                assert( nrenew>0 && nnonrenew>=0 );

                Inst_fillModeResUsage( inst, job->modes, nrenew, nnonrenew, durLine+3, idxStartRLocalR, idxStartRLocalNR, globalRCap );

                for ( int j=1; j<job->nModes; ++j ) {
                    s=fgets( line, LINE_SIZE, f);
                    assert( s );
                    readIntVectorStr( s, tok, nres+2, nres+2, durLine );
                    job->modes[j].duration = durLine[1];
                    job->minDuration = MIN( job->minDuration, job->modes[j].duration );
                    job->maxDuration = MAX( job->maxDuration, job->modes[j].duration );

                    Mode *m = job->modes + j;

                    Inst_fillModeResUsage( inst, m, nrenew, nnonrenew, durLine+2, idxStartRLocalR, idxStartRLocalNR, globalRCap );
                }

                inst->nMaxModes = MAX( inst->nMaxModes, job->nModes );
            }

            continue;
        }

        if (strstr( s, "RESOURCEAVAILABILITIES" )) {
            s=fgets( line, LINE_SIZE, f);
            assert( s );
            assert( nres == nrenew+nnonrenew );

#ifdef DEBUG
            /* checking if it is in the format R R (renewable)  followed by non renewable */
            char *pt = NULL;
            char *sp = s;
            for ( int i=0; i<nrenew; ++i ) {
                assert( (pt=strstr(sp,"R")) );
                sp = pt+1;
                pt = NULL;
            }
#endif

            s=fgets( line, LINE_SIZE, f);
            assert( s );

            //int resAV[ nrenew ];
            int *resAV = (int*)malloc( nres * sizeof(int));

            readIntVectorStr( s, tok, nres, nres, resAV );


            for ( int i=0; i<nrenew; ++i )
                if ( globalRCap[i] == -1 )
                    VInt_pushBack( capR, resAV[i] );

            {
                for ( int i=nrenew; i<nres; ++i)
                    VInt_pushBack( capNR, resAV[i] );
            }

            free(resAV);

            continue;
        }

    }

    Tok_free( &tok );

    free( globalRCap );

    /* validating data readed before */
    assert( nres == nrenew+nnonrenew );
    assert( currField == N_FIELDS );
    assert( nprojs == 1 );  /* should be only one project here */
    assert( ndoubly == 0 );

    fclose( f );
}


static void Inst_readProject_rcpsp( const char *dir, const char *projectFile, int projIdx,
                                    Instance *inst, VecInt *capR, VecInt *capNR )
{
    char fname[FILE_NAME_SIZE];
    sprintf( fname, "%s/%s", dir, projectFile );
    char line[LINE_SIZE];
    char *s;
    int nres = 0;

    FILE *f = fopen( fname, "r" );
    assert( f );

    Project *p = &(inst->projects[projIdx]);

    /* contents in the format   field(substr)   :  value  */
#define N_FIELDS 5
    const char field[N_FIELDS][STR_SIZE] = { "projects", "supersource", " renewable", "nonrenewable", "doubly" };
    int nprojs = 0, jobs = 0, nrenew = 0, nnonrenew = 0, ndoubly = 0;
    int *fieldPtr[ N_FIELDS ];
    fieldPtr[0] = &nprojs;
    fieldPtr[1] = &jobs;
    fieldPtr[2] = &nrenew;
    fieldPtr[3] = &nnonrenew;
    fieldPtr[4] = &ndoubly;
    int currField = 0;
    p->nJobs = 0;

    Tokenizer *tok = Tok_create();

    //    int totalModes = 0;
    int *globalRCap;

    while ( (s=fgets( line, LINE_SIZE, f)) ) {

        if (strstr( s, "RESOURCES" )) {
            int Rnum, Nnum;
            char vetCharRN[100];

            for(int h=0; h<3; ++h)
                fscanf(f, "%s", (char*)&vetCharRN);

            fscanf(f, "%d", &Rnum);

            nres += Rnum;

            for(int h=0; h<4; ++h)
                fscanf(f, "%s", (char*)&vetCharRN);

            fscanf(f, "%d", &Nnum);

            nres += Nnum;
        }
    }

    fseek(f, 0, SEEK_SET);

    ALLOCATE_VECTOR( globalRCap, int, nres );

    /* first: obtain project resources (required by: mrcpsp and rcpsp) */
    while ( (s=fgets( line, LINE_SIZE, f)) ) {
        if (strstr( s, "RESOURCEAVAILABILITIES" )) {
            //s=fgets( line, LINE_SIZE, f);
            //s=fgets( line, LINE_SIZE, f);

            int locRes[nres];
            int num;
            char charRN;
            char vetCharRN[nres];

            int z=0;
            while(z < nres) {
                fscanf(f, "%c", &charRN);
                fscanf(f, "%d", &num);

                if(charRN == ' ')
                    continue;

                ++z;

                if(charRN == 'R')
                    vetCharRN[num-1] = charRN;
                else
                    vetCharRN[num-1+(nres/2)] = charRN;
            }

            for(int i=0; i<nres; ++i)
                fscanf(f, "%d", &locRes[i]);

            for(int i=0; i<nres; ++i)
                if(vetCharRN[i] == 'R')
                    globalRCap[i] = locRes[i];
                else
                    globalRCap[i] = -1;
        }
    }

    fseek(f, 0, SEEK_SET);

    /* global resources will be the first renewable resources added */
    for ( int i=0; (i<nres); ++i )
        if ( globalRCap[i] != -1 )
            VInt_pushBack( capR, globalRCap[i] );

    int nGlobalResR = 0;
    for ( int i=0; (i<nres); ++i )
        if ( globalRCap[i] != -1 )
            nGlobalResR++;

    inst->nResRGlobal = nGlobalResR;

    while ( (s=fgets( line, LINE_SIZE, f)) ) {
        /* first fields */
        if (strstr( s, "rel.date" )) {
            int head[6];    /* pronr.  #jobs rel.date duedate tardcost  MPM-Time */

            for(int i=0; i<6; ++i)
                fscanf(f, "%d", &head[i]);

            p->releaseDate = head[2];
            p->dueDate = head[3];
            p->tardCost = head[4];
            p->mpmTime = head[5];
        }

        if ( currField<N_FIELDS ) {
            if (strstr( line, field[currField] )) {
                char *sp = strstr( s, ":" );
                assert( sp );
                ++sp;
                assert( *sp != '\0' );
                int r = sscanf( sp, "%d", fieldPtr[currField] );
                assert( r==1 );
                //currField++;
                currField= currField +1 + (r-r);
                continue;
            }
        }

        if (strstr( s, "#successors" )) {
            assert( p->nJobs==0 );
            assert( jobs >= 2 );
            p->nJobs = jobs;
            ALLOCATE_VECTOR_INI( p->jobs, Job, p->nJobs  );

            int succs[jobs+3];

            for ( int i=0; i<jobs; ++i ) {
                s=fgets( line, LINE_SIZE, f);
                assert( s );
                strRemoveDblSpaces( s );
                readIntVectorStr( s, tok, 3, jobs+3, succs );

#ifdef DEBUG
                const int jIdx   = succs[0];
                assert( jIdx==i+1 );
#endif

                const int jModes = succs[1];
                const int jnSucc = succs[2];

                assert( jModes>=1 );
                assert( jnSucc>=0 && jnSucc<=jobs-1 );

                Job *job = &(p->jobs[i]);
                job->index = i;
                job->idxOnProject = i;
                job->nModes = jModes;
                job->nSucc  = jnSucc;
                job->succ   = NULL;
                if ( job->nSucc )
                    ALLOCATE_VECTOR( job->succ, int, job->nSucc );

                ALLOCATE_VECTOR_INI( job->modes, Mode, job->nModes );
                for ( int j=0; (j<job->nModes); ++j )
                    job->modes[j].index = j;

                //                totalModes += job->nModes;
                if (job->nSucc)
                    memcpy( job->succ, &succs[3], sizeof(int)*job->nSucc );
                for ( int j=0; (j<job->nSucc); ++j ) {
                    job->succ[j]--;
                    assert( job->succ[j]>=1 && job->succ[j]<jobs );
                    p->jobs[job->succ[j]].nPred++;
                }
            }

            /* allocating pred */
            for ( int i=0; (i<p->nJobs); ++i ) {
                p->jobs[i].pred = NULL;
                if ( p->jobs[i].nPred )
                    ALLOCATE_VECTOR( p->jobs[i].pred, int, p->jobs[i].nPred );
                p->jobs[i].nPred = 0;
            }

            /* filling preds */
            for ( int i=0; (i<p->nJobs); ++i ) {
                for ( int j=0; (j<p->jobs[i].nSucc); ++j ) {
                    Job *jobSucc = &(p->jobs[p->jobs[i].succ[j]]);
                    jobSucc->pred[jobSucc->nPred]= i;
                    jobSucc->nPred++;
                }
            }

            continue;
        }

        if (strstr( s, "mode duration" )) {
            assert( countChar( s, 'R') == nrenew );
            assert( countChar( s, 'N') == nnonrenew );
            assert( firstOccurrence( s, 'N') > lastOccurrence( s, 'R' ) );

            int durLine[ nres+3 ];

            const int idxStartRLocalR  = VInt_size( capR );
            const int idxStartRLocalNR = VInt_size( capNR );
            //const int idxStartRLocalNR = 1;
            inst->idxIniResNProj[p->index] = idxStartRLocalNR;

            s=fgets( line, LINE_SIZE, f);
            assert( s );
            assert( strstr(s, "--") );

            for ( int i=0; i<p->nJobs; ++i ) {
                Job *job = &p->jobs[i];
                job->minDuration = INT_MAX;
                job->maxDuration = 0;
                s=fgets( line, LINE_SIZE, f);
                assert( s );

                readIntVectorStr( s, tok, nres+3, nres+3, durLine );
                job->modes[0].duration = durLine[2];
                assert( nrenew>0 && nnonrenew>=0 );

                Inst_fillModeResUsage( inst, job->modes, nrenew, nnonrenew, durLine+3, idxStartRLocalR, idxStartRLocalNR, globalRCap );

                for ( int j=1; j<job->nModes; ++j ) {
                    s=fgets( line, LINE_SIZE, f);
                    assert( s );
                    readIntVectorStr( s, tok, nres+2, nres+2, durLine );
                    job->modes[j].duration = durLine[1];
                    job->minDuration = MIN( job->minDuration, job->modes[j].duration );
                    job->maxDuration = MAX( job->maxDuration, job->modes[j].duration );

                    Mode *m = job->modes + j;

                    Inst_fillModeResUsage( inst, m, nrenew, nnonrenew, durLine+2, idxStartRLocalR, idxStartRLocalNR, globalRCap );
                }

                inst->nMaxModes = MAX( inst->nMaxModes, job->nModes );
            }


            continue;
        }

        if (strstr( s, "RESOURCEAVAILABILITIES" )) {
            s=fgets( line, LINE_SIZE, f);
            assert( s );
            assert( nres == nrenew+nnonrenew );

#ifdef DEBUG
            /* checking if it is in the format R R (renewable)  followed by non renewable */
            char *pt = NULL;
            char *sp = s;
            for ( int i=0; i<nrenew; ++i ) {
                assert( (pt=strstr(sp,"R")) );
                sp = pt+1;
                pt = NULL;
            }
#endif

            s=fgets( line, LINE_SIZE, f);
            assert( s );

            //int resAV[ nrenew ];
            int *resAV = (int*)malloc( nres * sizeof(int));

            readIntVectorStr( s, tok, nres, nres, resAV );

            for ( int i=0; i<nrenew; ++i )
                if ( globalRCap[i] == -1 )
                    VInt_pushBack( capR, resAV[i] );

            {
                for ( int i=nrenew; i<nres; ++i)
                    VInt_pushBack( capNR, resAV[i] );
            }

            free(resAV);

            continue;
        }

    }

    Tok_free( &tok );

    free( globalRCap );

    /* validating data readed before */
    assert( nres == nrenew+nnonrenew );
    assert( currField == N_FIELDS );
    assert( nprojs == 1 );  /* should be only one project here */
    assert( ndoubly == 0 );

    fclose( f );
}

int readIntVectorStr( const char *str, Tokenizer *tok, int minElements, int maxElements, int *el )
{
    char content[LINE_SIZE];
    strcpy( content, str );
    strRemoveDblSpaces( content );

    Tok_parse( tok, content, ' ');
    int nEl = Tok_nTokens(tok);

    assert( nEl<=maxElements );
    assert( nEl>=minElements );

    for ( int i=0; (i<nEl); ++i ) {
#ifdef DEBUG
        const int l = strlen( Tok_token(tok, i) );
        for ( int j=0; j<l; ++j )
            assert( isdigit( Tok_token(tok, i)[j] ) );
#endif

        el[i] = atoi( Tok_token(tok, i) );
    }

    return nEl;
}

static void Inst_freeJobVectorContents(Instance * inst, Job *jStart, Job *jEnd )
{
    Job *j = jStart;

    for (; j<jEnd; ++j ) {
        Mode *m = j->modes;
        Mode *mEnd = j->modes + j->nModes;
        for (; (m<mEnd); ++m ) {
            if (m->idxResR) {
                free( m->idxResR );
                free( m->useResR );
            }
            if (m->idxResN) {
                free( m->idxResN );
                free( m->useResN );
            }

        }

        free( j->modes );
        free (j->idxModesSort);

        if ( j->succ )
            free( j->succ );
        if ( j->pred )
            free( j->pred );
        if ( j->hasSucc )
            free( j->hasSucc );
        if ( j->hasPred )
            free( j->hasPred );
        if ( j->hasIndSucc )
            free( j->hasIndSucc );
        if ( j->hasIndPred )
            free( j->hasIndPred );


    }
}

void Inst_free( Instance **_inst )
{
    Instance *inst = *_inst;

    for(int i  = 0 ; i < inst->nJobs ; i ++) {
        free(inst->matMaxD[i]);
        VInt_free(&inst->paths[i]);
        const Job * job = Inst_job(inst,i);
        int nMode = Job_nModes(job);
        for(int m  = 0 ; m < nMode ; m ++) {
            free(inst->matMaxDJM[i][m]);
            for(int j  = 0 ; j < inst->nJobs ; j ++)
                free(inst->maxDurationPathInterMode[i][m][j]);
            free(inst->maxDurationPathInterMode[i][m]);
        }
        free(inst->maxDurationPathInterMode[i]);
        free(inst->maxDurationPath[i]);
        free(inst->matMaxDJM[i]);
    }

    free(inst->paths);
    free(inst->maxDurationPath);
    free(inst->maxDurationPathInterMode);
    free(inst->matMaxD);
    free(inst->matMaxDJM);

    for ( int i=0; i<inst->nProjects; ++i ) {
        Inst_freeJobVectorContents( inst, inst->projects[i].jobs, inst->projects[i].jobs + inst->projects[i].nJobs );
        free( inst->projects[i].jobs );
    }

    Inst_freeJobVectorContents( inst, inst->jobs, inst->jobs+inst->nJobs );



    free( inst->capR );
    free( inst->capN );
    free( inst->jobs );
    free( inst->projects );
    free( inst->idxIniResNProj );
    free( inst );

    *_inst = NULL;
}

int Inst_nJobs( const Instance *inst )
{
    assert( inst!=NULL );

    return inst->nJobs;
}

const Job *Inst_job( const Instance *inst, int idxJob )
{
    assert( inst!=NULL );
    assert( idxJob >= 0 );

    //   printf("idxJob %d inst->nJobs %d\n\n", idxJob , inst->nJobs);
    assert( idxJob < inst->nJobs );

    return inst->jobs + idxJob;
}

int Inst_nResR( const Instance *inst )
{
    assert( inst );

    return inst->nResR;
}

int Inst_nResN( const Instance *inst )
{
    assert( inst );

    return inst->nResN;
}

int Inst_capResR( const Instance *inst, int i )
{
    assert( inst );
    assert( i>=0 && i<inst->nResR );

    return inst->capR[i];
}

int Inst_capResN( const Instance *inst, int i )
{
    assert( inst );
    assert( i>=0 && i< Inst_nResN(inst) );

    return inst->capN[i];
}

int Mode_duration( const Mode *mode )
{
    assert( mode!=NULL );


    return mode->duration;
}

char Mode_isFeasible( const Instance *inst, const Mode *mode )
{
    assert( inst );
    assert( mode );
    //printf("\nINI\n")); getchar();
    for ( int i=0; (i<mode->nResR); ++i )
        if ( Inst_capResR( inst, mode->idxResR[i] ) < mode->useResR[i] )
            return False;
    //printf("\nMEIO\n"); getchar();
    for ( int i=0; (i<mode->nResN); ++i )
        if ( Inst_capResN( inst, mode->idxResN[i] ) < mode->useResN[i] )
            return False;

    //printf("\nFIM\n"); getchar();
    return True;
}

void Inst_jobsInfeasible(const Instance *inst)
{
    for(int j1 = 1 ; j1 < Inst_nJobs(inst) ; j1++) {
        const Job *job1 = Inst_job(inst,j1);
        int m1 = 0;
        const Mode *mode1 = Job_mode(job1,m1);
        for(int r = 0 ; r < Inst_nResR(inst) ; r++) {
            //printf("job %d  mode %d  resource %d\n", j1,m1,r);
            int aux = 0;
            int r1m = Mode_idxResROnMode(inst,mode1,r);
            if(r1m == -1) {
                //   printf("    job %d  mode %d  resource %d = %d\n", j1,m1,r, r1m); getchar();
                continue;
            }
            for(int j2 = 1 ; j2 < Inst_nJobs(inst); j2++) {
                if(j1==j2 || Job_hasIndSucc(job1,j2) || Job_hasIndSucc(job1,j2)) continue;
                const Job *job2 = Inst_job(inst,j2);
                int m2 = 0;
                const Mode *mode2= Job_mode(job2,m2);
                int r2m = Mode_idxResROnMode(inst,mode2,r);
                // printf("    job %d  mode %d  resource %d\n", j2,m2,r);
                if(r2m == -1) {
                    //printf("    job %d  mode %d  resource %d = %d\n", j2,m2,r,r2m); getchar();
                    aux = 0;
                    break;
                } else {
                    if(Mode_useResR(mode1,r1m)+Mode_useResR(mode2,r2m)> Inst_capResR(inst,r))
                        aux++;
                    else {
                        // printf("    job %d  mode %d and job %d mode % do not have conflicts in resource %d (cap %d, totaluses %d) \n", j1,m1,j2,m2,r,Inst_capResR(inst,r),Mode_useResR(mode1,r1m)+Mode_useResR(mode2,r2m)); getchar();
                        aux = 0;
                        break;
                    }
                }
            }
            if(aux>0) {
                printf("job %d in mode %d is enfeasible with all others to resource %d\n", j1,m1,r);//getchar();
                inst->jobmodeinfeasible[j1][m1] = 1;
                break;
            }
        }
    }

}

static void Mode_cpy( Mode *mTarget, const Mode *mSource)
{
    assert( mTarget );
    assert( mSource );
    assert( mTarget->idxResR == NULL );
    assert( mTarget->useResR == NULL );
    assert( mTarget->idxResN == NULL );
    assert( mTarget->useResN == NULL );

    memcpy( mTarget, mSource, sizeof(Mode) );

    if ( mSource->nResR ) {
        ALLOCATE_VECTOR( mTarget->idxResR, int, mSource->nResR );
        ALLOCATE_VECTOR( mTarget->useResR, int, mSource->nResR );

        memcpy( mTarget->idxResR, mSource->idxResR, sizeof(int)*mSource->nResR );
        memcpy( mTarget->useResR, mSource->useResR, sizeof(int)*mSource->nResR );
    } else {
        mTarget->idxResR = NULL;
        mTarget->useResR = NULL;
    }


    if ( mSource->nResN ) {
        ALLOCATE_VECTOR( mTarget->idxResN, int, mSource->nResN );
        ALLOCATE_VECTOR( mTarget->useResN, int, mSource->nResN );

        memcpy( mTarget->idxResN, mSource->idxResN, sizeof(int)*mSource->nResN );
        memcpy( mTarget->useResN, mSource->useResN, sizeof(int)*mSource->nResN );
    } else {
        mTarget->idxResN = NULL;
        mTarget->useResN = NULL;
    }
}

/* number of renewable resources used by this mode */
int Mode_nResR( const Mode *mode )
{
    assert( mode );
    return mode->nResR;
}

/* indexes of renewable resources used
 * with this mode */
int Mode_idxResR( const Mode *mode, int i )
{
    assert( mode );
    assert( i>=0 && i<mode->nResR );

    return mode->idxResR[i];
}


/* indexes of renewable resources of instance on mode i == index of resource on instance*/
int Mode_idxResROnMode( const Instance *inst, const Mode *mode, int i )
{
    assert( mode );
    assert( inst );
    assert( i>=0 && i< inst->nResR );

    for(int r = 0 ; r < mode->nResR; r++) {
        if(mode->idxResR[r] == i)
            return r;
    }

    return -1;

}

/* use of renewable resources */
int Mode_useResR( const Mode *mode, int i )
{
    assert( mode );
    assert( i>=0 && i<mode->nResR );

    return mode->useResR[i];
}

/* number of non renewable resources used by a mode */
int Mode_nResN( const Mode *mode )
{
    assert( mode );

    return mode->nResN;
}

/* to access indexes of local non-renewable resources used in a mode */
int Mode_idxResN( const Mode *mode, int i )
{
    assert( mode );
    assert( i>=0 && i<mode->nResN );

    return mode->idxResN[i];
}

/* to access the use of local non-renewable resources used in a mode */
int Mode_useResN( const Mode *mode, int i )
{
    assert( mode );
    assert( i>=0 && i<mode->nResN );

    return mode->useResN[i];
}

int Mode_index( const Mode *mode )
{
    assert( mode );

    return mode->index;
}

int Inst_nProjects( const Instance *inst )
{
    assert( inst!=NULL );

    return inst->nProjects;
}

int Inst_idxResNProj( const Instance *inst, int idxProj )
{
    assert( inst!=NULL );
    assert( idxProj >= 0);
    assert( idxProj < inst->nProjects);

    return inst->idxIniResNProj[idxProj];
}

const Project *Inst_project( const Instance *inst, int idxProject )
{
    assert( inst!=NULL );
    assert( idxProject >= 0 );
    assert( idxProject < inst->nProjects );

    return inst->projects + idxProject;
}

int Project_index( const Project *p )
{
    assert( p!=NULL );

    return p->index;
}

int Project_releaseDate( const Project *p )
{
    assert( p!=NULL );

    return p->releaseDate;
}

int Project_dueDate( const Project *p )
{
    assert( p!=NULL );

    return p->dueDate;
}

int Project_tardCost( const Project *p )
{
    assert( p!=NULL );

    return p->tardCost;
}

int Project_mpmTime( const Project *p )
{
    assert( p!=NULL );

    return p->mpmTime;
}

int Project_criticalPath( const Project *p )
{
    assert( p!=NULL );

    return p->criticalPath;
}

int Project_calc_criticalPath(Instance *inst)
{
    int critPathValue=0;

    int **vetVertices;

    vetVertices = (int**)malloc( inst->projects[0].nJobs * sizeof(int*) );

    for(int i = 0; i < inst->projects[0].nJobs; ++i)
        vetVertices[i] = (int*)malloc( 2 * sizeof(int) );

    /* initial vertex */
    vetVertices[0][0] = 0;
    vetVertices[0][1] = 0;

    for(int h = 1; h < inst->projects[0].nJobs; ++h) {
        vetVertices[h][0] = 0;
        vetVertices[h][1] = INT_MAX;
    }

    /* ida no grafo */
    for(int i = 1; i < inst->projects[0].nJobs; ++i) {
        for(int j = 0; j < inst->projects[0].jobs[i].nPred; ++j) {
            if( vetVertices[inst->projects[0].jobs[i].pred[j]][1] >= vetVertices[i][0]) {
                vetVertices[i][0] = vetVertices[inst->projects[0].jobs[i].pred[j]][1];
                vetVertices[i][1] = vetVertices[i][0] + inst->jobs[i].minDuration;
                //printf("\njob_%d [%d|%d|%d]", i, vetVertices[i][0], inst->jobs[i].minDuration, vetVertices[i][1]); getchar();
                critPathValue = vetVertices[i][1];
            }
        }
    }

    for(int i = 0; i < inst->projects[0].nJobs; ++i)
        free( vetVertices[i] );

    free( vetVertices );

    return critPathValue;
}

int Project_nJobs( const Project *p )
{
    assert( p!=NULL );

    return p->nJobs;
}

int Project_idxFirstJob( const Project *p )
{
    assert( p!=NULL );

    return p->idxFirstJob;
}

const Mode *Job_mode( const Job *job, int idxMode )
{
    assert( job!=NULL );
    assert( idxMode >= 0 );

    assert( idxMode < Job_nModes(job));

    return job->modes + idxMode;
}

int Job_project( const Job *job)
{
    assert( job!=NULL );

    return job->idxProject;
}

int Job_nModes( const Job *job )
{
    assert( job!=NULL );

    return job->nModes;
}

int Job_nInfeasModes( const Job *job )
{
    assert( job!=NULL );

    return job->nInfeasModes;
}

const Job *Project_job( const Project *p, int idxJob )
{
    assert( p!=NULL );
    assert( idxJob >= 0 );
    assert( idxJob < p->nJobs );

    return p->jobs + idxJob;
}

void Inst_setSumTPD(Instance *inst, const int value)
{
    inst->sumTPD = value;
}

int Inst_getSumTPD(const Instance *inst)
{
    return inst->sumTPD;
}

void Inst_setSumTMS(Instance *inst, const int value)
{
    inst->sumTMS = value;
}

int Inst_getSumTMS(const Instance *inst)
{
    return inst->sumTMS;
}


int Inst_getMaxDIJ(const Instance *inst, int i, int j)
{
    return inst->matMaxD[i][j];
}

int Inst_getMaxDIJM(const Instance *inst, int i, int j, int m)
{
    return inst->matMaxDJM[i][m][j];
}


int Inst_getMaxDIM(const Instance *inst, const Job *job,int m)
{
    return inst->maxDurationPath[Job_index(job)][m];
}


int Inst_getMaxDIMJM(const Instance *inst, const Job *job, int m,  int j2, int m2)
{
    return inst->maxDurationPathInterMode[Job_index(job)][m][j2][m2];
}

int Job_nSucc( const Job *job )
{
    assert( job!=NULL );

    return job->nSucc;
}

int Job_succ( const Job *job, int idxSucc )
{
    assert( job!=NULL );
    assert( idxSucc >= 0 );
    assert( idxSucc < job->nSucc );

    return job->succ[idxSucc];
}

int Job_nPred( const Job *job )
{
    assert( job!=NULL );

    return job->nPred;
}

int Job_pred( const Job *job, int idxPred )
{
    assert( job!=NULL );
    assert( idxPred >= 0 );
    assert( idxPred < job->nPred );

    return job->pred[idxPred];
}

int Job_hasIndPred(  const Job *job, int idxJob )
{
    assert( job!=NULL );
    assert( idxJob >= 0 );

    return job->hasIndPred[idxJob];
}

int Job_hasIndSucc(  const Job *job, int idxJob )
{
    assert( job!=NULL );
    assert( idxJob >= 0 );

    return job->hasIndSucc[idxJob];
}

int Job_hasPred( const Instance *inst, const Job *job, int idxJob )
{
    assert( job!=NULL );
    assert( idxJob >= 0 );

    assert( idxJob < Inst_nJobs(inst) );

    for(int p = 0; p < Job_nPred(job); ++p) {
        int idJob = Job_pred(job,p);
        if(idJob == idxJob) return 1;
        else  return 0;
    }
    return 0;
}

int Job_hasSucc( const Instance *inst, const Job *job, int idxJob )
{
    assert( job!=NULL );
    assert( idxJob >= 0 );
    assert( idxJob < Inst_nJobs(inst) );

    for(int s = 0; s < Job_nSucc(job); ++s) {
        int idJob = Job_succ(job,s);
        if(idJob == idxJob) return 1;
        // else   return 0;
    }
    return 0;
}

int Job_minDuration( const Job *job )
{
    assert( job!=NULL );

    return job->minDuration;
}

int Job_idxModeMinDuration( const Job *job )
{
    assert( job!=NULL );

    return job->idxModeMinDuration;
}

int Job_maxDuration( const Job *job )
{
    assert( job!=NULL );

    return job->maxDuration;
}

int Job_est( const Job *job )
{
    assert( job!=NULL );

    return job->est;
}

int Job_index( const Job *job )
{
    assert( job!=NULL );

    return job->index;
}

int Job_idxOnProject( const Job *job )
{
    assert( job!=NULL );

    return job->idxOnProject;
}

int Inst_getSizePath(const Instance *inst, int i)
{
    return VInt_size(inst->paths[i]);
}


int Inst_getValuePosPath(const Instance *inst, int i, int pos)
{
    return VInt_get(inst->paths[i], pos);
}

static void Job_cpy( const Instance *inst, Job *jTarget, const Job* jSource, int* nRemovedModes )
{
    assert( jTarget != NULL );
    assert( jTarget->modes == NULL );
    assert( jTarget->succ == NULL );
    assert( jTarget->pred == NULL );

    assert( jSource != NULL );
    assert( jSource->nModes > 0 );

    /* counting feasible jobs in source */
    int nFeasModes = 0, nInfeasModes=0;
    for ( int i=0; (i<Job_nModes(jSource)); ++i ) {
        nFeasModes ++;//= (int) Mode_isFeasible( inst, Job_mode(jSource,i) );
        if(!Mode_isFeasible( inst, Job_mode(jSource,i) )) nInfeasModes++;
    }

    assert(nFeasModes>0);
    /* copying modes */
    ALLOCATE_VECTOR_INI( jTarget->modes, Mode, nFeasModes );

    int idxMode = 0;
    for ( int i=0; (i<Job_nModes(jSource)); ++i ) {
        const Mode *modeSource = Job_mode( jSource, i );
        // if (!Mode_isFeasible( inst, modeSource )) {
        //   jTarget->idxModes[i] = -1;
        // continue;
        //}

        Mode_cpy( jTarget->modes+idxMode, modeSource );
        ++idxMode;
    }
    jTarget->nModes = nFeasModes;
    jTarget->nInfeasModes = nInfeasModes;
    nRemovedModes += Job_nModes(jSource) - nFeasModes;

    jTarget->nSucc = jSource->nSucc;
    jTarget->nPred = jSource->nPred;

    if (jSource->nSucc) {
        ALLOCATE_VECTOR( jTarget->succ, int, jSource->nSucc );
        COPY_VECTOR( jTarget->succ, jSource->succ, int, jSource->nSucc );
    } else
        jTarget->succ = NULL;

    if (jSource->nPred) {
        ALLOCATE_VECTOR( jTarget->pred, int, jSource->nPred );
        COPY_VECTOR( jTarget->pred, jSource->pred, int, jSource->nPred );
    } else
        jTarget->pred = NULL;
}

void Inst_print( const Instance *inst )
{
    printf("[Instance]:\n\nnResR=%d\t\tnResN=%d\n", Inst_nResR(inst), Inst_nResN(inst) );

    printf("\n*capR:\t");
    for ( int i=0; i<Inst_nResR( inst ); ++i )
        printf("%d\t", i );
    printf("\n\t");

    for ( int i=0; i<Inst_nResR( inst ); ++i )
        printf("%d\t", Inst_capResR( inst, i ) );

    printf("\n\n*capN:\t" );
    for ( int i=0; i<Inst_nResN( inst ); ++i )
        printf("%d\t", i );

    printf("\n\t");

    for ( int i=0; i<Inst_nResN( inst ); ++i )
        printf("%d\t", Inst_capResN( inst, i) );

    for ( int i=0; i<Inst_nProjects( inst ); ++i )
        Project_print( inst, i );

}

void Project_print(  const Instance *inst, int p )
{
    const Project *proj = Inst_project( inst, p );

    printf("\n\n[Project %d]: \n\nindex    crit. path    rel. Date    nJobs    dueDate    tardCost    mpmTime\n", p );
    printf(" %d\t  %d\t\t%d\t     %d\t      %d\t %d\t     %d", Project_index( proj ), Project_criticalPath( proj ), Project_releaseDate( proj ), Project_nJobs( proj ), Project_dueDate( proj ), Project_tardCost( proj ), Project_mpmTime( proj ) );

    printf("\n\n[Jobs P%d]:\n\nindex\tidxPrj\tnMods\tnSucc\tnPred\tminDur\tmaxDur\test\tsuccessors\n", p );
    int startProj=0;
    for ( int i=0; i< Inst_nJobs( inst ); ++i) {
        if( Job_project( Inst_job( inst, i ) )==p ) {
            startProj = i;
            break;
        }
    }

    for ( int i=startProj; i<(startProj+Project_nJobs( proj )); ++i)
        Job_print( inst, i );

    printf("\n[Modes/Jobs P%d]:\n\nidxJob\tidxMode\t  duration\tidxResR\tuseResR\tidxResN\tuseResN\n", p );
    for ( int j=startProj; j<(startProj+Project_nJobs( proj )); ++j ) {
        for ( int i=0; i<Job_nModes( Inst_job( inst, j) ); ++i ) {
            if ( i==0 )
                Mode_print( inst, Inst_job( inst, j), i, 0 );
            else
                Mode_print( inst, Inst_job( inst, j), i, 1 );
        }
    }

    printf("\n");
}

void Job_print( const Instance *inst, int i )
{
    const Job *job = Inst_job( inst, i );

    printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t", Job_index( job ), Job_idxOnProject( job ), Job_nModes( job ), Job_nSucc( job ), Job_nPred( job ), Job_minDuration( job ), Job_maxDuration( job ), Job_est( job ) );

    for ( int i=0; i<Job_nSucc( job ); ++i )
        printf("%d  ", Job_succ( job, i ) );

    printf("\n");
}

void Mode_print( const Instance *inst, const Job *job, int i, int nHid )
{

    const Mode* mode = Job_mode( job, i );

    if ( nHid==0 )
        printf("\n%d\t%d\t  %d\t\t", Job_index( job ), Mode_index( mode ), Mode_duration( mode ) );
    else
        printf("\n \t%d\t  %d\t\t", Mode_index( mode ), Mode_duration( mode ) );

    for ( int i=0; i<Mode_nResR( mode ); ++i )
        printf("%d ", Mode_idxResR( mode, i ) );
    printf("\t");
    for ( int i=0; i<Mode_nResR( mode ); ++i )
        printf("%d ", Mode_useResR( mode, i ) );
    printf("\t");
    for ( int i=0; i<Mode_nResN( mode ); ++i )
        printf("%d ", Mode_idxResN( mode, i ) );
    printf("\t");
    for ( int i=0; i<Mode_nResN( mode ); ++i )
        printf("%d ", Mode_useResN( mode, i ) );
}

void Inst_computeEST( Instance *inst )
{
    /* not assuming that dummy jobs are included, searching for jobs without dependencies */
    for ( int p=0 ; (p<Inst_nProjects(inst)) ; ++p ) {
        const Project *project = Inst_project( inst, p );
        for ( int j=Project_idxFirstJob(project) ; (j<Project_idxFirstJob(project)+Project_nJobs(project)) ; ++j ) {
            Job *job = (Job*) Inst_job( inst, j );
            if ( Job_nPred(job) > 0 )
                continue;
            job->est = Project_releaseDate( project );
            Inst_computeESTJob( inst, j );
        }
    }
}

void Inst_computeESTJob( Instance *inst, int jIdx )
{
    const Job *job = (Job*) Inst_job( inst, jIdx );

    for ( int i=0 ; (i<Job_nSucc(job)) ; ++i ) {
        const int jSuccIdx = Job_succ( job, i );
        Job *jSucc = (Job*) Inst_job( inst, jSuccIdx );
        jSucc->est = MAX( jSucc->est, job->est+Job_minDuration(job) );
        Inst_computeESTJob( inst, jSuccIdx );
    }
}

int Inst_nResRGlobal( const Instance *inst )
{
    return inst->nResRGlobal;
}

int somaMinDurationPathsByJob(const Instance *inst, int j, int nSuccs)
{
    assert(inst);

    const Job *job = Inst_job(inst,j);
    /*           printf("\nJob %d succs:", j);

       for(int s = 0; s< Inst_nJobs(inst) ; s++) {
               if(Job_hasIndSucc(job,s))
                   printf("%d, ", s);
       }*/
    // getchar();

    if(nSuccs==0) {
        //   printf("Dummy Job_minDuration(%d): %d\n", Job_index(job), Job_minDuration(job));
        return Job_minDuration(job);
    } else {
        int maxD=0;
        for(int s = 0; s< Job_nSucc(job) ; s++) {
            int index = Job_succ(job,s);
            const Job* jobSucc = Inst_job(inst,index);
            int nS = Job_nSucc( jobSucc);
            //      printf("JobSuc %d nS %d\n", index, nS);
            int soma = somaMinDurationPathsByJob(inst, index, nS);
            int value = Job_minDuration(jobSucc)+soma;
            //     printf("Job %d MinDur %d somaMinDurPath %d Value %d\n", Job_index(jobSucc), Job_minDuration(jobSucc), soma, value);
            maxD = MAX(maxD,value);
        }
        return maxD;
    }
}


int findJobOnPath(const Instance *inst, int s, int nSuccs, VecInt **paths, int j)
{
    assert(inst);

    const Job *job = Inst_job(inst,s);
    // VInt_pushBack(paths[j], s);

    if(nSuccs==0) {
        // printf(" %d ", Job_index(job));
        VInt_pushBack(paths[j], Job_index(job));
        return  1;
    } else {
        for(int ss = 0; ss < Job_nSucc(job) ; ss++) {
            int index = Job_succ(job,ss);
            const Job* jobSucc = Inst_job(inst,index);
            int nS = Job_nSucc( jobSucc);
            //   printf("%d ", index);

            VInt_pushBack(paths[j], s);

            findJobOnPath(inst, index, nS, paths, j);
            // return index;
        }
    }
    return 1;
}

void Inst_computeCompPaths( const Instance *inst, VecInt **paths)
{

    for(int j = 0 ; j < Inst_nJobs(inst) ; j++) {
        const Job *job = Inst_job(inst, j);
        int nSucc = Job_nSucc(job);
        //   printf("\nJob %d: \n ", j);
        findJobOnPath(inst, j, nSucc, paths,j);
        //VInt_pushBack(paths[j],value);
    }

}

int Inst_computeCompPathsByJob( const Instance *inst, int i, int m,  int d)
{
    const Job *job = Inst_job(inst,i);
    const Project *project = Inst_project(inst,Job_project(job));

    //int sumDurations = Job_getMaxDIM(job,m);
    //  int sumDurationsold = d+somaMinDurationPathsByJob(inst, i, Job_nSucc(job));

    int maxD = Project_releaseDate(project)+Project_criticalPath(project) - Inst_getMaxDIM(inst,job,m);
    // printf("sumDurations %d,sumDurationsold %d maxD %d \n", sumDurations, sumDurationsold, maxD);
    //printf("J %d D %d | maxD = %d => Project_releaseDate(project) %d + Project_criticalPath(project)%d - sumDurations %d \n ", i, d, maxD, Project_releaseDate(project), Project_criticalPath(project), sumDurations);
    return maxD;
}

int somaMinDurationPathsByJobAndInter(const Instance *inst, int j, int nSuccs, int inter, int durInter)
{
    assert(inst);

    const Job *job = Inst_job(inst,j);
    /*           printf("\nJob %d succs:", j);

       for(int s = 0; s< Inst_nJobs(inst) ; s++) {
               if(Job_hasIndSucc(job,s))
                   printf("%d, ", s);
       }*/
    // getchar();

    if(nSuccs==0) {
        //   printf("Dummy Job_minDuration(%d): %d\n", Job_index(job), Job_minDuration(job));
        return Job_minDuration(job);
    } else {
        int maxD=0;
        for(int s = 0; s< Job_nSucc(job) ; s++) {
            int index = Job_succ(job,s);
            const Job* jobSucc = Inst_job(inst,index);
            int nS = Job_nSucc( jobSucc);
            //      printf("JobSuc %d nS %d\n", index, nS);
            int soma = somaMinDurationPathsByJobAndInter(inst, index, nS, inter, durInter);
            int value =0;
            if(index == inter)
                value = durInter+soma;
            else
                value = Job_minDuration(jobSucc)+soma;
            //     printf("Job %d MinDur %d somaMinDurPath %d Value %d\n", Job_index(jobSucc), Job_minDuration(jobSucc), soma, value);
            maxD = MAX(maxD,value);
        }
        return maxD;
    }
}

int Inst_computeCompPathsByJobAndInter_old( const Instance *inst, int i, int d, int s, int ds )
{
    const Job *job = Inst_job(inst,i);
    const Project *project = Inst_project(inst,Job_project(job));

    int sumDurations = d+somaMinDurationPathsByJobAndInter(inst, i, Job_nSucc(job), s, ds);
    int maxD = Project_releaseDate(project)+Project_criticalPath(project) - sumDurations;
    printf("I sumDurations %d, maxD %d \n", sumDurations, maxD);
    //printf("J %d D %d | maxD = %d => Project_releaseDate(project) %d + Project_criticalPath(project)%d - sumDurations %d \n ", i, d, maxD, Project_releaseDate(project), Project_criticalPath(project), sumDurations);
    return maxD;
}

int Inst_computeCompPathsByJobAndInter( const Instance *inst, int i, int m, int d, int s, int m2, int ds, int dsmin )
{
    const Job *job = Inst_job(inst,i);
    // printf("job %d, m %d, s %d, m2 %d\n", i,m,s,m2);fflush(stdin);
    //int sumDurationsold = d+somaMinDurationPathsByJobAndInter(inst, i, Job_nSucc(job), s, ds);
    const Project *project = Inst_project(inst,Job_project(job));

    int maxD = Project_releaseDate(project)+Project_criticalPath(project) -  Inst_getMaxDIMJM(inst,job,m,s,m2);
    //   printf("I sumDurations %d, sumDurationsold %d, maxD %d \n", sumDurations,sumDurationsold, maxD);
    //printf("J %d D %d | maxD = %d => Project_releaseDate(project) %d + Project_criticalPath(project)%d - sumDurations %d \n ", i, d, maxD, Project_releaseDate(project), Project_criticalPath(project), sumDurations);
    return maxD;
}

void floydWarshallMax(const Instance *inst, int **D, int nJobs)
{
    assert(inst);

    for(int i = 0 ; i < nJobs ;  i++) {
        const Job* job = Inst_job(inst, i);
        for(int j = 0 ; j < nJobs ; j++) {
            if(Job_hasSucc(inst,job,j)) {
                //   printf("\ni %d,j %d, nJobs %d\n",i,j,nJobs);
                D[i][j] = Job_minDuration(job);
            }
        }
    }

    for(int k = 0; k < Inst_nJobs(inst) ; k++) {
        const Job * jobInter = Inst_job(inst, k);
        for(int p = 0; p < Inst_nJobs(inst) ; p++) {
            if(Job_hasIndPred(jobInter,p)) {
                const Job * jobPred = Inst_job(inst,p);
                for(int s = 0; s < Inst_nJobs(inst) ; s++) {
                    if(Job_hasIndSucc(jobPred,s) && Job_hasIndSucc(jobInter,s)) {
                        //  printf("(%d,%d,%d)", p,s);
                        if(D[p][s] < (D[p][k]+D[k][s]))
                            D[p][s] = (D[p][k]+D[k][s]);
                    }
                }
            }
        }
    }
}

void maxDistanceByModes(const Instance *inst, int **D, int *** DJM)
{

    for(int i  = 0 ; i < inst->nJobs ; i ++) {
        const Job *job = Inst_job(inst, i);
        int nMode = Job_nModes(job);
        for(int m  = 0 ; m < nMode ; m ++) {
            const Mode * mode = Job_mode(job,m);
            int plus = Mode_duration(mode) - Job_minDuration(job);
            for(int j  = 0 ; j < inst->nJobs ; j ++) {
                if(Job_hasIndSucc(job,j)) {
                    DJM[i][m][j] = D[i][j]+plus;
                    //                    printf(" D[%d][%d] = %d, DJM[%d][%d][%d] = %d plus %d\n", i, j, D[i][j], i,m,j, DJM[i][m][j], plus);
                }
            }
        }
    }
}


/*
void computeUnitsResourcesES(const Instance *inst, int nTimes, int **maxTJM )
{

    int nJobs = Inst_nJobs(inst);
    int nModes = Inst_nMaxModes(inst);

    for(int j = 0 ; j < nJobs ; j++){
        const Job* job = Inst_job(inst,j);
        ALLOCATE_VECTOR(inst->unitsResourceES[j],int***,nModes);
        for(int m = 0 ; m < Job_nModes(job) ; m++){
            const Mode *mode = Job_mode(job,m);
            ALLOCATE_VECTOR(inst->unitsResourceES[j][m],int**,nTimes);
            for(int s = 0 ; s < nTimes ; s++){
                ALLOCATE_VECTOR_INI(inst->unitsResourceES[j][m][s],int*,nTimes);
                for(int e = s+1 ; e < nTimes ; e++){
                    int part1 = e-s;
                    int part2 = MAX(0,Mode_duration(mode)-MAX(0,s-Job_est(job)));
                    int part3 = MAX(0,Mode_duration(mode)-MAX(0,maxTJM[j][m]));
                    inst->unitsResourceES[j][m][s][e] = MIN(part1,MIN(part2,part3));
                  //  if( inst->unitsResourceES[j][m][s][e]>0)
               //     printf("inst->unitsResourceES[%d][%d][%d][%d] %d\n", j,m,s,e, inst->unitsResourceES[j][m][s][e]);
                }

            }
        }
    }

}

int Inst_unitsResourceES( const Instance *inst, int j, int m, int s, int e)
{

    return inst->unitsResourceES[j][m][s][e];
}*/


int Inst_nMaxModes( const Instance *inst )
{
    return inst->nMaxModes;
}

int Inst_fileExists(char fileName[])
{
    FILE *fp;

    fp=fopen(fileName,"r");

    if(fp) {
        fclose(fp);
        return 1;
    } else
        return 0;

    return -1;
}

