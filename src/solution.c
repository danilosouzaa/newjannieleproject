
/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <string.h>
#include <assert.h>
#include "solution.h"
#include "instance.h"
#include "macros.h"
#include "rrusage.h"
//#include "neighborhood.h"
#include "node_heap.h"
#include <stdarg.h>


struct _Solution {

    Cost cost;
    Cost TPD;
    Cost TMS;

    int *startJobs;
    int *sequence;
    int *posJobs;
    unsigned int solutionHash;

    NodeHeap *priorities;

    RRUsage *rru;

    int *minT;
    int *origMinT;

    struct _ModeSet *modeSet;
    const struct _Instance *inst;
};

int Sol_read( Solution *sol, char *file )
{
    assert( sol != NULL );

    FILE *fp = fopen( file, "r" );
    int read=1;

    if ( fp == NULL ) {
        printf( "File was not opened(initSol). path: %s\n", file);
        read=0;
        return read;
    }

    int p, j, m, t;
    int cont = 0;

    unsigned int hash=0;

    while (fscanf(fp, "%d", &p) == 1
           && fscanf(fp, "%d", &j)
           && fscanf(fp, "%d", &m)
           && fscanf(fp, "%d", &t)) {

        const Job *job;
        if(p==0) job = Inst_job(sol->inst, j);
        else job = Inst_job(sol->inst, cont);

        const Mode *mode = Job_mode( job, m );

        Modes_modify(sol->modeSet,Job_index(job), m);

        RRU_allocate(sol->rru, mode, t);
        Sol_setStartJob(sol,Job_index(job),t);

        hash = ((hash << 5) - hash) + t+Mode_duration(mode);
        cont++;


    }

    sol->solutionHash = hash;

    Sol_fillSequence(sol);

    Modes_isntEmpty(sol->modeSet);

    Sol_calcCost(sol);

    fclose(fp);
    return read;

}

void Sol_fillSequence( Solution *sol )
{

    assert( sol != NULL );

    NodeHeap *seqStart = nh_create(Inst_nJobs( sol->inst ), INT_MAX_M);

    int ret;

    nh_reset(seqStart);

    for ( int i=0; (i<Inst_nJobs( sol->inst )); i++ ) {
        const Job *job = Inst_job(sol->inst,i);
        int value = sol->startJobs[Job_index(job)];
        nh_update(seqStart, i, value);
    }

    for ( int i=0; (i<Inst_nJobs( sol->inst )); i++ ) {
        nh_remove_first(seqStart, &ret);
        sol->sequence[i] = ret;
        sol->posJobs[ret] = i;
    }

    nh_free(&seqStart);

}

Solution *Sol_create( const Instance *inst )
{
    assert( inst != NULL );

    Solution* sol;

    ALLOCATE_INI( sol, Solution );
    ALLOCATE_VECTOR_INI( sol->startJobs, int, Inst_nJobs( inst ) );
    ALLOCATE_VECTOR_INI( sol->sequence, int, Inst_nJobs( inst ) );
    ALLOCATE_VECTOR_INI( sol->posJobs, int, Inst_nJobs( inst ) );
    ALLOCATE_VECTOR_INI( sol->minT, int, Inst_nJobs( inst ) );
    ALLOCATE_VECTOR_INI( sol->origMinT, int, Inst_nJobs( inst ) );

    sol->inst = inst;
    sol->priorities = nh_create(Inst_nJobs( inst ), INT_MAX_M);
    sol->modeSet = Modes_create( inst );
    sol->rru = RRU_create(inst);
    sol->solutionHash = 0;

    for (int j = 0; j<(Inst_nJobs( inst ) ); j++ ) {
        const Job* job = Inst_job(inst, j);
        sol->origMinT[j] =  Job_est(job);
    }

    memcpy( sol->minT, sol->origMinT, sizeof(int)*Inst_nJobs( inst ) );


    return sol;
}

void Sol_topSort( Solution *solution, int sequence[] )
{

    assert( solution != NULL );

    int ind, costInd;
    int ret;

    nh_reset(solution->priorities);

    for ( int i=0; (i<Inst_nJobs( solution->inst )); i++ ) {
        int j;
        j = sequence[i];
        const Job* job = Inst_job( solution->inst, j );
        nh_update(solution->priorities, Job_index(job), i+Job_nPred(job) * Inst_nJobs(solution->inst));
    }

    for ( int i=0; (i<Inst_nJobs( solution->inst )); i++ ) {

        // job which will be allocated
        nh_remove_first(solution->priorities, &ret);
        sequence[i] = ret;
        solution->posJobs[ret] = i;
        const Job* job = Inst_job(solution->inst, ret);
        for(int h = 0; h < Job_nSucc(job); ++h) {
            const Job* job2 =  Inst_job(solution->inst, Job_succ(job,h));
            ind = Job_index(job2);
            costInd = nh_get_dist(solution->priorities, ind);
            nh_update(solution->priorities, ind, costInd - Inst_nJobs(solution->inst));
        }
    }
}

void Sol_setStartJob(Solution *sol, int idxJob, int time)
{

    assert( sol != NULL );
    assert( idxJob >= 0 );
    assert( idxJob < Inst_nJobs(sol->inst) );

    sol->startJobs[idxJob] = time;
}

void Sol_rebuild( Solution *solution)
{

    assert(solution != NULL);

    Sol_clearBuildSequence(solution);

    int nJobs = Inst_nJobs( solution->inst );
    memcpy( solution->minT, solution->origMinT, sizeof(int)*nJobs );
    int *sequence = Sol_sequence(solution);
    ModeSet *ms = Sol_getModeSet(solution);

    for ( int i=0; (i<nJobs); i++ ) {

        const Job* job = Inst_job( solution->inst, sequence[i] );

        if(Job_nModes(job) >= 1) {
            int idMode = Modes_job(ms,Job_index(job));
            const Mode* mode = Job_mode( job, idMode );
            int modeDuration = Mode_duration(mode);
            int start = solution->minT[Job_index( job )];
            int end = start + modeDuration;

VALIDATE_STARTING_TIME:
            for ( int rr = 0; (rr<Mode_nResR(mode)); ++rr ) {
                const int useOfR = Mode_useResR(mode, rr);
                const int idR = Mode_idxResR(mode,rr);
                const int capResR = Inst_capResR( solution->inst, idR );
                for ( int tj=end-1; tj>=start; tj-- ) {
                    if ( RRU_usage(solution->rru,tj,idR) + useOfR > capResR ) {
                        start = tj+1;
                        end = start + modeDuration;
                        goto VALIDATE_STARTING_TIME;
                    }
                }
            }

            RRU_allocate(solution->rru,mode,start);
            Sol_setStartJob(solution, Job_index(job), start);

            for( int succ =0; (succ< Job_nSucc(job)); succ++ )
                solution->minT[Job_succ(job,succ)]= MAX(solution->minT[Job_succ(job,succ)], end );
        }
    }
    Sol_fillSequence(solution);
    Sol_calcCost(solution);

}

void Sol_clearBuildSequence(Solution *sol)
{
    assert(sol != NULL);

    Sol_setCost(sol,0);
    nh_reset(sol->priorities);
    RRU_clear(sol->rru);

    memcpy( sol->minT, sol->origMinT, sizeof(int)*Inst_nJobs( sol->inst ) );
}

RRUsage* Sol_getRRUsage( Solution *sol)
{
    return sol->rru;
}

void Sol_setMinTimeSucc(Solution *sol, const Job* job, int time)
{

    //for( int succ =0; (succ< Job_nSucc(job)); succ++ )
    //  sol->minT[Job_succ(job,succ)]= MAX(sol->minT[Job_succ(job,succ)], time );

    for( int succInd =0; (succInd< Inst_nJobs(sol->inst)); succInd++ ) {
        if(Job_hasIndSucc(job,succInd)) {
            sol->minT[succInd]= MAX(sol->minT[succInd], time );
            // printf(" EST %d, sol->minT[%d]: %d Time %d \n", Job_est(job), succInd, sol->minT[succInd], time );
        }
    }

}

int Sol_getMinTime2(Solution *sol, int idxJob)
{
    return sol->minT[idxJob];
}

void Sol_rebuild_opt( Solution *current, const Solution *solOLD)
{

    assert(current != NULL);
    assert(solOLD != NULL);
    int nJobs = Inst_nJobs( current->inst );

    int *sequence = Sol_sequence(current);
    ModeSet *ms = Sol_getModeSet(current);
    int *modes = Modes_getModes(ms);
    int *starts = Sol_startJobs(current);


    const int *seqOLD= Sol_sequence(solOLD);
    const ModeSet *msOLD = Sol_getModeSet(solOLD);
    const int *modeOLD = Modes_getModes(msOLD);
    const int *startOLD =   Sol_startJobs(solOLD);
    unsigned int hash = 0;

    Sol_setCost(current,0);
    memcpy( current->minT, current->origMinT, sizeof(int)*Inst_nJobs( current->inst ) );

    int j = 0;
    for(j = 0; (j<nJobs); j++) {
        const Job* job = Inst_job( current->inst, sequence[j] );
        if(seqOLD[j] != sequence[j] || modeOLD[Job_index(job)] != modes[Job_index(job)])
            break;
        current->minT[Job_index(job)] = solOLD->minT[Job_index(job)];
    }
    if(j < nJobs*0.25) {
        RRU_clear(current->rru);
        RRU_copy_part( current->rru, j, seqOLD, startOLD, modeOLD);
    } else if(j < nJobs)
        RRU_clear_opt(current->rru, j, seqOLD, startOLD, modeOLD, nJobs);

    for(int i = 0; i < j; i++) {
        const Job* job = Inst_job( current->inst, sequence[i] );
        Sol_setPosJob(current, Job_index(job), i);
        int idMode = modes[Job_index(job)];
        const Mode* mode = Job_mode( job, idMode );
        for( int succ =0; (succ< Job_nSucc(job)); succ++ )
            current->minT[Job_succ(job,succ)] = MAX(current->minT[Job_succ(job,succ)], starts[Job_index( job )]+Mode_duration(mode) );

        hash = ((hash << 5) - hash) + starts[Job_index( job )]+Mode_duration(mode);
    }

    for ( int i=j; (i<nJobs); i++ ) {

        const Job* job = Inst_job( current->inst, sequence[i] );
        Sol_setPosJob(current, Job_index(job), i);
        int idMode = modes[Job_index(job)];
        const Mode* mode = Job_mode( job, idMode );
        int start = current->minT[Job_index( job )];
        int end = start + Mode_duration(mode);
VALIDATE_STARTING_TIME:
        for ( int rr = 0; (rr<Mode_nResR(mode)); ++rr ) {
            const int useOfR = Mode_useResR(mode, rr);
            const int idR = Mode_idxResR(mode,rr);

            for ( int tj=end-1; tj>=start; tj-- ) {
                if ( RRU_usage(current->rru,tj,idR) + useOfR > Inst_capResR( current->inst, idR )) {
                    start = tj+1;
                    end = start + Mode_duration(mode);
                    goto VALIDATE_STARTING_TIME;
                }
            }
        }

        hash = ((hash << 5) - hash) + end;

        RRU_allocate(current->rru,mode,start);
        Sol_setStartJob(current, Job_index(job), start);
        for( int succ =0; (succ< Job_nSucc(job)); succ++ )
            current->minT[Job_succ(job,succ)]= MAX(current->minT[Job_succ(job,succ)], end );


    }

    current->solutionHash = hash;

    Sol_calcCost(current);

}

int Sol_getMinTime(Solution *sol, const Job* job, const Mode* mode, int minTime)
{

    int modeDuration = Mode_duration(mode);
    int start = minTime;
    int end = start + modeDuration;

VALIDATE_STARTING_TIME:
    for ( int rr = 0; (rr<Mode_nResR(mode)); ++rr ) {
        const int useOfR = Mode_useResR(mode, rr);
        const int idR = Mode_idxResR(mode,rr);
        const int capResR = Inst_capResR( sol->inst, idR );
        for ( int tj=end-1; tj>=start; tj-- ) {
            //printf("\n tj %d start %d idR %d: %d + %d (%d) > %d\n", tj, start, idR, RRU_usage(sol->rru,tj,idR) , useOfR, ( RRU_usage(sol->rru,tj,idR) + useOfR), capResR);
            if ( RRU_usage(sol->rru,tj,idR) + useOfR > capResR ) {
                start = tj+1;
                end = start + modeDuration;
                goto VALIDATE_STARTING_TIME;
            }
        }
    }

    for( int succ =0; (succ< Job_nSucc(job)); succ++ )
        sol->minT[Job_succ(job,succ)]= MAX(sol->minT[Job_succ(job,succ)], end );


    return start;
}

ModeSet* Sol_getModeSet( const Solution *sol )
{
    assert( sol != NULL);

    return sol->modeSet;
}

const Instance *Sol_inst( const Solution *sol )
{
    assert( sol!=NULL );

    return sol->inst;
}

int *Sol_sequence( const Solution *sol )
{
    assert( sol!= NULL );

    return sol->sequence;
}

int Sol_getPosJob(Solution *sol, int idxJob)
{

    assert( sol != NULL );
    assert( idxJob >= 0 );
    assert( idxJob < Inst_nJobs(sol->inst) );

    return sol->posJobs[idxJob];
}

Cost Sol_getCost( const Solution *sol )
{
    assert( sol != NULL);

    return sol->cost;
}

Cost Sol_getTPD( const Solution *sol )
{
    assert( sol != NULL);

    return sol->TPD;
}

Cost Sol_getTMS( const Solution *sol )
{
    assert( sol != NULL);

    return sol->TMS;
}

void Sol_setCost( Solution *sol, Cost cost )
{
    assert( sol != NULL);

    sol->cost = cost;

}

unsigned int Sol_getSolutionHash( const Solution *sol)
{
    assert( sol != NULL );

    return sol->solutionHash;
}

int Sol_getStartTime( const Solution *sol, int idxJob)
{
    assert( sol != NULL );
    assert( idxJob >= 0 );
    assert( idxJob < Inst_nJobs(sol->inst) );

    return sol->startJobs[idxJob];
}

int Sol_getSequence( const Solution *sol, int idx)
{

    assert( sol != NULL );
    assert( idx >= 0 );
    assert( idx < Inst_nJobs(sol->inst) );

    return sol->sequence[idx];
}

int Sol_getMode( const Solution *sol, int job )
{
    assert( sol != NULL );
    assert( job  >= 0 && job < Inst_nJobs(sol->inst) );

    return Modes_job( sol->modeSet, job );
}

int* Sol_startJobs( const Solution *sol)
{
    assert( sol != NULL );

    return sol->startJobs;
}

void Sol_calcCost( Solution *sol)
{
    assert( sol != NULL );

    sol->TPD = 0;
    sol->TMS = 0;
    int nJobProject = 0;

    for (int p = 0; p < Inst_nProjects(sol->inst); p++) {
        const Project *proj = Inst_project(sol->inst, p);
        nJobProject += Project_nJobs(proj);
        const Job *job = Inst_job(sol->inst, nJobProject-1);
        sol->TPD += Sol_getStartTime( sol, Job_index(job)) - Project_releaseDate(proj) - Project_criticalPath(proj);  // MS - CPD
        sol->TMS = MAX(sol->TMS, Sol_getStartTime(sol, Job_index(job)));
    }

    Cost cost = sol->TPD * 100000 + sol->TMS;
    Sol_setCost(sol, cost);

}

void Sol_cpy( Solution *target, const Solution *sol )
{

    assert( target != NULL );
    assert( sol != NULL );

    target->cost = sol->cost;
    target->TMS = sol->TMS;
    target->TPD = sol->TPD;
    target->solutionHash = sol->solutionHash;

    int nJobs = Inst_nJobs(sol->inst);

    COPY_VECTOR( target->sequence, sol->sequence, int, nJobs );
    COPY_VECTOR( target->startJobs, sol->startJobs, int, nJobs );
    COPY_VECTOR( target->posJobs, sol->posJobs, int, nJobs );
    COPY_VECTOR( target->origMinT, sol->origMinT, int, nJobs );
    COPY_VECTOR( target->minT, sol->minT, int, nJobs );

    Modes_cpy(target->modeSet, sol->modeSet);
    RRU_cpy(target->rru, sol->rru);


}

void Sol_write( const Solution *sol, char *file )
{
    assert( sol != NULL );

    FILE *fp = fopen( file, "w" );
    if ( fp == NULL ) {
        printf( "File was not opened. : path: %s\n", file);
        exit( 0 );
    }

    int nProj = Inst_nProjects(sol->inst);
    int idxOnInstance = 0;
    for( int p = 0; p < nProj; ++p) {
        const Project *proj = Inst_project( sol->inst, p );
        for( int j = 0; j < Project_nJobs( proj ); ++j, ++idxOnInstance )
            fprintf( fp, "%d %d %d %d\n", p, j, Modes_job( sol->modeSet, idxOnInstance ), sol->startJobs[idxOnInstance] );
    }

    fclose( fp );
}

void Sol_print( const Solution *sol )
{
    assert( sol != NULL );

    int nProj = Inst_nProjects( sol->inst );
    int idxOnInstance = 0;
    for( int p = 0; p < nProj; ++p ) {
        const Project *proj = Inst_project( sol->inst,p );
        for( int j = 0; j < Project_nJobs( proj ); j++, idxOnInstance++ )
            printf( "%d %d %d %d\n", p, j, Modes_job( sol->modeSet, idxOnInstance ), sol->startJobs[idxOnInstance] );
    }

    printf( "\nCost Instance NR Proj %ld, inf %d \n", Modes_cost( sol->modeSet ), Modes_inf( sol->modeSet ) );
    printf( "\nCost Instance %ld \n", Sol_getCost( sol ));

}

void Sol_setPosJob(Solution *sol, int idxJob, int value)
{

    assert( sol != NULL );
    assert( idxJob >= 0 );
    assert( idxJob < Inst_nJobs(sol->inst) );

    sol->posJobs[idxJob] = value;
}

void Sol_free( Solution **_solution )
{

    Solution *solution = *_solution;

    Modes_free( &solution->modeSet );
    RRU_free( &solution->rru );
    nh_free(&solution->priorities );

    free( solution->startJobs );
    free( solution->posJobs );
    free( solution->sequence );
    free( solution->minT );
    free( solution->origMinT );

    free( solution );

    *_solution = NULL;
}

/* Function of mutation to GA*/
void Sol_firstMutation( Solution *sol )
{

    //    int * sequence = Sol_sequence(sol);

    //int idx = RAND; Inst_nJobs(Sol_inst(sol)));
    //int newpos = rand(Inst_nJobs(Sol_inst(sol)));

    // function of neighborhood to change one job on sequence
    // int_
    // printf("Job %d, Old mode %d New mode %d\n", idx, Sol_getMode(sol,idx), newmode);
    // int factivel = Modes_modifyAndVerify(modeset, idx, newmode));


}

/* Function of crossover to GA*/
void Sol_firstCrossover( Solution *sol )
{

}

Cost Sol_evaluate( Solution *sol)
{
    assert( sol != NULL );

    sol->TPD = 0;
    sol->TMS = 0;
    int nJobProject = 0;

    for (int p = 0; p < Inst_nProjects(sol->inst); p++) {
        const Project *proj = Inst_project(sol->inst, p);
        nJobProject += Project_nJobs(proj);
        const Job *job = Inst_job(sol->inst, nJobProject-1);
        sol->TPD += Sol_getStartTime( sol, Job_index(job)) - Project_releaseDate(proj) - Project_criticalPath(proj);  // MS - CPD
        sol->TMS = MAX(sol->TMS, Sol_getStartTime(sol, Job_index(job)));
    }

    Cost cost = sol->TPD * 100000 + sol->TMS;
    Sol_setCost(sol, cost);

    return cost;
}

/*int *Sol_randomSequence( const Solution *sol)
{

    assert( sol != NULL );

    int *sequence;
    ALLOCATE_VECTOR_INI( sequence, int, Inst_nJobs(sol->inst) );

    for(int i = 0; i < Inst_nJobs(sol->inst); ++i)
        sequence[i] = i;

    for(int s = 0; s < Inst_nJobs(sol->inst); ++s) {

        int r = rand() % Inst_nJobs(sol->inst);

        int temp = sequence[s];
        sequence[s] = sequence[r];
        sequence[r] = temp;
    }

    return sequence;
}
*/













/*
void Sol_clear( Solution *sol )
{
    memset(sol->startJobs, 0, sizeof(int)*Inst_nJobs( sol->inst ));
    memset(sol->posJobs, 0, sizeof(int)*Inst_nJobs( sol->inst ));
    Sol_setCost(sol,0);
    nh_reset(sol->priorities);
    RRU_clear(sol->rru);
    Modes_isntEmpty(sol->modeSet);

    for (int j = 0; j<(Inst_nJobs( sol->inst ) ); j++ ) {
        const Job* job = Inst_job(sol->inst, j);
        sol->origMinT[j] = Job_est(job);
    }

    memcpy( sol->minT, sol->origMinT, sizeof(int)*Inst_nJobs( sol->inst ) );

}

void Sol_setModesByProj( Solution *sol, const ModeSet *modeSet )
{
    assert( modeSet != NULL );
    assert( sol != NULL );

    Modes_cpyToSolByProj( sol->modeSet, modeSet );
}


void Sol_setMode( Solution *sol, int job, int mode )
{
    assert( sol != NULL );

    Modes_modifyCount( sol->modeSet, job, mode );

}


RRUsage *Sol_rru( const Solution *sol )
{

    assert( sol != NULL);

    return sol->rru;

}


void Sol_setSequence(Solution *sol, int idx, int idJob)
{

    assert( sol != NULL );
    assert( idx >= 0 );
    assert( idx < Inst_nJobs(sol->inst) );
    assert( idJob >= 0 );
    assert( idJob < Inst_nJobs(sol->inst) );

    sol->sequence[idx] = idJob;

}



void Sol_reconstruct( Solution *solution, int sequenceOLD[], int startJobsOLD[], int modesOLD[] )
{

    assert(solution != NULL);

    //  int *sequence = Sol_sequence(solution);
    // int *startJobs = Sol_startJobs(solution);
    ModeSet *ms = Sol_getModeSet(solution);
    //    int *modes = Modes_getModes(ms);
    Sol_setCost(solution,0);
    Modes_clear(solution->modeSet);
    RRU_clear(solution->rru);
    //  int j = 0;
    //  for(j = 0 ; (j<Inst_nJobs(solution->inst)) ; j++) {
    //      if(sequenceOLD[j] != sequence[j] || modesOLD[sequenceOLD[j]] != modes[sequence[j]]) {
    //          RRU_clear_opt(solution->rru, j, sequence, startJobs, modes, Inst_nJobs(solution->inst));
    //          break;
    //      }
    //  }

    for ( int i=0; (i<Inst_nJobs( solution->inst )); i++ ) {

        const Job* job = Inst_job( solution->inst, sequenceOLD[i] );
        Sol_setSequence(solution, i, Job_index(job));
        Sol_setPosJob(solution, Job_index(job), i);

        int start = startJobsOLD[Job_index(job)];
        Sol_setStartJob(solution, Job_index(job), start);

        int idMode = modesOLD[Job_index(job)];
        const Mode*mode = Job_mode(job,idMode);
        Modes_modify( ms, Job_index(job), idMode);
        RRU_allocate(solution->rru,mode,start);

    }
    Modes_isntEmpty(ms);
    Sol_calcCost(solution);

}

void Sol_build( Solution *solution )
{

    int ret, ind, costInd;

    for ( int i=0; (i<Inst_nJobs( solution->inst )); i++ ) {
        const Job* job = Inst_job( solution->inst, i );
        nh_update(solution->priorities, Job_index(job), i+Job_nPred(job) * Inst_nJobs(solution->inst));
    }

    for ( int i=0; (i<Inst_nJobs( solution->inst )); i++ ) {

        nh_remove_first(solution->priorities, &ret);
        const Job* job = Inst_job(solution->inst, ret);
        for(int h = 0; h < Job_nSucc(job); ++h) {
            const Job* job2 =  Inst_job(solution->inst, Job_succ(job,h));
            ind = Job_index(job2);
            costInd = nh_get_dist(solution->priorities, ind);
            nh_update(solution->priorities, ind, costInd - Inst_nJobs(solution->inst));
        }

        int idMode = Modes_job(solution->modeSet, Job_index(job));

        const Mode* mode = Job_mode( job, idMode );
        int start = solution->minT[Job_index( job )];
        int end = start + Mode_duration(mode);

        // resource check and updating t
VALIDATE_STARTING_TIME:

        for ( int rr = 0; (rr<Mode_nResR(mode)); ++rr ) {
            const int useOfR = Mode_useResR(mode, rr);
            const int idR = Mode_idxResR(mode,rr);
            for ( int tj=end-1; tj>=start; tj-- ) {
                if ( RRU_usage(solution->rru,tj,idR) + useOfR > Inst_capResR( solution->inst, idR )) {
                    start = tj+1;
                    end = start + Mode_duration(mode);
                    goto VALIDATE_STARTING_TIME;
                }
            }
        }



        // end of resource check

        RRU_allocate(solution->rru,mode,start);
        Sol_setStartJob(solution, Job_index(job), start);

        for( int succ =0; (succ< Job_nSucc(job)); succ++ )
            solution->minT[Job_succ(job,succ)]= MAX(solution->minT[Job_succ(job,succ)], end );

    }
    Sol_fillSequence(solution);
    Sol_calcCost(solution);

}


void Sol_buildBySequence( Solution *solution, int sequence[],...)
{

    assert( solution != NULL);

    ModeSet *ms = Sol_getModeSet(solution);

    va_list arg;
    va_start(arg, sequence);
    int validTopSort = va_arg(arg,int);
    va_end(arg);

    if(validTopSort != 1)
        Sol_topSort(solution, sequence);

    for ( int i=0; (i<Inst_nJobs( solution->inst )); i++ ) {

        // job which will be allocated
        const Job* job = Inst_job( solution->inst, sequence[i] );
        int idMode = Modes_job(ms,Job_index(job));
        const Mode* mode = Job_mode( job, idMode );
        int start = solution->minT[Job_index( job )];
        int end = start + Mode_duration(mode);

        // resource check and updating t
VALIDATE_STARTING_TIME:

        for ( int rr = 0; (rr<Mode_nResR(mode)); ++rr ) {
            const int useOfR = Mode_useResR(mode, rr);
            const int idR = Mode_idxResR(mode,rr);
            for ( int tj=end-1; tj>=start; tj-- ) {
                if ( RRU_usage(solution->rru,tj,idR) + useOfR > Inst_capResR( solution->inst, idR )) {
                    start = tj+1;
                    end = start + Mode_duration(mode);
                    goto VALIDATE_STARTING_TIME;
                }
            }
        }
        // end of resource check

        RRU_allocate(solution->rru,mode,start);
        Sol_setStartJob(solution, Job_index(job), start);

        for( int succ =0; (succ< Job_nSucc(job)); succ++ )
            solution->minT[Job_succ(job,succ)]= MAX(solution->minT[Job_succ(job,succ)], end );

    }
    Sol_fillSequence(solution);
    Sol_calcCost(solution);

}


void Sol_changeSequence( Solution *sol, int pos, int job, int sequence[] )
{
    assert( sol != NULL);
    assert( pos > 0);
    assert( pos < Inst_nJobs( sol->inst ) );
    assert( job > 0);
    assert( job < Inst_nJobs( sol->inst ) );

    sequence[pos] = job;
    Sol_setPosJob(sol,pos,job);
}


int* Sol_getModes( const Solution *sol )
{
    assert( sol != NULL);

    ModeSet *ms = sol->modeSet;

    return Modes_getModes(ms);
}


int Sol_getJobOnSequence( const Solution *sol, int idx)
{

    assert( sol != NULL );
    assert( idx >= 0 );
    assert( idx < Inst_nJobs(sol->inst) );

    return sol->sequence[idx];
}

void Sol_cpySequence(const Solution *sol, int sequenceOLD[])
{
    assert( sol != NULL);
    COPY_VECTOR( sequenceOLD, sol->sequence, int, Inst_nJobs( sol->inst ) );
}

void Sol_cpyStartJobs(const Solution *sol, int startJobsOLD[])
{

    assert( sol != NULL);

    COPY_VECTOR( startJobsOLD, sol->startJobs, int, Inst_nJobs( sol->inst ) );

}
*/


/*
void Sol_rebuild_opt_old( Solution *current, const Solution *solOLD)
{

    assert(current != NULL);
    assert(solOLD != NULL);

    int nJobs = Inst_nJobs( current->inst );

    int *sequence = Sol_sequence(current);
    ModeSet *ms = Sol_getModeSet(current);
    int *modes = Modes_getModes(ms);
    int *starts = Sol_startJobs(current);


    const int *seqOLD= Sol_sequence(solOLD);
    const ModeSet *msOLD = Sol_getModeSet(solOLD);
    const int *modeOLD = Modes_getModes(msOLD);
    const int *startOLD =   Sol_startJobs(solOLD);

    Sol_setCost(current,0);

    memcpy( current->minT, current->origMinT, sizeof(int)*Inst_nJobs( current->inst ) );


    int j = 0;


    for(j = 0; (j<nJobs); j++) {
        const Job* job = Inst_job( current->inst, sequence[j] );
        if(seqOLD[j] != sequence[j] || modeOLD[Job_index(job)] != modes[Job_index(job)]) {
            RRU_clear_opt(current->rru, j, seqOLD, startOLD, modeOLD, nJobs);
            break;
        }
        current->minT[Job_index(job)] = solOLD->minT[Job_index(job)];

    }

    for(int i = 0; i < j; i++) {
        const Job* job = Inst_job( current->inst, sequence[i] );
        Sol_setPosJob(current, Job_index(job), i);
        int idMode = modes[Job_index(job)];
        const Mode* mode = Job_mode( job, idMode );
        for( int succ =0; (succ< Job_nSucc(job)); succ++ )
            current->minT[Job_succ(job,succ)] = MAX(current->minT[Job_succ(job,succ)], starts[Job_index( job )]+Mode_duration(mode) );
    }

    for ( int i=j; (i<nJobs); i++ ) {

        const Job* job = Inst_job( current->inst, sequence[i] );
        Sol_setPosJob(current, Job_index(job), i);
        int idMode = modes[Job_index(job)];
        const Mode* mode = Job_mode( job, idMode );
        int start = current->minT[Job_index( job )];
        int end = start + Mode_duration(mode);

VALIDATE_STARTING_TIME:
        for ( int rr = 0; (rr<Mode_nResR(mode)); ++rr ) {
            const int useOfR = Mode_useResR(mode, rr);
            const int idR = Mode_idxResR(mode,rr);

            for ( int tj=end-1; tj>=start; tj-- ) {
                if ( RRU_usage(current->rru,tj,idR) + useOfR > Inst_capResR( current->inst, idR )) {
                    start = tj+1;
                    end = start + Mode_duration(mode);
                    goto VALIDATE_STARTING_TIME;
                }
            }
        }

        RRU_allocate(current->rru,mode,start);
        Sol_setStartJob(current, Job_index(job), start);
        for( int succ =0; (succ< Job_nSucc(job)); succ++ )
            current->minT[Job_succ(job,succ)]= MAX(current->minT[Job_succ(job,succ)], end );

    }

    //Sol_fillSequence(current);
    Sol_calcCost(current);

}
*/
