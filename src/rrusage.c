
/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 *
 * rrusage: controls the usage of renewable resources in a solution
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rrusage.h"
#include "assert.h"
#include "vec_int.h"
#include "macros.h"
#include "memory.h"

#define INI_CAP 2048

struct _RRUsage {

    /* number of renewable resources */
    int nResR;

    /* already allocated timeslots */
    int tsCap;
    int tsUsed;

    /* use resor*/
    int **useR;
    const struct _Instance *inst;
};

RRUsage *RRU_create( const Instance *inst )
{
    assert( inst!=NULL );

    RRUsage* rru;

    ALLOCATE_INI( rru, RRUsage );
    rru->inst = inst;
    rru->tsCap = INI_CAP;
    rru->nResR = Inst_nResR( inst );
    ALLOCATE_VECTOR( rru->useR, int*, rru->nResR );
    for ( int i=0; (i<rru->nResR); ++i )
        ALLOCATE_VECTOR_INI( rru->useR[i], int, rru->tsCap );
    rru->tsUsed = 0;

    return rru;
}

int RRU_usage( const RRUsage *rru, int time, int idxRR)
{
    assert( rru!=NULL );
    assert( idxRR >= 0 );
    assert( idxRR < Inst_nResR( rru->inst ) );

    return rru->useR[idxRR][time];

}

void RRU_clear( RRUsage *rru)
{

    assert( rru != NULL);

    for ( int i=0; (i<Inst_nResR(rru->inst)); ++i )
        memset( rru->useR[i], 0, (rru->tsUsed)*sizeof(int) );

    rru->tsUsed = 0;
}

void RRU_clear_opt( RRUsage *rru, int j, const int seq[], const int starts[], const int modes[], int nJobs)
{
    assert( rru != NULL);
    assert( j >= 0 && j <= Inst_nJobs( rru->inst ) );

    for ( int i=j; (i < nJobs); i++ ) {

        const Job* job = Inst_job( rru->inst, seq[i] );
        const Mode* modeOl = Job_mode( job, modes[seq[i]] );

        for ( int r=0; (r<Mode_nResR(modeOl)); ++r ) {
            int use = Mode_useResR(modeOl,r);
            for(int jj = starts[ seq[i]]; jj < starts[seq[i]]+Mode_duration(modeOl); jj++) {
                //printf("\nMode OLD %d, Mode_duration %d, UsageR%d OLD at time %d is %d, Cap %d\n", Mode_index(modeOl), Mode_duration(modeOl), Mode_idxResR(modeOl,r),jj,rru->useR[Mode_idxResR(modeOl,r)][jj],  Inst_capResR( rru->inst, Mode_idxResR(modeOl,r) ));
                rru->useR[Mode_idxResR(modeOl,r)][jj] -= use;
                //printf("\nMode OLD %d,  Mode_duration %d,UsageR%d Novo desalocado at time %d is %d - %d, Cap %d\n", Mode_index(modeOl), Mode_duration(modeOl),  Mode_idxResR(modeOl,r),jj,rru->useR[Mode_idxResR(modeOl,r)][jj], use, Inst_capResR( rru->inst, Mode_idxResR(modeOl,r) ));
                // if(rru->useR[Mode_idxResR(modeOl,r)][jj] < 0 ) {
                //  printf("\nMode OLD %d,  Mode_duration %d,UsageR%d Novo desalocado at time %d is %d - %d, Cap %d\n", Mode_index(modeOl), Mode_duration(modeOl),  Mode_idxResR(modeOl,r),jj,rru->useR[Mode_idxResR(modeOl,r)][jj], use, Inst_capResR( rru->inst, Mode_idxResR(modeOl,r) ));
                //getchar();
                //}

            }
        }

    }

}

void RRU_allocate( RRUsage *rru, const Mode *mode, int start )
{


    assert( rru != NULL);
    assert( mode != NULL);

    if ( rru->tsUsed+Mode_duration(mode)+1 >= rru->tsCap ) {
        const int oldCap = rru->tsCap;

        rru->tsCap = MAX( 2*rru->tsCap, rru->tsUsed+Mode_duration(mode)+1 );

        for ( int i=0; (i<Inst_nResR(rru->inst)); ++i ) {
            int *v = xrealloc( rru->useR[i], sizeof(int)*rru->tsCap );
            rru->useR[i] = v;
            memset( v+oldCap, 0, (rru->tsCap-oldCap)*sizeof(int) );
        }
    }

    for ( int i=0; (i<Mode_nResR(mode)); ++i ) {
        const int idxR = Mode_idxResR(  mode, i  );
        const int use = Mode_useResR( mode, i );

        /* computing resource usage */
        int j, endTime = start+Mode_duration(mode);
        for ( j=start; (j<endTime); ++j ) {
            //            printf("\nMode %d, Mode_duration %d, UsageR%d OLD at time %d is %d, Cap %d\n", Mode_index(mode), Mode_duration(mode), idxR,j,rru->useR[idxR][j],  Inst_capResR( rru->inst, idxR ));
            rru->useR[idxR][j] += use;
            //            printf("\nMode %d, Mode_duration %d, UsageR%d novo deslocadoat time %d is %d + %d, Cap %d\n", Mode_index(mode), Mode_duration(mode), idxR,j,rru->useR[idxR][j], use,  Inst_capResR( rru->inst, idxR ));

            assert( rru->useR[idxR][j] <= Inst_capResR( rru->inst, idxR ) );
        }
    }

    rru->tsUsed = MAX( rru->tsUsed, start+Mode_duration(mode) );
}

void RRU_cpy( RRUsage *target, const RRUsage *rru)
{

    assert( target != NULL);
    assert( rru != NULL);

    if(target->tsUsed <= rru->tsUsed) {
        if(target->tsCap < rru->tsCap) {
            for ( int i=0; (i<Inst_nResR(rru->inst)); ++i ) {
                int *v  = xrealloc( target->useR[i], sizeof(int)*rru->tsCap );
                target->useR[i] = v;
                COPY_VECTOR( target->useR[i], rru->useR[i], int, rru->tsUsed );
            }
        } else {
            for ( int i=0; (i<Inst_nResR(rru->inst)); ++i )
                COPY_VECTOR( target->useR[i], rru->useR[i], int, rru->tsUsed );
        }
    } else {
        for ( int i=0; (i<Inst_nResR(rru->inst)); ++i ) {
            int *v =  target->useR[i];
            int tsUsed = rru->tsUsed;
            CLEAR_VECTOR( v+tsUsed, int, (target->tsUsed-rru->tsUsed) );
            COPY_VECTOR( target->useR[i], rru->useR[i], int, rru->tsUsed );
        }
    }

    target->tsUsed = rru->tsUsed;
    target->tsCap = MAX( rru->tsCap, target->tsCap );
}

void RRU_copy_part( RRUsage *rru, int j, const int seq[], const int starts[], const int modes[])
{
    assert( rru != NULL);
    fflush(stdout);

    assert( j >= 0 && j <= Inst_nJobs( rru->inst ) );

    for ( int i=0; (i < j); i++ ) {

        const Job* job = Inst_job( rru->inst, seq[i] );
        const Mode* modeOl = Job_mode( job, modes[seq[i]] );

        for ( int r=0; (r<Mode_nResR(modeOl)); ++r ) {
            int use = Mode_useResR(modeOl,r);
            for(int jj = starts[ seq[i]]; jj < starts[seq[i]]+Mode_duration(modeOl); jj++) {
                //  printf("\nMode OLD %d, Mode_duration %d, UsageR%d OLD at time %d is %d, Cap %d\n", Mode_index(modeOl), Mode_duration(modeOl), Mode_idxResR(modeOl,r),jj,rru->useR[Mode_idxResR(modeOl,r)][jj],  Inst_capResR( rru->inst, Mode_idxResR(modeOl,r) ));
                rru->useR[Mode_idxResR(modeOl,r)][jj] += use;

            }
        }

    }

}

void RRU_free( RRUsage **_rru )
{
    RRUsage *rru = *_rru;
    for ( int i=0; (i<rru->nResR); ++i )
        free( rru->useR[i] );

    free( rru->useR);
    free( rru );
    *_rru = NULL;
}






















void RRU_deallocate( RRUsage *rru, const Mode *mode, int start )
{
    for ( int i=0; (i<Mode_nResR(mode)); ++i ) {
        const int idxR = Mode_idxResR(  mode, i  );
        const int use = Mode_useResR( mode, i );

        // computing resource usage
        int j, endTime = start+Mode_duration(mode);
        for ( j=start; (j<endTime); ++j ) {
            rru->useR[idxR][j] -= use;
            // printf("\nTime %d, userR %d, Cap %d\n", j, rru->useR[idxR][j], Inst_capResR( rru->inst, idxR ));
            assert( rru->useR[idxR][j] <= Inst_capResR( rru->inst, idxR ) );
        }
    }

}


/*
int RRU_find( const RRUsage *rru, const Mode *mode, int startingTime )
{

    assert( rru!=NULL );
    assert( mode!=NULL );

    // TODO: check if it is optimized

    for ( int i=0; (i<Mode_nResR(mode)); ++i ) {
        const int idxR = Mode_idxResR(  mode, i  );
        const int cap = Inst_capResR( rru->inst, idxR );
        const int use = Mode_useResR( mode, i );
        // first time where this resource is available
        int j;
        for ( j=startingTime; (j<rru->tsCap); ++j )
            if ( rru->useR[idxR][startingTime]+use <= cap ) {
                startingTime = j;
                break;
            }
    }

    return startingTime;
}
*/
