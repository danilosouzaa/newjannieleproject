
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "top_sort.h"
#include "macros.h"

struct _TopSort {
    const Instance *inst;

    int *jobsTS;
    int nJobsTS;

    int *degree;
};

TopSort *TopSort_create( const Instance *inst )
{
    TopSort *ts;
    ALLOCATE_INI( ts, TopSort );
    ts->inst = inst;

    ALLOCATE_VECTOR( ts->jobsTS, Inst_nJobs(inst), Inst_nJobs( inst ) );
    ALLOCATE_VECTOR_INI( ts->degree, Inst_nJobs(inst), Inst_nJobs( inst ) );

    return ts;
}

void TopSort_computeProjectTS( TopSort *ts, int idxProject )
{
    const Instance *inst = ts->inst;

    int *degree = ts->degree;
    FILL( degree, 0, Inst_nJobs(inst), 0 );
    int *jobsTS = ts->jobsTS;

    const Project *project = Inst_project( inst, idxProject );
    const int nJobsProj = Project_nJobs( project );
    const int idxFirstJobP = Project_idxFirstJob( project );

    int stk[nJobsProj], nStk=0;

    stk[nStk++] = idxFirstJobP;

    for ( int j=idxFirstJobP ; (j<idxFirstJobP+nJobsProj) ; ++j ) {
        const Job *job = Inst_job( inst, j );
        degree[ j ] += Job_nPred( job );
    }

    ts->nJobsTS = 0;
    while ( nStk ) {
        int nextJob = stk[--nStk];
        jobsTS[ts->nJobsTS++] = nextJob;
        const Job *jNext = Inst_job( inst, nextJob);

        for ( int i=0 ; (i<Job_nSucc(jNext)) ; ++i ) {
            const int idxSucc = Job_succ( jNext, i );
            degree[ idxSucc ]--;
            if (!degree[ idxSucc ])
                stk[nStk++] = idxSucc;
        }
    }

    /* checking */
    assert( ts->nJobsTS == Project_nJobs( project ) );
}

int TopSort_nTS( const TopSort *ts )
{
    return ts->nJobsTS;
}

int *TopSort_ts( const TopSort *ts )
{
    return ts->jobsTS;
}

void TopSort_free( TopSort **_topSort )
{
    TopSort *ts = *_topSort;

    free( ts->degree );
    free( ts->jobsTS );
    free( ts );
    *_topSort = NULL;

}



