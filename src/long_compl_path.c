#include <stdlib.h>
#include <stdio.h>
#include "long_compl_path.h"
#include "macros.h"
#include "top_sort.h"

struct _LongestComplPath {
    const Instance *inst;
    TopSort *topSort;
    int *lp;

    int *tlp;
};

LongestComplPath *LongCP_create( const Instance *inst )
{
    LongestComplPath *lcp;
    ALLOCATE_INI( lcp, LongestComplPath );
    lcp->inst = inst;
    lcp->topSort = TopSort_create( inst );
    int nJobs = Inst_nJobs(inst);
    ALLOCATE_VECTOR_INI( lcp->lp, int, nJobs);
    ALLOCATE_VECTOR_INI( lcp->tlp, int, nJobs);

    return lcp;
}

void LongCP_solve( LongestComplPath *lcp )
{
    const Instance *inst = lcp->inst;

    int *tlp = lcp->tlp;
    int *lp = lcp->lp;

    const int nProjects = Inst_nProjects( inst );
    for ( int ip=0 ; (ip<nProjects) ; ip++ ) {
        TopSort_computeProjectTS( lcp->topSort, ip );
        const int idxFirstJobProject = Project_idxFirstJob( Inst_project( inst, ip ) );
        const int nJobsProject = Project_nJobs( Inst_project( inst, ip ) );
        const int idxLastJobProject = idxFirstJobProject+nJobsProject-1;

        const int *ts = TopSort_ts( lcp->topSort );
        const int nTS = TopSort_nTS( lcp->topSort );

        for ( int j1=0 ; j1<nTS ; ++j1 ) {
            FILL( tlp, idxFirstJobProject, idxFirstJobProject+nTS, 0 );

            for ( int ij=j1+1 ; (ij<nTS) ; ++ij ) {
                const int j = ts[ij];
                const Job *job = Inst_job( inst, j );
                for ( int ipr=0 ; (ipr<Job_nPred(job)) ; ++ipr  ) {
                    const int idxPred = Job_pred( job, ipr );
                    const Job *jPred = Inst_job( inst, idxPred );
                    const int minDPred = Job_minDuration( jPred );
                    tlp[j] = MAX( tlp[j], tlp[idxPred]+minDPred );
                } /* all predecessors of this job */
            } /* all jobs of this project */

            lp[ts[j1]] = tlp[idxLastJobProject];

        }
    } /* all projects */
}


const int *LongCP_get( const LongestComplPath *lcp )
{
    return lcp->lp;
}

void LongCP_free( LongestComplPath **_lcp )
{
    LongestComplPath *lcp = *_lcp;

    TopSort_free( &lcp->topSort );
    free( lcp->tlp );
    free( lcp->lp );
    free( lcp );

    *_lcp = NULL;
}
