/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <time.h>
#include <assert.h>
#include "ms_solver_lahc.h"
#include "mode_set.h"
#include "neighborhood.h"
#include "macros.h"

struct _MSSolverLAHC {
    ModeSet *currentModes;
    ModeSet *bestModes;
    Neighborhood *neighborhood;

    int lfa;
    Cost *f;

    const struct _Instance *inst;
};
/*
MSSolverLAHC *MS_createByProj( const Instance *inst, int sizeList, int firstJob, int lastJob )
{
    MSSolverLAHC* msLahc;

    ALLOCATE_INI( msLahc, MSSolverLAHC );

    msLahc->inst = inst;

    msLahc->neighborhood = Neighbor_createMS(inst);

    msLahc->lfa = sizeList;

    int nResN = Inst_nResN(inst)/Inst_nProjects(inst);
    msLahc->currentModes = Modes_createForMSP(inst, firstJob, lastJob, nResN);
    msLahc->bestModes = Modes_createForMSP(inst, firstJob, lastJob, nResN);

    ALLOCATE_VECTOR(msLahc->f, Cost, msLahc->lfa);

    return msLahc;
}*/
/*
MSSolverLAHC *MS_create( const Instance *inst, int sizeList)
{
    MSSolverLAHC* msLahc;

    ALLOCATE_INI( msLahc, MSSolverLAHC );

    msLahc->inst = inst;

    msLahc->neighborhood = Neighbor_createMS(inst);

    msLahc->lfa = sizeList;

    msLahc->currentModes = Modes_createForMS(inst);
    msLahc->bestModes = Modes_createForMS(inst);

    ALLOCATE_VECTOR(msLahc->f, Cost, msLahc->lfa);

    return msLahc;
}

*/
void MS_free( MSSolverLAHC **_msLahc )
{
    MSSolverLAHC *msLahc = *_msLahc;

    Modes_free( &msLahc->currentModes);
    Modes_free( &msLahc->bestModes);
    Neighbor_free( &msLahc->neighborhood );

    free( msLahc->f );
    free( msLahc );

    *_msLahc = NULL;
}

void MS_run(MSSolverLAHC *msLahc, int iterat)
{

    /*filling the cost vector with the value of the current solution*/

    for(int i = 0; i< (int) msLahc->lfa; ++i)
        msLahc->f[i] = Modes_cost(msLahc->currentModes);

    int v, I=0, iteratRem=0;
    Cost oldCost;
    int *modes;
    int repeat=1;

    ALLOCATE_VECTOR(modes, int, Modes_nJobs(msLahc->currentModes));

    Modes_fillModes(msLahc->currentModes, modes);

    while( repeat==1 || I < iteratRem) {
        oldCost = Modes_cost(msLahc->currentModes);

        Neighbor_callStocMS( msLahc->neighborhood, msLahc->currentModes);

        v = I % msLahc->lfa;
        if( Modes_inf(msLahc->bestModes) == 0 && repeat==1 ) {
            repeat = 0;
            iteratRem = iterat + I;
        }

        if( Modes_cost(msLahc->currentModes) < msLahc->f[v] || Modes_cost(msLahc->currentModes) < oldCost ) {
            if( Modes_cost(msLahc->currentModes) < Modes_cost(msLahc->bestModes)) {
                //printf("\nImproved from %ld, to %ld\n", Modes_cost(msLahc->bestModes), Modes_cost(msLahc->currentModes));
                Modes_clear(msLahc->bestModes);
                Modes_cpy(msLahc->bestModes, msLahc->currentModes);
                if ( repeat==0 )
                    iteratRem = iterat + I;
            }
            Modes_fillModes(msLahc->currentModes, modes);
        } else
            Modes_reconstructsModes( msLahc->currentModes, modes );

        //   printf("passou\n");
        msLahc->f[v] = Modes_cost(msLahc->currentModes);
        I++;

    }

    free( modes );

    // printf("\n\n Best cost %ld, inf %d\n", Modes_cost(msLahc->bestModes), Modes_inf(msLahc->bestModes) );

}

void MS_runByProj(const Instance *inst, Solution *sol, int LFA, int iterat)
{
    int firstJob = 0, lastJob =0;
    for(int p = 0; p <Inst_nProjects(inst); ++p) {
        const Project* proj = Inst_project(inst,p);
        lastJob = firstJob + Project_nJobs(proj)-1;
        MSSolverLAHC *msLahc = MS_createByProj(inst,LFA, firstJob, lastJob);
        MS_run(msLahc, iterat);
        firstJob = lastJob+1;
        Sol_setModesByProj(sol, MS_bestModes(msLahc));

        MS_free( &msLahc );
    }
}

ModeSet *MS_bestModes( const MSSolverLAHC *solver )
{
    assert(solver != NULL);

    return solver->bestModes;
}

Neighborhood *MS_getNeighborhood( const MSSolverLAHC *solver )
{
    assert(solver != NULL);

    return solver->neighborhood;
}
