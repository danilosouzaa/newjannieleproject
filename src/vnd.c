/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <time.h>
#include <assert.h>
#include <string.h>
#include <omp.h>
#include "test.h"
#include "vnd.h"
#include "neighborhood.h"

struct _VND {

    Solution *iniSol;
    Solution *bestSol;

    int nThread;

    int parallel;

    const struct _Instance *inst;

};

void VND_checkArgs(VND *vnd, char **argv, int argc)
{

    for(int n=0 ; n < argc; n++) {

        if (strcmp(argv[n],"-parallel") == 0) {
            n++;
            vnd->parallel = atoi(argv[n]);
            continue;
        }

        if (strcmp(argv[n],"-nThread") == 0) {
            n++;
            vnd->nThread = atoi( argv[n] );
            continue;
        }


    }
}

VND *VND_create( const Instance *inst, Solution *sol, char **argv, int argc)
{

    assert(inst != NULL);
    assert(sol != NULL);

    VND* vnd;

    ALLOCATE_INI( vnd, VND );

    vnd->iniSol = sol;
    vnd->bestSol = Sol_create(inst);
    Sol_cpy( vnd->bestSol, vnd->iniSol );
    vnd->parallel = 0;
    vnd->nThread = 1;

    VND_checkArgs(vnd, argv, argc);

    vnd->inst = inst;

    return vnd;
}

int VND_getParallel(VND *vnd)
{
    assert(vnd != NULL);
    return vnd->parallel;

}

void VND_runDet(VND *vnd, Neighborhood *neighborhood, double timeRem, int nNeighbor, int firstImprovement, Test *test)
{

    assert(neighborhood != NULL);
    assert(vnd != NULL);

    double _time = 0;
    clock_t tStart = clock();

    int kN = 0;
    int improved = 0;

    Solution *current = Sol_create(vnd->inst);
    Solution **currentT;
    ALLOCATE_VECTOR(currentT, Solution*, vnd->nThread);

    for (int i = 0; i < vnd->nThread; i++)
        currentT[i]= Sol_create(Sol_inst( vnd->iniSol));

    while(kN < nNeighbor) {

        _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

        if(!VND_getParallel(vnd)) {
            Sol_cpy(current, vnd->iniSol);
            improved = Neighbor_callDetLS( neighborhood,
                                           vnd->iniSol,
                                           vnd->bestSol,
                                           current,
                                           //nChangesModes,
                                           //nTimesJobOnModes,
                                           kN,
                                           firstImprovement,
                                           timeRem - _time,
                                           test);
            Neighbor_setNullLastJobModify(neighborhood);
        } else {
            for (int i = 0; i < vnd->nThread; i++)
                Sol_cpy(currentT[i], vnd->iniSol);

            improved = Neighbor_callDetLS_Parallel( neighborhood,
                                                    vnd->iniSol,
                                                    vnd->bestSol,
                                                    currentT,
                                                    // nChanges,
                                                    kN,
                                                    firstImprovement,
                                                    timeRem - _time,
                                                    vnd->nThread,
                                                    test);
            Neighbor_setNullLastJobModify(neighborhood);
        }

        if(improved && kN) kN = 0;
        else kN++;
    }

    Sol_free(&current);

    for(int i = 0 ; i < vnd->nThread; i++)
        Sol_free(&currentT[i]);

    free(currentT);

}

Solution *VND_getBestSol(VND *vnd)
{
    assert( vnd != NULL );

    return vnd->bestSol;
}

void VND_free( VND **_vnd )
{

    VND *vnd = *_vnd;
    Sol_free( &vnd->iniSol );
    Sol_free( &vnd->bestSol );

    free( vnd );

    *_vnd = NULL;

}
