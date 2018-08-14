/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "instance.h"
#include "ms_solver_lahc.h"
#include "vnd.h"
#include "neighborhood.h"
#include "solution.h"
#include "macros.h"
#include "time.h"
#include "test.h"


int main( int argc, char **argv )
{
    clock_t t_begin, t_end;
    double _time = 0;
    t_begin = clock();

    if ( argc< 6 ) {
        fprintf( stderr, "ERROR! Enter path, instance, _nameFile, firstImproved(0,1), time, _nameSolini\n");
        fflush( stderr );
        exit( EXIT_FAILURE );
    }

    int tmp = atoi( argv[5]);

    srand(100000);

    Instance *inst = Inst_create( argv[1], argv[2] );

    Solution *sol = Sol_create(inst);

    Sol_read(sol,  argv[6] );

    Cost initFO = Sol_getCost(sol);

    Test *test = Test_create(14);

    Neighborhood *neigh = Neighbor_create(inst, argv, argc);

    /*Local Search*/
    VND *vnd = VND_create(inst, sol, argv, argc);

    _time = ( (double) tmp - (double)(clock()-t_begin )/CLOCKS_PER_SEC );

    VND_runDet(vnd, neigh, _time, 14, atoi(argv[4]), VND_getNChangesModes(vnd), VND_getNTimesJobOnModes(vnd), test);

    /*Writing solution*/
    Sol_write( sol,   argv[3] );

    t_end = clock();
    _time = ( (double)( t_end-t_begin )/CLOCKS_PER_SEC );

#ifdef DEBUG
    printf( "\nCost Instance: %ld |||| totalTime = %f | Infact = %d\n", Sol_getCost(sol), _time, Modes_inf( Sol_getModeSet( sol ) ) );
    fflush(stdin);
    fflush(stdout);
#endif // DEBUG

    Test_setTotalTime(test, _time);
    Test_writeResultVND(test, argv, initFO, Neighbor_getAssortment(neigh), sol, atoi(argv[4]));

    Neighbor_free(&neigh);
    VND_free( &vnd );
    Test_free( &test);
    Inst_free( &inst );

    return (0);
}

