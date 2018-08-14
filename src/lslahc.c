/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <stdlib.h>
#include <stdio.h>
#include "test.h"
#include "instance.h"
#include "lahc.h"
#include "macros.h"
#include "time.h"

int main( int argc, char **argv )
{
    clock_t t_begin, t_end;
    double _time = 0;
    t_begin = clock();


    if ( argc < 6 ) {
        fprintf( stderr, "ERROR! Enter path, instance, _nameFile, time, _nameSolIni, numSol");
        fflush( stderr );
        exit( EXIT_FAILURE );
    }


    if(Test_fileExists(argv[3]))
        exit(0);


    int tmp = atoi( argv[4] );

    int seed = 100000;

    srand(seed);

    //Instance *inst = Inst_create( argv[1], argv[2] );
    Instance *inst = Inst_read( argv, argc );

    Solution *sol = Sol_create(inst);

    /*read initial solution*/
    Sol_read(sol, argv[5]);

    Cost initFO = 0;
    initFO = Sol_getCost(sol);


    Neighborhood *neigh = Neighbor_create(inst, argv, argc);
    Test *test = Test_create(Neighbor_nNeighborhood(neigh), inst);

    /*Local Search*/
    LAHC *lahc = LAHC_create(inst, sol, neigh, argv, argc);

    _time = ( (double) tmp - (double)(clock()-t_begin )/CLOCKS_PER_SEC );
    LAHC_run(lahc, neigh, _time, test, argv[2]);

    /*Writing solution*/
    Sol_write( sol,   argv[3] );

    t_end = clock();
    _time = ( (double)( t_end-t_begin )/CLOCKS_PER_SEC );

#ifdef DEBUG
    printf( "\nCost Instance: %ld  TPD %ld TMS %ld   |||| totalTime = %f \n", Sol_getCost(sol), Sol_getTPD(sol), Sol_getTMS(sol),  _time );
#endif // DEBUG

    Test_setTotalTime(test, _time);

    //Test_readAnalisysNeigh(test, LAHC_getSW(lahc), LAHC_getItUp(lahc), LAHC_getPSW(lahc), argv[2]);
    //Test_readAnalisysNeighToInstance(test, LAHC_getSW(lahc), LAHC_getItUp(lahc), LAHC_getPSW(lahc), argv[2]);

    Test_writeResultLAHC(test, argv, initFO, neigh, sol, atoi(argv[6]), LAHC_getLfa(lahc), seed, LAHC_getNCostList(lahc), LAHC_getNDiversification(lahc), LAHC_getNStayDiversification(lahc), LAHC_getNWOImprove(lahc), LAHC_getPerc(lahc), LAHC_getSW(lahc), LA_getLearningRate(LAHC_getLA(lahc)), LAHC_getItUp(lahc), LAHC_getPSW(lahc), LAHC_getTM(lahc), LAHC_getTJ(lahc), LAHC_getRM(lahc), LAHC_getRJ(lahc));

    Test_free( &test);
    LAHC_free( &lahc );
    Neighbor_free(&neigh);
    Inst_free( &inst );

    return EXIT_SUCCESS;
}

