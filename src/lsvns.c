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
#include "vns.h"
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


    if ( argc < 8 ) {
        fprintf( stderr, "ERROR! Enter path, instance, _nameFile, time, _nameSolIni, firstImproved(0,1), type(1BVNS,2VNSRNA), numSol");
        fflush( stderr );
        exit( EXIT_FAILURE );
    }

    if(Test_fileExists(argv[3]))
        exit(0);

    int tmp = atoi( argv[4] );
    int type = atoi( argv[7] );

    int seed = 100000;

    srand(seed);
    Instance *inst = Inst_create( argv[1], argv[2] );

    Solution *sol = Sol_create(inst);

    /*read initial solution*/
    Sol_read(sol, argv[5]);

    Cost initFO = 0;
    initFO = Sol_getCost(sol);

    Test *test = Test_create(14, inst);

    Neighborhood *neigh = Neighbor_create(inst, argv, argc);

    /*Local Search*/
    VNS *vns = VNS_create(inst, sol, neigh, argv, argc);

    _time = ( (double) tmp - (double)(clock()-t_begin )/CLOCKS_PER_SEC );

    /*default vnd, best improved*/

    if(type == 1)
        //0.717369 BI // 0.723171 FI
        VNS_run_general_vnd(vns, neigh, _time, 14, atoi(argv[6]), test, argv, argc );
    if(type == 2)
        //0.903
        VNS_run_RNA(vns, neigh, _time, 14, test, argv, argc );
    if(type == 3)
        VNS_run_RNA_shake(vns, neigh, _time, 14, test, argv, argc );
    if(type == 4)
        VNS_run_RNA_shake2(vns, neigh, _time, 14, test, argv, argc );
    if(type == 5)
        VNS_run_reduced(vns, neigh, _time, 14, atoi(argv[6]), test, argv, argc );
    if(type ==6 ) {
        //0.930
        VNS_run_LAHC_shake2(vns, neigh, _time, 14, test, argv, argc );
    }
    if(type ==7 )
        VNS_run_LAHC_smartshake(vns, neigh, _time, 14, argv, argc );

    /*Writing solution*/
    Sol_write( sol,   argv[3] );

    t_end = clock();
    _time = ( (double)( t_end-t_begin )/CLOCKS_PER_SEC );

#ifdef DEBUG
    printf( "\nCost Instance: %ld |||| totalTime = %f | Infact = %d\n", Sol_getCost(sol), _time, Modes_inf( Sol_getModeSet( sol ) ) );
    fflush(stdin);
    fflush(stdout);
#endif // DEBUG

    //Test_setTotalTime(test, _time);
    Test_writeResultVNS(test, argv, initFO, Neighbor_getAssortment(neigh),  atoi(argv[6]), atoi(argv[7]), sol, VNS_getLfa(vns),  VNS_getItRNA(vns), VNS_getItLAHC(vns), VNS_getNMoves(vns), VNS_getNSizeSamplingShake(vns),  VNS_getPerc(vns), VNS_getPercRS(vns), VNS_getDivTM(vns), VNS_getDivTJ(vns), VNS_getDivRM(vns), VNS_getDivRJ(vns), atoi(argv[8]));

    VNS_free( &vns );
    Neighbor_free( &neigh );
    Test_free( &test);
    Inst_free( &inst );

    return (0);
}

