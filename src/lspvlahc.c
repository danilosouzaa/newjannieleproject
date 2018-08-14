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
#include <omp.h>


int main( int argc, char **argv )
{

    double startT = omp_get_wtime();

    if ( argc < 5 ) {
        fprintf( stderr, "ERROR! Enter path, instance, _nameFile, time, _nameSolIni, seed, nThreads");
        fflush( stderr );
        exit( EXIT_FAILURE );
    }

    // if(Test_fileExists(argv[3]))
    //   exit(0);

    int seed = atoi( argv[6] );

    srand(seed);

    const int nthreads = atoi( argv[7] );
    omp_set_num_threads( nthreads );

    Instance *inst = Inst_create( argv[1], argv[2] );
    Solution *globalBest = Sol_create(inst);

    //read initial solution
    Sol_read(globalBest, argv[5]);


    Cost initFO = Sol_getCost(globalBest);


#ifdef DEBUG
    printf( "time spent before LAHC = %.2f seconds.\n", omp_get_wtime()-startT );
    printf("\nstarting LAHC with cost: %ld ...\n\n", initFO);
#endif // DEBUG

    double timeLimit = (double)(atoi(argv[4]));
    const double timeRem = timeLimit - ( omp_get_wtime() - startT );

    #pragma omp parallel
    {
        Neighborhood *neigh = Neighbor_create(inst, argv, argc);

        Test *test = Test_create(14,inst);

        Solution *sol = Sol_create(inst);
        Sol_cpy(sol, globalBest);

        /*Local Search*/
        LAHC *lahc = LAHC_create(inst, sol, neigh, argv, argc);
        LAHC_run_parallel(lahc, neigh, timeRem, globalBest, test);

        Test_setTotalTime(test, omp_get_wtime() - startT);
        Test_writeResultLAHC_thread(test, argv, initFO, neigh, sol, LAHC_getLfa(lahc), atoi(argv[8]), seed, omp_get_thread_num());

        LAHC_free( &lahc );
        Test_free( &test);
        Neighbor_free(&neigh);

    }

#ifdef DEBUG
    printf("\nsearch concluded in %.3f seconds.\n", omp_get_wtime()-startT );
    fflush( stdout );
#endif // DEBUG

    /*Writing solution*/
    Sol_write( globalBest,   argv[3] );


#ifdef DEBUG
    printf( "\nCost Instance: %ld  TPD %ld TMS %ld   |||| totalTime = %f \n", Sol_getCost(globalBest), Sol_getTPD(globalBest), Sol_getTMS(globalBest), omp_get_wtime()-startT );
#endif // DEBUG

    Sol_free(&globalBest);
    Inst_free( &inst );


    return EXIT_SUCCESS;
}

