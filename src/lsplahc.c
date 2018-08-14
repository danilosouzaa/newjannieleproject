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

    int seed = atoi( argv[6] );

    srand(seed);

    const int nthreads = atoi( argv[7] );
    omp_set_num_threads( nthreads );

    Instance *inst = Inst_create( argv[1], argv[2] );
    Solution *bestSol = Sol_create(inst);

    //read initial solution
    Sol_read(bestSol, argv[5]);


    Cost initFO = Sol_getCost(bestSol);



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
        Sol_cpy(sol, bestSol);

        /*Local Search*/
        LAHC *lahc = LAHC_create(inst, sol, neigh, argv, argc);
        LAHC_run(lahc, neigh, timeRem, test,argv[2]);
        #pragma omp critical
        {
            int th_id = omp_get_thread_num();
#ifdef DEBUG
            printf("\nThread: %d Cost %ld \n", th_id, Sol_getCost(sol) );
#endif // DEBUG
            if (Sol_getCost(sol) < Sol_getCost(bestSol)  )
                Sol_cpy(bestSol, sol);
        }

        LAHC_free( &lahc );
        Test_free( &test);
        Neighbor_free(&neigh);
    }

#ifdef DEBUG
    printf("\nsearch concluded in %.3f seconds.\n", omp_get_wtime()-startT );
    fflush( stdout );
#endif // DEBUG

    /*Writing solution*/
    Sol_write( bestSol,   argv[3] );


#ifdef DEBUG
    printf( "\nCost Instance: %ld  TPD %ld TMS %ld   |||| totalTime = %f \n", Sol_getCost(bestSol), Sol_getTPD(bestSol), Sol_getTMS(bestSol), omp_get_wtime()-startT );
#endif // DEBUG


    Sol_free(&bestSol);
    Inst_free( &inst );


    return EXIT_SUCCESS;
}

