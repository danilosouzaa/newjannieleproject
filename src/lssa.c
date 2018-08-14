/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <stdlib.h>
#include <stdio.h>
#include "instance.h"
#include "sa.h"
#include "neighborhood.h"
#include "solution.h"
#include "macros.h"
#include "time.h"
#include "test.h"
#include <omp.h>    // Setar -lgomp em Projects > Build Options> Linker settings > Other linker options

int main( int argc, char **argv )
{
    //clock_t t_begin, t_end;
    double t_begin, t_end;
    double _time = 0;
    //t_begin = clock();
    t_begin = omp_get_wtime();


    if ( argc< 5 ) {
        fprintf( stderr, "ERROR! Enter path, instance, _nameFile, time, _nameSolIni");
        fflush( stderr );
        exit( EXIT_FAILURE );
    }



    int tmp = atoi( argv[4]);

    srand(100000);

    Instance *inst = Inst_create( argv[1], argv[2] );

    Solution *sol = Sol_create(inst);

    Sol_read(sol, argv[5]);

    Cost initFO = 0;
    initFO = Sol_getCost(sol);

    Test *test = Test_create(14);

    Neighborhood *neigh = Neighbor_create(inst, argv, argc);

    /*Local Search*/
    SA *sa = SA_create(inst, neigh, sol, argv, argc);

    Test_setT(test, SA_getT(sa));
    Test_setSAmax(test, SA_getSAmax(sa));
    Test_setAlpha(test, SA_getAlpha(sa));

    _time = ( (double) tmp - (double)(clock()-t_begin )/CLOCKS_PER_SEC );


    SA_run(sa, neigh, _time, test);

    /*Writing solution*/
    Sol_write( sol,   argv[3] );

    t_end = omp_get_wtime();

    _time = (double)( t_end-t_begin );

#ifdef DEBUG
    printf( "\nCost Instance: %ld  TPD %ld TMS %ld   |||| totalTime = %f \n", Sol_getCost(sol), Sol_getTPD(sol), Sol_getTMS(sol),  _time );
#endif // DEBUG

    Test_setTotalTime(test, _time);

    Test_writeResultSA(test, argv, initFO, Neighbor_getAssortment(neigh), sol);


    Test_free( &test);
    Neighbor_free(&neigh);
    SA_free( &sa );
    Inst_free( &inst );


    return EXIT_SUCCESS;
}

