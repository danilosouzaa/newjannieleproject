/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <stdlib.h>
#include <stdio.h>
#include "instance.h"
#include "constructive.h"
#include "ms_solver_lahc.h"
#include "vnd.h"
#include "neighborhood.h"
#include "solution.h"
#include "macros.h"
#include "time.h"

int main( int argc, char **argv )
{
    clock_t t_begin, t_end;
    double _time = 0;
    t_begin = clock();

    int LFA = atoi( argv[4] );
    int ITERAT = atoi( argv[5] );
    int tmp = atoi( argv[5]);

    if ( argc<5 ) {
        fprintf( stderr, "ERROR! Enter path, instance, _nameFile, imp, time");
        fflush( stderr );
        exit( EXIT_FAILURE );
    }


    Instance *inst = Inst_create( argv[1], argv[2] );

    Solution *sol = Sol_create(inst);

    /*Creating initial solution*/
    Constructive *cons = Cons_create(inst, sol, argv);
    Cons_run(cons);

    /*Writing solution*/
    Sol_write( sol,   argv[3] );

    t_end = clock();
    _time = ( (double)( t_end-t_begin )/CLOCKS_PER_SEC );
    printf( "\nCost Instance: %ld |||| totalTime = %f | Infact = %d\n", Sol_getCost(sol), _time, Modes_inf( Sol_getModeSet( sol ) ) );

    Cons_free( &cons );
    Sol_free( &sol );
    Inst_free( &inst );


    return EXIT_SUCCESS;
}

