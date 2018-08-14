/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <stdlib.h>
#include <stdio.h>
#include "instance.h"
#include "ms_solver_lahc.h"
#include "neighborhood.h"
#include "solution.h"
#include "macros.h"
#include "time.h"

//#define TIME 20
//#define LFA 10

int main( int argc, char **argv )
{
    clock_t t_begin, t_end;
    double _time = 0;
    t_begin = clock();

    /* random seed */
    srand(time(NULL));

    int LFA = atoi( argv[4] );
    int ITERAT = atoi( argv[5] );
    double vetIntens[2];
    Cost *vetCost;
    int *vetInf;

    if ( argc<8 ) {
        fprintf( stderr, "ERROR! Enter path, instance, _nameFile, LFA, ITERATIONS, intens0, intens1");
        fflush( stderr );
        exit( EXIT_FAILURE );
    }

    for( int i=0; i<2; ++i )
        vetIntens[i] = (double)( atof( argv[i+6] ) );

    Instance *inst = Inst_create( argv[1], argv[2] );

    //Inst_print( inst );

    Solution *sol = Sol_create(inst);

    int firstJob = 0, lastJob =0;

    ALLOCATE_VECTOR( vetCost, Cost, Inst_nProjects( inst ) );
    ALLOCATE_VECTOR( vetInf, int, Inst_nProjects( inst ) );


    for(int p = 0; p <Inst_nProjects(inst); ++p) {
        const Project* proj = Inst_project(inst,p);
        printf("\n\nProject %d\n", p);
        lastJob = firstJob + Project_nJobs(proj)-1;
        MSSolverLAHC *msLahc = MS_create(inst,LFA, firstJob, lastJob);

        for ( int i=0; i<2; ++i )
            MS_setNeighborhoodPriority( msLahc, i, vetIntens[i] );

        MS_run(msLahc, ITERAT);
        firstJob = lastJob+1;
        Sol_setModes(sol, MS_bestModes(msLahc));

        vetCost[p] = Modes_cost( MS_bestModes( msLahc ) );
        vetInf[p] = Modes_inf( MS_bestModes( msLahc ) );

        MS_free( &msLahc );
    }

    Sol_build(sol);

    Sol_write( sol,  argv[3] );

    t_end = clock();
    _time = ( (double)( t_end-t_begin )/CLOCKS_PER_SEC );
    printf( "\nCost Instance: %ld |||| totalTime = %f | Infact = %d\n", Sol_getCost(sol), _time, Modes_inf( Sol_getModeSet( sol ) ) );

    Results_write( inst, argv[2], LFA, ITERAT, vetIntens, 2, Modes_cost( Sol_getModeSet( sol ) ), Modes_inf( Sol_getModeSet( sol ) ), Inst_nProjects( inst ), vetCost, vetInf, _time );

    Sol_free( &sol );
    Inst_free( &inst );
    free(vetCost);
    free(vetInf);


    return EXIT_SUCCESS;
}

