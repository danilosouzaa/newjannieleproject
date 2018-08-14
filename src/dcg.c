/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <omp.h>
#include "instance.h"
#include "proj_decomp.h"
#include "mip_compact.h"

int main( int argc, char **argv )
{
    if (argc<4) {
        fprintf( stderr, "usage:\n\tdcg instanceDir instance tpdSum\n");
        exit( EXIT_FAILURE );
    }

    Instance *inst = Inst_create( argv[1], argv[2] );

    /* tpd sum from best known solution, used to prune the time horizon     */
    const int tpdSum = atoi( argv[3] );

    /* getting optimal TPDs */

    for ( int i=0 ; i<Inst_nProjects( inst ) ; ++i ) {
        double start = omp_get_wtime();
        MIPCompact *mipP = MipC_create( inst, i, tpdSum );
        MipC_setMaxSeconds( mipP, 3600 );
        MipC_solve( mipP );
        double end = omp_get_wtime();
        const double sec = end-start;
        int tpd = 9999;
        if ( MipC_hasSolution( mipP ) )
            tpd = MipC_TPD( mipP );
        printf("%s %d %.2f %d\n", argv[2], i, sec, tpd );
    }

    Inst_free( &inst );

    return EXIT_SUCCESS;
}

