/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include "instance.h"
#include "solution.h"
#include "macros.h"
#include "mip_compact.h"
#include "time.h"
#include "lp.h"

int main( int argc, char **argv )
{

    double startT = omp_get_wtime();

    if ( argc < 9 ) {
        fprintf( stderr, "ERROR! Enter path, instance, _nameFile, time, _nameSolIni, numSol, sumTPD, code(1: exact 2: AsContinuousCuts, 3: heuristic)");
        fflush( stderr );
        exit( EXIT_FAILURE );
    }

    printf("%s\n",argv[3]);
    fflush(stdout);
    fflush(stderr);
    int tmp = atoi( argv[4] );

    int seed = 100000;

    srand(seed);
    Instance *inst = Inst_read( argv, argc );

    //Inst_print(inst);
    //getchar();
    int sumtpd = atoi( argv[7] );
    double _time;
    _time = ( (double) tmp - (omp_get_wtime()-startT) );
    if(_time <= 0) {
        printf( "\nMipCreate: Time is over %f \n", _time);
        exit(0);
    }

    int code = atoi( argv[8] );
    MIPCompact *mipC = NULL;

    if(code == 1 )
        mipC = MipC_runExact(inst, argv, argc, sumtpd, ( (double) tmp - (omp_get_wtime()-startT) ));
    else if(code == 2 )
        mipC = MipC_run(inst, argv, argc, sumtpd, ( (double) tmp - (omp_get_wtime()-startT) ));
    else if(code == 3 )
        mipC = MipC_runHeuristc(inst, argv, argc, sumtpd, ( (double) tmp - (omp_get_wtime()-startT) ));

    /*Writing solution*/

    int cg = MipC_getCutCG( mipC );
    int cp = MipC_getCutPrec( mipC );
    int cr = MipC_getCutRR( mipC );
    int cc = MipC_getCutCLIQUE( mipC );
    double slk = MipC_getSlack( mipC );
    double mc = MipC_getMaxCut(mipC);
    int mi = MipC_getMaxInstant(mipC);
    int ju = MipC_getJump(mipC);

    Sol_write( MipC_getSol(mipC), argv[3] );
    char filelp[256]=" ";
    if(cg) {
        sprintf( filelp, "%s_cg_%d_cp_%d_cr_%d_cc_%d_l_%d_mc_%f_rc_%f_slk_%f_mi_%d_ju_%d.lp",  argv[3], cg,cp,cr,cc, MipC_getLifting(mipC),  mc, MipC_getMaxReducedCost(mipC), slk,mi,ju);
        printf("\n%s_cg_%d_cp_%d_cr_%d_cc_%d_l_%d_mc_%f_rc_%f_slk_%f_mi_%d_ju_%d\n",  argv[3], cg,cp,cr,cc, MipC_getLifting(mipC),  mc, MipC_getMaxReducedCost(mipC), slk,mi,ju);
    } else {
        sprintf( filelp, "%s_cg_%d_cp_%d_cr_%d_cc_%d_l_%d_mc_%f_rc_%f_slk_%f.lp",  argv[3], cg,cp,cr,cc, MipC_getLifting(mipC),  mc, MipC_getMaxReducedCost(mipC), slk);
        printf("\n%s_cg_%d_cp_%d_cr_%d_cc_%d_l_%d_mc_%f_rc_%f_slk_%f\n",  argv[3], cg,cp,cr,cc, MipC_getLifting(mipC),  mc, MipC_getMaxReducedCost(mipC), slk );
    }
    lp_write_lp(MipC_mip(mipC),filelp);

    char fileName[300] = {"results.txt"};
    for(int n = 0; n < argc ; n++)
        if (strcmp(argv[n], "-filename") == 0) {
            strcpy(fileName, argv[n+1]);
            break;
        }

    FILE *fp = fopen( fileName, "a" );
    if ( fp == NULL ) {
        printf( "File was not opened. : path: %s\n", fileName);
        exit( 0 );
    }

    MipC_writeArgs(mipC, fp, argv, argc,  startT, sumtpd);

    fclose( fp );

    MipC_free(&mipC);

    Inst_free( &inst );

    lp_close_env();
#ifdef DEBUG
    printf("Total time taken by CPU: %f on instance %s \n", (double) omp_get_wtime()-startT, argv[2] );
    //   getchar();
#endif // DEBUG
    return 0;
}

