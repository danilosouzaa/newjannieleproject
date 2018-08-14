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
#include "parameters.h"
#include "mip_compact.h"
#include "time.h"
#include "lp.h"


int Test_fileExists(char fileName[])
{

    FILE *fp;

    //char mainPath[500] = "";
    //strcat(mainPath, fileName);
    //printf("\nmainPath: %s", mainPath);

    fp=fopen(fileName,"r");

    if(fp) {
        fclose(fp);

        printf("\nErro in fileExists!");

        return 1;
    } else
        return 0;


    return -1;

}



int Write_filesValidation(Instance *inst, int argc, char **argv )
{

    FILE *fp;

    char filevalidation[500] = "";
    strcat(filevalidation, argv[2]);
    strcat(filevalidation, ".txt");

    char pathandfile[500] = "";
    strcat(pathandfile, argv[3]);
    strcat(pathandfile, argv[2]);
    printf("\npathandfile: %s", pathandfile);

    fp=fopen(filevalidation,"a");

    if(fp) {

        fprintf(fp,"1\n");
        const Project *proj = Inst_project(inst,0);
        fprintf(fp,"%d\n", Project_releaseDate(proj));
        fprintf(fp,"%d\n", Project_criticalPath(proj));
        fprintf(fp,"%s\n", pathandfile);
        fprintf(fp,"%d\n", Inst_nResR(inst));
        for(int r = 0 ;r < Inst_nResR(inst) ; r++)
            fprintf(fp,"-1 ");
        fclose(fp);
    }

    exit(0);

}

int main( int argc, char **argv )
{
    double startT = omp_get_wtime();

    if ( argc < 7 ) {
        fprintf( stderr, "ERROR! Enter path, instance, _nameFileSol, time, sumTPD, sumTMS");
        fflush( stderr );
        exit( EXIT_FAILURE );
    }


    int tmp = atoi( argv[4] );

    if(Test_fileExists(argv[3]))
        exit(0);
    else {
        printf("%s\n",argv[3]);
        fflush(stdout);
    }

    int seed = 100000;

    srand(seed);

    Instance *inst = Inst_read( argv, argc );
   // Write_filesValidation(inst, argc, argv );

    Solution *sol = Sol_create(inst);
    Inst_print(inst); //getchar();
    Parameters *par = Par_create(inst, argv, argc);
    Results *res = Res_create();


    int sumtpd = atoi( argv[5] );
    int sumtms = atoi( argv[6] );
    double _time = ( (double) tmp - (omp_get_wtime()-startT) );
    MIPCompact *mipC = MipC_create(inst, par,  -1, sumtpd, sumtms, _time);
    _time = ( (double) tmp - (omp_get_wtime()-startT) );
    //printf("mipC->tpdSum %d\n",MipC_getSumTPD(mipC));
    MipC_cutting_plane(mipC,sol,res,_time);
    //MipC_cutting_plane_cgraph(mipC,sol,res,_time);
   // MipC_linear_relaxation( mipC,sol,res,_time);
    //MipC_solve( mipC, _time);

#ifdef DEBUG
    printf(" Search - best objective %g ", MipC_getBestObj(mipC));
#endif

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

    Res_writeResults(par, MipC_mip(mipC), MipC_getSol(mipC), res, MipC_getBestObj(mipC), fp, argv, argc,  startT, sumtpd,sumtms);


    fclose( fp );

    MipC_free( &mipC );
    free(mipC);
    Inst_free( &inst );
    Sol_free( &sol );
    Res_free( &res);

#ifdef DEBUG
    printf("Total time taken by CPU: %f on instance %s \n", (double) omp_get_wtime()-startT, argv[2] );
    //   getchar();
#endif // DEBUG
    return 0;
}

