/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <stdlib.h>
#include <stdio.h>
#include "instance.h"
#include "macros.h"
#include "info_instance.h"

int main( int argc, char **argv )
{

  //  Instance *inst = Inst_create( argv[1], argv[2] );
    //Instance *inst = Inst_read( argv, argc );
    //int sumJobProj =0;
    char filelp[256];
    fprintf(filelp, "%s.txt", argv[2]);

    FILE *fp = fopen( filelp, "a" );
    if ( fp == NULL ) {
        printf( "File was not opened. : path: %s\n", filelp);
        exit( 0 );
    }

    int rd=0, cp=0;
    fprintf(fp, "1\n");
   /* for(int p = 0 ; p < Inst_nProjects(inst) ; p++) {
        Project *proj = Inst_project(inst,p);
        rd = Project_releaseDate(proj);
        cp = Project_criticalPath(proj);
    }*/
    fprintf(fp, "%d\n", rd);
    fprintf(fp, "%d\n", cp);
    fprintf(fp, "%s%s\n", argv[1], argv[2]);
    fprintf(fp, "4\n");
    fprintf(fp, "-1 -1 -1 -1");
    fprintf(fp, " \n");


    exit(0);
   /* char filelp[256]="InfoInstanceCP.txt";

    FILE *fp = fopen( filelp, "a" );
    if ( fp == NULL ) {
        printf( "File was not opened. : path: %s\n", filelp);
        exit( 0 );
    }

    fprintf(fp, "%s ", argv[2]);
    for(int p = 0 ; p < Inst_nProjects(inst) ; p++) {
        Project *proj = Inst_project(inst,p);
        fprintf(fp, " ; %d ; %d", Project_criticalPath(proj), Project_releaseDate(proj));
    }
    fprintf(fp, " \n");*/

      /*Writing solution*/



  /*  printf( "instance, nProjects, minJobsProj, maxJobsProj, avgJobsProj, starts, avgStartProj, ends, avgEndProj, nJobs,  minDuration, maxDuration, avgDuration, minModes, maxModes, avgModes, minNumRRModes, maxNumRRModes, avgNumRRModes, minNumNRModes, maxNumNRModes, avgNumNRModes, minPrec, maxPrec, avgPrec, nRR, nNR, minConsumptionRR, maxConsumptionRR,  avgConsumptionRR, minConsumptionNR, maxConsumptionNR, avgConsumptionNR, maxCapRR, minCapRR, avgCapRR, maxCapNR, minCapNR, avgCapNR\n");


    Inst_print(inst);

    InfoInst *iinst = IInst_create();
    IInst_calculate(inst, iinst);
  //  IInst_writeSumTPD(iinst,argv);
    IInst_wrieInf(iinst,argv);
    IInst_free(&iinst);

    Inst_free( &inst );*/

    return EXIT_SUCCESS;
}
