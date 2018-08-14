/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include "parameters.h"
#include "macros.h"
#include "lp.h"
#include <limits.h>
#include <assert.h>
#include <string.h>

Parameters *Par_create(Instance *inst,  char **argv, int argc)
{
    Parameters *par;
    ALLOCATE_INI( par, Parameters );

    ALLOCATE_VECTOR_INI( par->name, char, 256 );
    sprintf(par->name, "%s", argv[3]);
    ALLOCATE_VECTOR_INI( par->name_danilo, char, 256 );
    sprintf(par->name_danilo, "%s", argv[2]);
    par->inst = inst;
    par->lifting = 0;


    par->cutPrec = 0;
    par->cutCGCPU = 0;
    par->cutRR = 0;
    par->cutCLIQUE = 0;
    par->cutCG = 0;
    par->cutCGGPU = 0;
    par->cutCGGPUR2 = 0;
    par->cutDefaultCBC = 0;
    ALLOCATE_VECTOR_INI( par->roundCuts, int, LP_CUT_TYPES);
    FILL( par->roundCuts,0,LP_CUT_TYPES,1);

    par->maxinstant = 0;
    par->maxinstant = INT_MAX;
    par->jump = 0;

    par->maxNode = 0;
    par->maxReducedCost = -1.0;
    par->maxcuts = Inst_nJobs(inst);
    par->slack = 0.0;
    par->typepsp = 0; //0: single mode 1: multimode 2:multiproj

    Par_checkArgs(par, argv, argc);

    return par;
}


void Par_checkArgs(Parameters *par, char **argv, int argc)
{

    assert(par != NULL);

    for(int n = 1; n< argc ; n++) {
        //  printf("argc %d n %d \n",argc, n);fflush(stdout);
        if (strcmp(argv[n],"-rcpsp") == 0) {
            par->typepsp =  0;
            //            printf("par->typepsp %d inside\n",par->typepsp); fflush(stdout);
            continue;
        }
        if (strcmp(argv[n],"-mrcpsp") == 0) {
            par->typepsp =  1;
            //            printf("par->typepsp %d inside\n",par->typepsp); fflush(stdout);
            continue;
        }
        if (strcmp(argv[n],"-mmrcmpsp") == 0) {
            par->typepsp =  2;
            //            printf("par->typepsp %d inside\n",par->typepsp); fflush(stdout);
            continue;
        }
        if (strcmp(argv[n],"-continuous") == 0) {
            n++;
            par->continuous =  atoi(argv[n]);
            //   printf("par->continuous %d inside",par->continuous); fflush(stdout);
            continue;
        }
        if (strcmp(argv[n],"-cutRR") == 0) {
            n++;
            par->cutRR =  atoi(argv[n]);
            //   printf("par->cutRR %d inside\n",par->cutRR); fflush(stdout);
            continue;
        }
        if (strcmp(argv[n],"-cutCLIQUE") == 0) {
            n++;
            par->cutCLIQUE =  atoi(argv[n]);
            //printf("par->cutCLIQUE %d inside\n",par->cutCLIQUE); fflush(stdout);
            continue;
        }
        if (strcmp(argv[n],"-cutODDHOLES") == 0) {
            n++;
            par->cutODDHOLES=  atoi(argv[n]);
            //   printf("par->cutODDHOLES %d inside\n",par->cutODDHOLES); fflush(stdout);
            continue;
        }
        if (strcmp(argv[n],"-cutPREC") == 0) {
            n++;
            par->cutPrec =  atoi(argv[n]);
            //  printf("par->cutPrec %d inside\n",par->cutPrec); fflush(stdout);
            continue;
        }
         if (strcmp(argv[n],"-cutJS") == 0) {
            n++;
            par->cutJS =  atoi(argv[n]);
            //  printf("par->cutJS %d inside\n",par->cutJS); fflush(stdout);
            continue;
        }
        if (strcmp(argv[n],"-cutCGCPU") == 0) {
            n++;
            par->cutCGCPU =  atoi(argv[n]);
              //printf("par->cutCGCPU %d inside\n",par->cutCGCPU); fflush(stdout);
            continue;
        }
        if (strcmp(argv[n],"-cutCG") == 0) {
            n++;
            par->cutCG =  atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-cutCGGPU") == 0) {
            n++;
            par->cutCGGPU =  atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-cutCGGPUR2") == 0) {
            n++;
            par->cutCGGPUR2 =  atoi(argv[n]);
            //   printf("par->cutCGGPUR2 %d inside\n",par->cutCGGPUR2); fflush(stdout);
            continue;
        }
        if (strcmp(argv[n],"-cutDefaultCBC") == 0) {
            n++;
            par->cutDefaultCBC =  atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-gomory") == 0) {
            n++;
            par->roundCuts[LPC_GOMORY] =  atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-reduce") == 0) {
            n++;
            par->roundCuts[LPC_REDUCE] =  atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-mir") == 0) {
            n++;
            par->roundCuts[LPC_MIR] =  atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-twomir") == 0) {
            n++;
            par->roundCuts[LPC_TWO_MIR] =  atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-landp") == 0) {
            n++;
            par->roundCuts[LPC_L_AND_P] =  atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-zerohalf") == 0) {
            n++;
            par->roundCuts[LPC_ZERO_HALF] =  atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-knapsack") == 0) {
            n++;
            par->roundCuts[LPC_KNAPSACK] =  atoi(argv[n]);
            continue;
        }

        if (strcmp(argv[n],"-flow") == 0) {
            n++;
            par->roundCuts[LPC_FLOW] =  atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-maxNode") == 0) {
            n++;
            par->maxNode =  atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-maxReducedCost") == 0) {
            n++;
            par->maxReducedCost =  atof(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-lifting") == 0) {
            n++;
            par->lifting =  atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-slack") == 0) {
            n++;
            par->slack =  atof(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-maxcuts") == 0) {
            n++;
            double value = atof(argv[n]);
            par->maxcuts =  value ;
            continue;
        }
        if (strcmp(argv[n],"-maxinstant") == 0) {
            n++;
            par->maxinstant =  atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-jump") == 0) {
            n++;
            par->jump =  atoi(argv[n]);
            continue;
        }
    }
}


void Par_free( Parameters **_par )
{
    Parameters *par = *_par;

    free(par->name);
    free(par->name_danilo);
    free(par->roundCuts);
    free( par );
    *_par = NULL;
}


int Par_getCutJS( Parameters *par )
{
    return par->cutJS;
}


int Par_getCutPrec( Parameters *par )
{
    return par->cutPrec;
}

int Par_getCutCGCPU( Parameters *par )
{
    return par->cutCGCPU;
}

int Par_getCutCG( Parameters *par )
{
    return par->cutCG;
}

int Par_getTypePSP( Parameters *par )
{
    return par->typepsp;
}

int Par_getCutCGGPU( Parameters *par )
{
    return par->cutCGGPU;
}

int Par_getCutCGGPUR2( Parameters *par )
{
    return par->cutCGGPUR2;
}

int Par_getCutRR( Parameters *par )
{
    return par->cutRR;
}

int Par_getCutCLIQUE( Parameters *par )
{
    return par->cutCLIQUE;
}

int Par_getCutODDHOLES( Parameters *par )
{
    return par->cutODDHOLES;
}

int Par_getMaxInstant( Parameters *par )
{
    return par->maxinstant;
}

int Par_getJump( Parameters *par )
{
    return par->jump;
}

int Par_getCutDefaultCBC( Parameters *par )
{
    return par->cutDefaultCBC;
}

double Par_getSlack(  Parameters *par )
{
    return par->slack;
}

double Par_getMaxCut(  Parameters *par )
{
    return par->maxcuts;
}


char *Par_getName(  Parameters *par )
{
    return par->name;
}

char *Par_getNameDanilo(  Parameters *par )
{
    return par->name_danilo;
}


int Par_getContinuous( Parameters *par)
{
    return par->continuous;
}

int Par_getLifting( Parameters *par)
{
    return par->lifting;
}

double Par_getMaxReducedCost( Parameters *par)
{
    return par->maxReducedCost;
}
