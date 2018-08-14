
/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * mip_compact: routines for the creation of compact MIP formulations
 *              for the whole problem or for just one project
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef PARAMETERS
#define  PARAMETERS

#include "instance.h"

struct _Parameters {
    /*if the code will be running on exact or continuous mode */
    int continuous;
    char *name;
    char *name_danilo;

    /* cuts activated*/
    int cutJS;
    int cutPrec;
    int cutCGCPU;
    int cutRR;
    int cutCLIQUE;
    int cutODDHOLES;
    int cutDefaultCBC;
    int cutCG;
    int cutCGGPU;
    int cutCGGPUR2;
    int *roundCuts;

    /* parameters to CG cut*/
    int mininstant;
    int maxinstant;
    int jump;

    double maxReducedCost;
    int maxNode;
    int lifting;
    double slack;
    double maxcuts;
    int  typepsp;

    const Instance *inst;
} typedef Parameters;

//typedef struct _Parameters Parameters;
int Par_getTypePSP( Parameters *par );
Parameters *Par_create(Instance *inst,  char **argv, int argc );
void Par_free( Parameters **par );
void Par_checkArgs(Parameters *par, char **argv, int argc);


int Par_getCutPrec( Parameters *par );
int Par_getCutJS( Parameters *par );
int Par_getCutCGCPU( Parameters *par );
int Par_getCutCG( Parameters *par );
int Par_getCutCGGPU( Parameters *par );
int Par_getCutCGGPUR2( Parameters *par );
int Par_getCutRR( Parameters *par );
int Par_getCutCLIQUE( Parameters *par );
int Par_getCutODDHOLES( Parameters *par );
int Par_getMaxInstant( Parameters *par );
int Par_getJump( Parameters *par );
int Par_getCutDefaultCBC( Parameters *par );
double Par_getSlack(  Parameters *par );
char *Par_getName(  Parameters *par );
//danilo
char *Par_getNameDanilo(  Parameters *par );
double Par_getMaxCut(  Parameters *par );
int Par_getContinuous( Parameters *par);
int Par_getLifting( Parameters *par);
double Par_getMaxReducedCost( Parameters *par);
#endif // PARAMETERS

