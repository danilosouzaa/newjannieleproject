/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * mip_compact: routines for the creation of compact MIP formulations
 *              for the whole problem or for just one project
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef MIP_COMPACT
#define  MIP_COMPACT

#include "instance.h"
#include "solution.h"
#include "lp.h"
#include "parameters.h"
#include "results.h"
#include "cut_pool.h"

#define MIPC_ALL_PROJECTS -1

typedef struct _MIPCompact MIPCompact;


//-1 idxProj para todos, tpdSum e tmsSum passar o tpd e tms para todos os projetos.
MIPCompact *MipC_create( Instance *inst, Parameters *par, const int idxProj, const int tpdSum, const int tmsSum, double timeLeft );

/* parses a name in the format prefix(idx1,idx2,...) */
int parseName( const char *name, char *prefix, int *idx );

void MipC_solve( MIPCompact *mipp, double timeLeft );
double MipC_getSlack( MIPCompact *mipC );
double MipC_getMaxCut( MIPCompact *mipC );
int MipC_getRemoveCuts( MIPCompact *mipC );

/* solver the LP relaxation of the mip */
int MipC_linear_relaxation( MIPCompact *mipC, Solution *sol, Results *res,  double timeLeft);

/* solver the LP relaxation of the mip with cuts*/
int MipC_cutting_plane( MIPCompact *mipC, Solution *sol, Results *res,  double timeLeft);

/* tries to improve the current formulation by cutting the
   current fractional solution, returns the number of cuts added */
void MipC_setVariablesZero(MIPCompact *mipC);
int MipC_getMaxInstant( MIPCompact *mipC );
int MipC_getJump( MIPCompact *mipC );
int MipC_getMaxConstraint( MIPCompact *mipC );


void MipC_checkArgs(MIPCompact *mipC, char **argv, int argc);
void MipC_writeLP( const MIPCompact *mipp, const char *fileName );
char MipC_hasSolution( MIPCompact *MipP );
void MipC_allocateSol(MIPCompact *mipC);
int MipC_saveSol(MIPCompact *mipC);
int MipC_isInteger( MIPCompact *mipC);

/* gets and sets*/
int MipC_getnRounds();
int MipC_getMaxElements(int idx);
int MipC_getSumAllVarWithConflicts(int idx);
int MipC_getSumAllElements(int idx);
Solution* MipC_getSol( MIPCompact *mipC );
int MipC_getNCut( MIPCompact *mipC, int idx );
double MipC_getTCut( MIPCompact *mipC, int idx );
int MipC_getContinuous(MIPCompact *mipC);
int MipC_getLifting(MIPCompact *mipC);
double MipC_getBestObj( MIPCompact *mipC );
double MipC_getCurrentObj( MIPCompact *mipC );
double MipC_getBestPossibleObj( MIPCompact *mipC);
int MipC_getMaxCutElements(int idx);
int MipC_getSumElementsCut(int idx);
double MipC_TPD( MIPCompact *MipP );
LinearProgram* MipC_mip( MIPCompact *mipC);
const Instance* MipC_inst( MIPCompact *mipC );
void MipC_setInitialSolution( MIPCompact *mipp, const Solution *_bestSol );
void MipC_setMaxSeconds( MIPCompact *MipP, int maxSeconds );
//double computeSumXAllModesSecondPart( const Instance *inst,  int j, int ij, int t, int ***TJ, int **nCont, int start, const double ***x );
//double computeSumXAllModes( const Instance *inst, int j, int ***TJ, int **nCont, int start, int end, const double ***x );
void MipC_parseParameters( int argc, const char **argv );
Parameters *MipC_getPar(MIPCompact *mipC);

void MipC_help( );
void MipC_printConfig( );

void MipC_free( MIPCompact **mipp );

#endif // MIP_PROJECT

