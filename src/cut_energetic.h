

/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * mip_compact: routines for the creation of compact MIP formulations
 *              for the whole problem or for just one project
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef CUT_ENERGETIC
#define  CUT_ENERGETIC
#include "cut_pool.h"

struct _CutE {

    VecInt **cutElem; /*index of cuts elements */
    VecDbl **cutCoef; /* coefficient of cuts elements */
    VecDbl *cutrhs; /* rhs of cuts */
    VecDbl *cutviolation; /* violation of cuts */
    VecInt *cutsense; /* sense of cuts */
    VecInt *cutnelem; /* size of cuts */
    VecInt *cutdominated; /* if cut is dominated or not (1,0)*/
    VecStr *cutname; /* name of cuts */

    int nMaxWindowWithCutEnergeticViolated; /*max tw found violated cuts */
    double nMaxWindowWithCutEnergeticViolatedPerc; /*max tw in Perc found violated cuts */

    int nS; /*start time window*/
    int nE; /*end time window*/


    int nMinWindowWithCutEnergeticViolated; /*min tw found violated cuts */
    double nMinWindowWithCutEnergeticViolatedPerc; /*min tw in Perc found violated cuts */

    int nmS; /*start time window*/
    int nmE; /*end time window*/

    int nAllocations;
    const Instance* inst;
};

typedef struct _CutE CutE;
void CutE_add_cuts_energetic_guloso_parallel( CutE *cen,  LinearProgram *lp, const int *origCols,  LinearProgram *origLP, const Instance *inst, double timeLeft, int continuous,  int **maxTJM,  double maxcuts, int nround );
void CutE_add_cuts_energetic_parallel(  CutE *cE, LinearProgram *lp, const int *origCols, LinearProgram *origLP, const Instance *inst, double timeLeft, int continuous,  int **maxTJM,  double maxcuts,  int nround );
void CutE_add_cuts_energetic_guloso_interval (CutE *cen,  LinearProgram *lp, const int *origCols,  LinearProgram *origLP, const Instance *inst, double timeLeft, int continuous,  int **maxTJM,  double maxcuts, int nround, int s, int e );
CutE *CutE_create(const Instance *inst);
void CutE_free( CutE **cutE );

#endif // CUT_ODDHOLES

