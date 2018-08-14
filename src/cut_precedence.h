

/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * mip_compact: routines for the creation of compact MIP formulations
 *              for the whole problem or for just one project
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef CUT_PRECEDENCE
#define  CUT_PRECEDENCE
#include "cut_pool.h"
#include "gurobi_c.h"

struct _CutPR {

    VecInt **cutElem; /*index of cuts elements */
    VecDbl **cutCoef; /* coefficient of cuts elements */
    VecDbl *cutrhs; /* rhs of cuts */
    VecDbl *cutviolation; /* violation of cuts */
    VecInt *cutsense; /* sense of cuts */
    VecInt *cutnelem; /* size of cuts */
    VecStr *cutname; /* name of cuts */
    VecInt *cutdominated; /* if cut is dominated or not (1,0)*/
    int nAllocations;
    //double timeseparation;
    //int nAddCuts;
    const Instance* inst;
};

typedef struct _CutPR CutPR;
CutPR *CutPR_create(const Instance *inst);

//void CutPR_add_cuts_precedence_parallel( CutPR *cpr, LinearProgram *lp, const int *origCols, LinearProgram *origLP, const Instance *inst, double timeLeft, int continuous, int lifting, int **maxTJM, int *maxTJ, double maxcuts, int nround );
void CutPR_add_cuts_precedence_parallel( CutPR *cpr, LinearProgram *lp, const int *origCols, LinearProgram *origLP, const Instance *inst, double timeLeft, int continuous, int lifting, double maxcuts, int nround );

void CutPR_add_cuts_precedence_callback(  GRBmodel *model,   void *cbdata, const Instance *inst, double *timeLeft, int lifting, double maxcuts, int  *nCuts);
//void CutPR_add_cuts_precedence_callback(  GRBmodel *model,  void *cbdata, const Instance *inst, double timeLeft, int lifting, double maxcuts, int *nCuts);
void CutPR_free( CutPR **cutPR );

#endif // CUT_PRECEDENCE

