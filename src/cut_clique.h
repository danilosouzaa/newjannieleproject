

/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * mip_compact: routines for the creation of compact MIP formulations
 *              for the whole problem or for just one project
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef CUT_CLIQUE
#define  CUT_CLIQUE
#include "build_cgraph.h"
#include "clique_separation.h"
#include "clique_elite_set.h"
#include "cut_pool.h"



struct _CutCL {

    VecInt **cutElem; /*index of cuts elements */
    VecDbl **cutCoef; /* coefficient of cuts elements */
    VecDbl *cutrhs; /* rhs of cuts */
    VecDbl *cutviolation; /* violation of cuts */
    VecInt *cutsense; /* sense of cuts */
    VecInt *cutdominated; /* if cut is dominated or not (1,0)*/
    VecInt *cutnelem; /* size of cuts */
    VecStr *cutname; /* name of cuts */
    int nAllocations;
    // double timeseparation;
    // int nAddCuts;
    const Instance* inst;
};

typedef struct _CutCL CutCL;

void CutCL_add_cuts_conflicts_clique_parallel(  const CGraph *cgraph, CutCL *ccl, LinearProgram *lp, const Instance *inst, double timeLeft,  double maxcuts,  int nround );
void CutCL_add_cuts_clique_callback(  GRBmodel *model,   void *cbdata, const Instance *inst, double *timeLeft, int lifting, double maxcuts, int  *nCuts);
//void CutCL_add_cuts_clique_callback(  GRBmodel *model,   void *cbdata, const Instance *inst, double timeLeft, int lifting, double maxcuts, int  *nCuts);

CutCL *CutCL_create(const Instance *inst);

void CutCL_free( CutCL **cutCL );

#endif // CUT_CLIQUE

