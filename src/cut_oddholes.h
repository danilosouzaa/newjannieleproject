

/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * mip_compact: routines for the creation of compact MIP formulations
 *              for the whole problem or for just one project
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef CUT_ODDHOLES
#define  CUT_ODDHOLES
#include "build_cgraph.h"
#include "cut_pool.h"

struct _CutOH {

    VecInt **cutElem; /*index of cuts elements */
    VecDbl **cutCoef; /* coefficient of cuts elements */
    VecDbl *cutrhs; /* rhs of cuts */
    VecDbl *cutviolation; /* violation of cuts */
    VecInt *cutsense; /* sense of cuts */
    VecInt *cutnelem; /* size of cuts */
    VecInt *cutdominated; /* if cut is dominated or not (1,0)*/
    VecStr *cutname; /* name of cuts */
    int nAllocations;
    const Instance* inst;
};

typedef struct _CutOH CutOH;

void CutOH_add_cuts_conflicts_odd_parallel( const CGraph *cgraph, CutOH *coh,  LinearProgram *lp, const Instance *inst, double timeLeft,  double maxcuts, int nround );
CutOH *CutOH_create(const Instance *inst);
void CutOH_free( CutOH **cutOH );

#endif // CUT_ODDHOLES

