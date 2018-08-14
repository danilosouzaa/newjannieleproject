
/* GOAL PSP Solver II
 * Resource-Constrained Project Scheduling Problems (RCPSP)
 *
 * Develop as part of the D.Sc. thesis of Araujo, Janniele A. S., with collaboration
 *                                   of Santos, H.G.
 */


#ifndef CUT_COVER
#define  CUT_COVER
#include "cut_pool.h"

struct _CutC {

    VecInt ***cutElem; /*index of cuts elements */
    VecDbl ***cutCoef; /* coefficient of cuts elements */
    VecDbl **cutrhs; /* rhs of cuts */
    VecDbl **cutviolation; /* violation of cuts */
    VecInt **cutsense; /* sense of cuts */
    VecInt **cutnelem; /* size of cuts */
    VecStr **cutname; /* name of cuts */
    VecInt **cutdominated; /* if cut is dominated or not (1,0)*/
    int *nAllocations;

    const Instance* inst;
};

typedef struct _CutC CutC;

CutC *CutC_create(const Instance *inst);

void CutC_add_cuts_cover_parallel( CutC *cc, LinearProgram *lp, const int *origCols, LinearProgram *origLP, const Instance *inst, double timeLeft, int continuous, int lifting, int nround);
void CutC_add_cuts_cover_model_parallel( CutC *ccr, LinearProgram *lp, const int *origCols,  LinearProgram *origLP, const Instance *inst, double timeLeft, int continuous, int lifting,  int nround );

void CutC_free( CutC **cutC );

#endif // CUT_COVER

