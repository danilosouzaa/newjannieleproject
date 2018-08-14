
/* GOAL PSP Solver II
 * Resource-Constrained Project Scheduling Problems (RCPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G.
 */

#ifndef CUT_CG
#define  CUT_CG
#include "clique_separation.h"
#include "bron_kerbosch.h"
#include "cut_clique.h"
#include "cut_pool.h"
#include "build_cgraph.h"


struct _CutCG {

    VecInt *idxElemPP; /*index preprocessed of elements*/
    VecDbl *xfElemPP;  /*xf values of elements with index preprocessed*/
    VecDbl *rowrhs; /*rhs of selected rows*/
    VecInt *rowtype; /*type of selected rows*/
    VecInt *rowtimeslotmin; /*started of time slot  of selected rows*/
    VecInt *rowtimeslotmax; /*end of time slot  of selected rows*/
    VecInt **rowElem;
    VecDbl **rowCoef;
    VecInt **cutElem;
    VecDbl **cutCoef;
    int nAllocations;
    VecDbl *cutrhs;
    VecDbl *cutviolation;
    VecInt *cutnelem;
    VecInt *cutsense;
    VecInt *cutdominated; /* if cut is dominated or not (1,0)*/
    VecStr *cutname;


    // int nAddCuts;
    //double timeseparation;

    /*to cg original*/
    int cont;
    int contElem;
    int nCols;

    int *elements;
    int *idxelements;

    VecInt *idxelementsvec; /*original index of elements*/
    VecInt **elemrow;
    VecDbl **coefelemrows;
    VecStr *rname;
    VecStr *nameelements;
    VecDbl *rrhs;

    int contUR, contURN, contUM, contUCI, contUCP, nJumps;

    const Instance* inst;
};

typedef struct _CutCG CutCG;

CutCG *CutCG_create(const Instance *inst, int nRow, int nElem, int nCols);
void CutCG_add_cuts_model_CG_parallel( CutCG *ccg, const CGraph *cgraph, int cutType, LinearProgram *lp, const Instance *inst, double timeLeft,   double maxcuts, int nround, int horizon);
//void CutCG_add_cuts_model_CG_parallel( CutCG *ccg,int cutType, LinearProgram *lp, const Instance *inst, double timeLeft,   double maxcuts, int nround, int horizon);
//void CutCG_add_cuts_model_CG_parallel( CutCG *ccg,int cutType,LinearProgram *lp, const Instance *inst, double timeLeft, int continuous,  double maxcuts, int nround);
LinearProgram * lp_create_cgsep(const double *xf, int nr, int nelem, int *idxelements, VecStr *nameelements, VecInt **elemrow, VecDbl **coefelemrows, VecDbl *rrhs, VecStr *rname, double timeleft);

CutCG *CutCG_creat_and_identify_rows( CGraph *cgraph, LinearProgram *lp, const int *origCols, LinearProgram *origLP, const Instance *inst, double timeLeft, int continuous,  int **maxTJM,  double maxcuts, int mininstant, int maxinstant, int jump, int nround);
CutCG *CutCG_create_and_identify_rows(  LinearProgram *lp, const int *origCols, LinearProgram *origLP, const Instance *inst, double timeLeft, int continuous,  int **maxTJM,  double maxcuts, int mininstant, int maxinstant, int jump, int nround);
CutCG *CutCG_create_and_identify_rows_all(  LinearProgram *lp, const int *origCols, LinearProgram *origLP, const Instance *inst, double timeLeft,   int mininstant, int maxinstant, int jump, int nround);

//CutCG *CutCG_create_and_identify_rows_all_conflitcts(  LinearProgram *lp, const int *origCols, LinearProgram *origLP, const Instance *inst, double timeLeft, int continuous,  int **maxTJM,  double maxcuts, int mininstant, int maxinstant, int jump, int nround);
CutCG *CutCG_create_and_identify_rows_all_conflitcts(  LinearProgram *lp, const int *origCols, LinearProgram *origLP, const Instance *inst, double timeLeft, int mininstant, int maxinstant, int jump, int nround, int *horizon);
CutCG *CutCG_create_and_identify_rows_all_conflitcts_cgraph(  LinearProgram *lp, const CGraph *cgraph, const int *origCols, LinearProgram *origLP, const Instance *inst, double timeLeft, int mininstant, int maxinstant, int jump, int nround, int *horizon);

//void CutCG_add_cuts_CG_parallel( CutCG *ccg,int cutType,LinearProgram *lp, const Instance *inst, double timeLeft, int continuous,  double maxcuts, int nround);
void CutCG_add_cuts_CG_parallel( CutCG *ccg, const CGraph *cgraph, int cutType,LinearProgram *lp, const Instance *inst, double timeLeft, int continuous,  double maxcuts, int nround, int horizon);
//void CutCG_add_cuts_CG_parallel( CutCG *ccg,int cutType,LinearProgram *lp, const Instance *inst, double timeLeft, int continuous,  double maxcuts, int nround, int horizon);

void CutCG_print(CutCG *ccg);
void CutCG_printCut(CutCG *ccg);

VecDbl *CutCG_getrowrhs(CutCG *ccg);
VecDbl *CutCG_getxfElemPP(CutCG *ccg);

void CutCG_free( CutCG **cutCG );

#endif // CUT_CG

