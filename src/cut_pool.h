
/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * mip_compact: routines for the creation of compact MIP formulations
 *              for the whole problem or for just one project
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef CUT_POOL
#define  CUT_POOL

#include "macros.h"
#include "vec_str.h"
#include "vec_int.h"
#include "vec_double.h"
#include "dict_int.h"
#include "tokenizer.h"
#include "vec_str.h"
#include "vint_set.h"
#include "str_utils.h"
#include "cgraph.h"
#include "lp.h"
#include "results.h"
#include "parameters.h"

#include <limits.h>
#include <omp.h>
#include <assert.h>
#include <string.h>

typedef struct _CutPool CutPool;

typedef struct {
    int j;
    int m;
    int t;
    float value;
    int idx;
} IntTriple;
/*
typedef struct {
    int a;
    int b;
} IntPair;
*/
typedef struct {
    int a;
    double b;
} IntDblPair;


typedef struct {
    int a;
    double b;
    double c;
} IntDblTuple;


typedef struct {
    int *idx;
    double *coef;
    double rhs;
    int sense;
    int sizec;
    int row;
    int endc;
    int ncut;

} IdxCoefCut;

int parseName( const char *name, char *prefix, int *idx );
int binary_search_cut (int elem, const int *idxB,  int l, int r);
void CutP_setIdRow( CutPool *cutP, int key, int idRow, int value);
IdxCoefCut *CutP_getCut(CutPool *cutP, int key, int idxR);
void CutP_setNameRow( CutPool *cutP, int key, int idCut, const char *value);
/*Creates the hash to keep all sets of cuts.
For each cut in hashCuts the following vector positions indicate:
- position 0 indicates whether or not the cut is active (1 active, -1 inactive).
- position 1 type of cut (LPC_RR, LPC_PREC, LPC_CLIQUE, LPC_CG, LPC_ODDHOLES, LPC_CGGPUR2, LPC_CGGPU, LPC_CGCPU,LPC_JS)
- position 2 value of the rhs of cut
- position 3 the sense of cut
- position 4 the number of elements of cut
- from position 5 the index of the elements first and then the coefficient.*/
CutPool *CutP_create( const Instance *inst,  LinearProgram *mip);

/* Separates the enabled cuts */
int CutP_separation(CutPool *cutP, Parameters *par, Results *res, int cr, int cp, int cc, int co, int cg, int cggpu, int cggpur2,  int** maxTJM, int* maxTJ, double startT, double timeLeft);
int CutP_separationParallel(CutPool *cutP, FILE *fp2, const Parameters *par,  Results *res, int cr, int cp, int cc, int co, int cg, int cggpur2, int cgcpu, int** maxTJM, int* maxTJ, double startT, double timeLeft);
int CutP_separationParallelCgraph(CutPool *cutP, FILE *fp2, CGraph *cgraph, const Parameters *par,  Results *res, int cr, int cp, int cc, int co, int cg, int cggpur2, int cgcpu,   int** maxTJM, int* maxTJ, double startT, double timeLeft);

IntDblTuple cutP_findElement( CutPool *cutP, int key, int c, const int *idx, const double *coe, double rhs, int sense, int lp_cols, int type);
IntDblTuple CutP_compElemRepeatedCoef( CutPool *cutP, int key, int c, const int* idx, double* coe, double *rhs, int sense, int lp_cols, int lp_rows, int type); // cutPool

void CutP_removeCut( CutPool *cutP, LinearProgram* lp, double maxslack);
void CutP_addCut( CutPool *cutP, LinearProgram* lp, int continuous, double maxslack, int nround);
int CutP_addSeparatedCuts( CutPool *cutP, Results *res, LinearProgram* lp,  VecInt **cutsElem, VecDbl **cutsCoef, VecDbl *cutsRhs, VecDbl *cutsViolation, VecStr *cutsName, VecInt *cutsSense, VecInt *cutDominated, int continuous, int nround, int cut, double timerem); // double timeseparation, int naddcut)
CGraph *CutP_compute_conflicts_create_graph( LinearProgram *lp, const int *origCols, LinearProgram *origLP, const Instance *inst, double timeLeft);
CGraph *CutP_compute_conflicts_create_complete_graph( LinearProgram *lp, const int *origCols, LinearProgram *origLP, const Instance *inst, double timeLeft);


void CutP_printHash( CutPool *cutP);
VecDbl** CutP_getHC( CutPool *cutP );
int* CutP_getnHC( CutPool *cutP );
void CutP_free( CutPool **cutP );
double CutP_getTotalTimeSeparation( CutPool *cutP );

int hash_code_vint( int n, const int *v, int hashSize );

int CutP_dominanceBetweenTwo( LinearProgram *lp, int idxcutA,  int *idxA, double *coefA, double rhs, int sense, int sizeA,  int idxcutB,  int *idxB, double *coefB, double rhsB, int senseB, int sizeB, VecInt *dominatedcuts);
int CutP_dominance(CutPool *cutP, const int *idxA, double *coefA, double rhs, int sense, int typeA, int sizeA);
int CutP_maxDivisorCommonRec(int m, int n);
int CutP_maxDivisorCommonVector(int coefs[], int nElem);
void CutP_quick_sort(IntDblPair *a, int n);
void CutP_quick_sort_vec ( int *idx, double *coe, int n);
void CutP_quick_sort_vec_by_double ( int *idx, double *coe, int n);

double CutP_model_lift_cgraph( const Instance *inst, const CGraph *cgraph, const char* cutname, int cutnelem, int *cutElem, double *cutCoef, double rhs, LinearProgram *lp, double timeLeft, int horizon);
double CutP_model_lift( const Instance *inst, const char* cutname, int cutnelem, int *cutElem, double *cutCoef, double rhs, LinearProgram *lp, double timeLeft, int horizon);

#endif // CUT_POOL

