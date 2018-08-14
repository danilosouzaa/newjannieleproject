/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * mip_compact: routines for the creation of compact MIP formulations
 *              for the whole problem or for just one project
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */
#include "cut_oddholes.h"
#include "oddhs.h"

#define VERBOSE  1
#define PERCMAXCLIQUE  0.5

#define MIN_VIOL 0.02


CutOH *CutOH_create(const Instance *inst)
{

    CutOH *coh;
    ALLOCATE_INI(coh,CutOH);

    coh->cutrhs = VDbl_create();
    coh->cutviolation = VDbl_create();
    coh->cutnelem = VInt_create();
    coh->cutsense = VInt_create();
    coh->cutname = VStr_create(256);
    coh->cutdominated = VInt_create();
    coh->nAllocations = 0;

    coh->inst = inst;

    return coh;
}


void CutOH_free( CutOH **_cutOH )
{

    CutOH *cutOH = *_cutOH;

    for(int nr = 0 ; nr < cutOH->nAllocations ; nr++) {
        VInt_free(&cutOH->cutElem[nr]);
        VDbl_free(&cutOH->cutCoef[nr]);
    }

    free(cutOH->cutElem);
    free(cutOH->cutCoef);

    VStr_free(&cutOH->cutname);
    VDbl_free(&cutOH->cutrhs);
    VInt_free(&cutOH->cutnelem);
    VInt_free(&cutOH->cutsense);
    VInt_free(&cutOH->cutdominated);
    VDbl_free(&cutOH->cutviolation);

    free( cutOH );
    *_cutOH = NULL;

}

void CutOH_add_cuts_conflicts_odd_parallel( const CGraph *cgraph, CutOH *coh,  LinearProgram *lp, const Instance *inst, double timeLeft,  double maxcuts, int nround )
{

    //    double tinit = omp_get_wtime();
    // double iniOH = omp_get_wtime();

    int ncols = lp_cols(lp);
    int nCut =  lp_rows(lp);

    OddHoleSep *oddhs = oddhs_create();
    oddhs_search_odd_holes( oddhs, ncols, lp_x(lp), lp_reduced_cost(lp), cgraph );
    int j = 0, i=0;
    // int nRows = lp_rows(lp)*PERCMAXCLIQUE;
    // if(nRows<maxcuts) nRows = maxcuts;
    int nodds = oddhs_get_odd_hole_count(oddhs);


    // double timeOH = omp_get_wtime()-iniOH;
    // printf("timeOH %f \n", timeOH); fflush(stdout);

    if(nodds>0) {
        //REALLOCATE_VECTOR(coh->cutElem,VecInt*,nodds);
        //REALLOCATE_VECTOR(coh->cutCoef,VecDbl*,nodds);

        ALLOCATE_VECTOR(coh->cutElem,VecInt*,nodds);
        ALLOCATE_VECTOR(coh->cutCoef,VecDbl*,nodds);

        while(j < nodds) {

            const int *oddEl = oddhs_get_odd_hole( oddhs, j );
            const int size = oddhs_get_odd_hole( oddhs, j+1 ) - oddEl;
            const double rhs = oddhs_rhs( size );
            const double vi = oddhs_viol( size, oddEl, lp_x(lp));
            int cutSize = 0;

            if ( vi < MIN_VIOL ) {
                j++;
                continue;
            }

            // if(i>0) {
            coh->cutElem[i] = VInt_create();
            coh->cutCoef[i] = VDbl_create();
            coh->nAllocations++;
            //}

            for(int sel = 0 ; sel<size ; sel++) {
                VInt_pushBack(coh->cutElem[i],oddEl[sel]);
                VDbl_pushBack(coh->cutCoef[i],1.0);
                //oddIdx[cutSize] = oddEl[j];
                //oddCoef[cutSize] = 1.0;
                cutSize++;
            }



            const int centerSize = oddhs_get_nwc_doh( oddhs, j );
            const int *centerIdx = oddhs_get_wc_doh( oddhs, j );

            for(int sel = 0; sel < centerSize; sel++) {
                VInt_pushBack(coh->cutElem[i],centerIdx[sel]);
                VDbl_pushBack(coh->cutCoef[i], rhs);
                //oddIdx[cutSize] = centerIdx[j];
                //oddCoef[cutSize] = rhs;
                cutSize++;
            }



            char nameCutodd[STR_SIZE];
            sprintf( nameCutodd, "cutODDHOLES(%d)#%d#%d", i,nround,lp_rows(lp)+nCut++);
            // printf(  "cutODDHOLES(%d)#%d#%d", i,nround,lp_rows(lp)+nCut++);
            int *cutIdx = VInt_getPtr(coh->cutElem[i]);
            double *cutCoef = VDbl_getPtr(coh->cutCoef[i]);
            VStr_pushBack(coh->cutname, nameCutodd);
            VDbl_pushBack(coh->cutrhs, rhs);
            VInt_pushBack(coh->cutsense, 0);
            VInt_pushBack(coh->cutdominated,0);
            VInt_pushBack(coh->cutnelem, cutSize);
            VDbl_pushBack(coh->cutviolation, vi);

            // printf(" rhs %f\n vi %f", VDbl_get(coh->cutrhs,i), vi); getchar();

            CutP_quick_sort_vec(cutIdx,cutCoef, cutSize);

            j++;
            i++;
        }
    }

    oddhs_free( &oddhs );



    /*  double timeinitdomresource = omp_get_wtime();
      for(int ncdA = 0; ncdA < i ; ncdA++){
          for(int ncdB = ncdA+1; ncdB < i ; ncdB++){
              if(VInt_get(coh->cutdominated,ncdA) == 1) break;
              if(VInt_get(coh->cutdominated,ncdB) == 1) continue;
              CutP_dominanceBetweenTwo( lp, ncdA, VInt_getPtr(coh->cutElem[ncdA]), VDbl_getPtr(coh->cutCoef[ncdA]), VDbl_get(coh->cutrhs,ncdA), VInt_get(coh->cutsense,ncdA), VInt_get(coh->cutnelem,ncdA), ncdB, VInt_getPtr(coh->cutElem[ncdB]), VDbl_getPtr(coh->cutCoef[ncdB]), VDbl_get(coh->cutrhs,ncdB), VInt_get(coh->cutsense,ncdB), VInt_get(coh->cutnelem,ncdB),  coh->cutdominated);
              //   getchar();
          }
      }
      printf("time compute dominance odd holes %f\n",omp_get_wtime()-timeinitdomresource);
    */


#undef MAX_IDX
}

