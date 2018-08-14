/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <assert.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <omp.h>
#include "vnd.h"
#include "vns.h"
#include "neighborhood.h"
#include "ms_solver_mip.h"
#include "instance.h"
#include "memory.h"
#include "macros.h"

#define N_NEIGHBORHOOD_LS 14
#define BETA 0

struct _Neighborhood {

    int nNeighborhood;
    long double *intensity;
    long double *intensityAux;
    int *assortment;
    int *neighborMinK;
    int *neighborMaxK;
    int it;
    int penaltyChangeMode;
    int lastN;
    int timePrint;
    int *nLastJModify;
    int *posNLastJModify;
    int contLastJ;
    int moveSW;
    int log;
    int sw;
    int online;
    int learning;
    int uniform;
    int typeinstance;

    int *lastJ;
    int *newModes;

    long double *I;
    long double *EQ;
    long double *TI;
    long double *TE;
    long double *TIV;
    long double *FI;
    long double *FE;

    long double *normFI;
    long double *normFE;

    long double *PI;
    long double *PE;

    long double maxFI;
    long double minFI;

    long double maxFE;
    long double minFE;

    long double intervalI;
    long double intervalE;

    //MSM_Solver **mss;

    int nThread;
    int nStages;

    const struct _Instance *inst;
};

void Neighbor_checkArgs(Neighborhood *neighbor, char **argv, int argc)
{

    assert( neighbor!=NULL );

    double intens;
    int assort;
    for(int n=0 ; n < argc; n++) {

        if (strcmp(argv[n],"-intensity") == 0) {
            n++;
            int idx = atoi(argv[n]) -1;
            n++;
            intens = atof( argv[n] );
            neighbor->intensity[idx] = intens;
            continue;
        }
        if (strcmp(argv[n],"-intensityAux") == 0) {
            n++;
            int idx = atoi(argv[n]) -1;
            n++;
            intens = atof( argv[n] );
            neighbor->intensityAux[idx] = intens;
            continue;
        }
        if (strcmp(argv[n],"-assort") == 0) {
            n++;
            int idx = atoi(argv[n]);
            n++;
            assort = atoi( argv[n] );
            neighbor->assortment[idx] = assort;
            continue;
        }
        if (strcmp(argv[n],"-minK") == 0) {
            n++;
            int idx = atoi(argv[n]);
            n++;
            neighbor->neighborMinK[idx] = atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-maxK") == 0) {
            n++;
            int idx = atoi(argv[n]);
            n++;
            neighbor->neighborMaxK[idx] = atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-penaltyCM") == 0) {
            n++;
            neighbor->penaltyChangeMode = atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-nThread") == 0) {
            n++;
            neighbor->nThread = atoi( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-nStages") == 0) {
            n++;
            neighbor->nStages = atoi( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-mmrcmpsp") == 0) {
            n++;
            neighbor->typeinstance = 2;
            continue;
        }
        if (strcmp(argv[n],"-mrcpsp") == 0) {
            n++;
            neighbor->typeinstance = 1;
            continue;
        }
        if (strcmp(argv[n],"-rcpsp") == 0) {
            n++;
            neighbor->typeinstance = 0;
            continue;
        }
        if (strcmp(argv[n],"-log") == 0) {
            n++;
            neighbor->log = atoi( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-sw") == 0) {
            n++;
            neighbor->sw = atoi( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-online") == 0) {
            n++;
            neighbor->online = atoi( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-learning") == 0) {
            n++;
            neighbor->learning =  atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-uniform") == 0) {
            n++;
            neighbor->uniform =  atoi(argv[n]);
            continue;
        }

    }
}

Neighborhood *Neighbor_create( const Instance* inst, char ** argv, int argc )
{

    assert( inst!=NULL );

    Neighborhood *neighborhood;
    ALLOCATE_INI(neighborhood, Neighborhood);

    neighborhood->nNeighborhood = N_NEIGHBORHOOD_LS;

    long double *intensity;
    ALLOCATE_VECTOR_INI(intensity, long double, N_NEIGHBORHOOD_LS);
    long double *intensityAux;
    ALLOCATE_VECTOR_INI(intensityAux,long double, N_NEIGHBORHOOD_LS);
    int * assortment;
    ALLOCATE_VECTOR_INI(assortment,int, N_NEIGHBORHOOD_LS);
    int * lastJ;
    ALLOCATE_VECTOR_INI(lastJ,int, 4);
    int * newModes;
    ALLOCATE_VECTOR_INI(newModes,int, 4);
    int * nLastJModify;
    ALLOCATE_VECTOR_INI(nLastJModify, int, Inst_nJobs(inst)*Inst_nJobs(inst)*Inst_nJobs(inst));
    int * posNLastJModify;
    ALLOCATE_VECTOR_INI(posNLastJModify,int, Inst_nJobs(inst)*Inst_nJobs(inst)*Inst_nJobs(inst));


    /*  ALLOCATE_VECTOR(neighborhood->mss, MSM_Solver*, Inst_nProjects(inst));
      for(int i = 0; i < Inst_nProjects(inst) ; i++) {
          const Project *project = Inst_project( inst, i );
          neighborhood->mss[i] = MSM_create( inst, Project_idxFirstJob(project), Project_idxFirstJob(project)+Project_nJobs(project)-1 );
          MSM_solve( neighborhood->mss[i] );
      }
    */
    /*number of improvements*/
    long double *I = (long double *) calloc (neighborhood->nNeighborhood, sizeof(long double));
    neighborhood->I = I;

    /*number of sideway moves *0.1 (cookies) */
    long double *EQ = (long double *) calloc (neighborhood->nNeighborhood, sizeof(long double));
    neighborhood->EQ = EQ;

    /*total time spent with improvements*/
    long double *TI = (long double *) calloc (neighborhood->nNeighborhood, sizeof(long double));
    neighborhood->TI = TI;

    /*total time spent with a sideway moves*/
    long double *TE = (long double *) calloc (neighborhood->nNeighborhood, sizeof(long double));
    neighborhood->TE = TE;

    /*total time spent without improve solution*/
    long double *TIV = (long double *) calloc (neighborhood->nNeighborhood, sizeof(long double));
    neighborhood->TIV = TIV;


    /*value of improvements by time*/
    long double *FI= (long double *) calloc (neighborhood->nNeighborhood, sizeof(long double));
    neighborhood->FI = FI;
    /*value of improvements plus sideway moves by time*/
    long double *FE= (long double *) calloc (neighborhood->nNeighborhood, sizeof(long double));
    neighborhood->FE = FE;

    /*values of FI normalized*/
    long double *normFI= (long double *) calloc (neighborhood->nNeighborhood, sizeof(long double));
    neighborhood->normFI = normFI;
    /*values of FE normalized*/
    long double *normFE= (long double *) calloc (neighborhood->nNeighborhood, sizeof(long double));
    neighborhood->normFE = normFE;

    /*variables to calculate the interval to normalize*/
    neighborhood->minFI = (long double)INT_MAX_M;
    neighborhood->minFE = (long double)INT_MAX_M;
    neighborhood->maxFI = 0;
    neighborhood->maxFE = 0;
    neighborhood->timePrint = 0;
    neighborhood->log =0;
    neighborhood->online =0;
    neighborhood->learning =0;
    neighborhood->uniform =0;

    /*0 to identify fake sideways moves and
      1 to identify real sideways moves
      of the LAST neighborhood applied*/
    neighborhood->moveSW = 0;

    assortment[0] = seqSwapProj;
    assortment[1] = seqCompOnExtrem;
    assortment[2] = seqCompactProj;
    assortment[3] = seqShiftProj;
    assortment[4] = seqSwapJobFILS;
    assortment[5] = changeTwoMode;
    assortment[6] = changeOneMode;
    assortment[7] = changeThreeMode;
    assortment[8] = changeFourMode;
    assortment[9] = seqShiftJob;
    assortment[10] = seqInsertJobFILS;
    assortment[11] = seqSwapJob;
    assortment[12] = seqInvert;
    assortment[13] = seqMoveProj;

    /*sort: analysis 2 neighborhood*/
    /*assortment[0] = changeOneMode;
    assortment[1] = seqShiftProj;
    assortment[2] = changeTwoMode;
    assortment[3] = seqShiftJob;
    assortment[4] = seqCompactProj;
    assortment[5] = changeThreeMode;
    assortment[6] = seqSwapJob;
    assortment[7] = seqSwapProj;
    assortment[8] = changeFourMode;
    assortment[9] = seqCompOnExtrem;
    assortment[10] = seqInvert;
    assortment[11] = seqSwapJobFILS;
    assortment[12] = seqInsertJobFILS;
    assortment[13] = seqMoveProj;
    */
    int *minK, *maxK;
    ALLOCATE_VECTOR_INI(minK,int, N_NEIGHBORHOOD_LS);
    ALLOCATE_VECTOR_INI(maxK,int, N_NEIGHBORHOOD_LS);

    minK[seqInvert-1] = 3.0;
    minK[seqShiftJob-1] = 1.0;
    minK[seqSwapJob-1] = 1.0;
    minK[seqShiftProj-1] = 1.0;
    minK[seqSwapProj-1] = 1.0;
    minK[seqCompactProj-1] = 1.0;
    minK[changeOneMode-1] = 1.0;
    minK[changeTwoMode-1] = 1.0;
    minK[changeThreeMode-1] = 1.0;
    minK[changeFourMode-1] = 1.0;
    minK[seqSwapJobFILS-1] = 1.0;
    minK[seqInsertJobFILS-1] = 1.0;
    minK[seqCompOnExtrem-1] = 1.0;
    minK[seqMoveProj-1] = 1.0;

    maxK[seqInvert-1] = 10;
    maxK[seqShiftJob-1] = 10;
    maxK[seqSwapJob-1] = 10;
    maxK[seqShiftProj-1] = 10;
    maxK[seqSwapProj-1] = 10;
    maxK[seqCompactProj-1] = 10;
    maxK[changeOneMode-1] = 10;
    maxK[changeTwoMode-1] = 10;
    maxK[changeThreeMode-1] = 10;
    maxK[changeFourMode-1] = 10;
    maxK[seqSwapJobFILS-1] = 10;
    maxK[seqInsertJobFILS-1] = 10;
    maxK[seqCompOnExtrem-1] = 10;
    maxK[seqMoveProj-1] = 10;

    neighborhood->neighborMinK = minK;
    neighborhood->neighborMaxK = maxK;

    neighborhood->intensity = intensity;
    neighborhood->intensityAux = intensityAux;
    neighborhood->assortment = assortment;
    neighborhood->it = 0;

    neighborhood->lastN = 0;
    neighborhood->penaltyChangeMode = 0;
    neighborhood->lastJ = lastJ;
    neighborhood->newModes = newModes;
    neighborhood->nLastJModify = nLastJModify;
    neighborhood->posNLastJModify = posNLastJModify;

    neighborhood->nThread = 1;
    neighborhood->nStages = 2;
    neighborhood->contLastJ = 0;
    neighborhood->sw = 0;
    neighborhood->inst = inst;

    Neighbor_checkArgs(neighborhood, argv,argc);

    if(neighborhood->online || neighborhood->learning || neighborhood->uniform) {
        /* padrao to online and learning*/
        intensity[seqInvert-1] = 1.000;
        intensity[seqShiftJob-1] = 1.000;
        intensity[seqSwapJob-1] =1.000;
        intensity[seqShiftProj-1] =1.000;
        intensity[seqSwapProj-1] = 1.000;
        intensity[seqCompactProj-1] = 1.000;
        intensity[changeOneMode-1] = 1.000;
        intensity[changeTwoMode-1] = 1.000;
        intensity[changeThreeMode-1] = 1.000;
        intensity[changeFourMode-1] = 1.000;
        intensity[seqSwapJobFILS-1] = 1.000;
        intensity[seqInsertJobFILS-1] = 1.000;
        intensity[seqCompOnExtrem-1] = 1.000;
        intensity[seqMoveProj-1] = 1.000;
    } else {
        /*without normalized time*/
        if(neighborhood->sw) {
            /* intensities default to two stage inside LAHC*/
            intensity[seqInvert-1] = 0.114898;
            intensity[seqShiftJob-1] = 0.239016;
            intensity[seqSwapJob-1] = 0.178385;
            intensity[seqShiftProj-1] = 0.574908;
            intensity[seqSwapProj-1] = 0.994617;
            intensity[seqCompactProj-1] =  0.949450;
            intensity[changeOneMode-1] = 0.373147;
            intensity[changeTwoMode-1] = 0.337556 ;
            intensity[changeThreeMode-1] = 0.282499;
            intensity[changeFourMode-1] = 0.228154;
            intensity[seqSwapJobFILS-1] = 0.374059;
            intensity[seqInsertJobFILS-1] = 0.122834;
            intensity[seqCompOnExtrem-1] = 1.000000;
            intensity[seqMoveProj-1] = 0.262191;

            intensityAux[seqInvert-1] = 0.228947;
            intensityAux[seqShiftJob-1] = 0.436167;
            intensityAux[seqSwapJob-1] = 0.320362;
            intensityAux[seqShiftProj-1] = 0.270712;
            intensityAux[seqSwapProj-1] = 0.075088;
            intensityAux[seqCompactProj-1] =  0.247841;
            intensityAux[changeOneMode-1] = 0.192918;
            intensityAux[changeTwoMode-1] =  0.042769;
            intensityAux[changeThreeMode-1] = 0.026710;
            intensityAux[changeFourMode-1] = 0.002802;
            intensityAux[seqSwapJobFILS-1] = 0.021347;
            intensityAux[seqInsertJobFILS-1] = 0.005274;
            intensityAux[seqCompOnExtrem-1] = 0.252793;
            intensityAux[seqMoveProj-1] = 1.000000;

        } else {
            intensity[seqInvert-1] = 0.079260;
            intensity[seqShiftJob-1] = 0.148346;
            intensity[seqSwapJob-1] = 0.118582;
            intensity[seqShiftProj-1] = 0.539655;
            intensity[seqSwapProj-1] = 1.000000;
            intensity[seqCompactProj-1] = 0.934844;
            intensity[changeOneMode-1] = 0.314266;
            intensity[changeTwoMode-1] = 0.323262;
            intensity[changeThreeMode-1] = 0.274779;
            intensity[changeFourMode-1] = 0.228228;
            intensity[seqSwapJobFILS-1] = 0.378151;
            intensity[seqInsertJobFILS-1] = 0.124178;
            intensity[seqCompOnExtrem-1] =0.987266 ;
            intensity[seqMoveProj-1] = 0.000000;

            intensityAux[seqInvert-1] = 0.226802;
            intensityAux[seqShiftJob-1] =0.226771 ;
            intensityAux[seqSwapJob-1] = 0.211450;
            intensityAux[seqShiftProj-1] = 0.277710;
            intensityAux[seqSwapProj-1] = 0.039467;
            intensityAux[seqCompactProj-1] = 0.006135;
            intensityAux[changeOneMode-1] = 0.658818;
            intensityAux[changeTwoMode-1] = 0.273215;
            intensityAux[changeThreeMode-1] = 0.137227;
            intensityAux[changeFourMode-1] = 0.031540;
            intensityAux[seqSwapJobFILS-1] = 1.000000;
            intensityAux[seqInsertJobFILS-1] = 0.247072;
            intensityAux[seqCompOnExtrem-1] = 0.147583;
            intensityAux[seqMoveProj-1] = 0.000000;
        }
    }

    /*converting intensities at probabilities of Intensity normalized*/
    long double *PI = (long double *) calloc (neighborhood->nNeighborhood, sizeof(long double));
    long double *PE = (long double *) calloc (neighborhood->nNeighborhood, sizeof(long double));

    long double sumIntensities = 0;
    for(int i=0; i< neighborhood->nNeighborhood; ++i) {
        sumIntensities += intensity[i];// +BETA;
    }

    PE[seqInvert-1] = (intensity[seqInvert-1])/sumIntensities;
    PE[seqShiftJob-1] = (intensity[seqShiftJob-1])/sumIntensities;
    PE[seqSwapJob-1] = (intensity[seqSwapJob-1])/sumIntensities;
    PE[seqShiftProj-1] =(intensity[seqShiftProj-1])/sumIntensities;
    PE[seqSwapProj-1] = (intensity[seqSwapProj-1])/sumIntensities;
    PE[seqCompactProj-1] = (intensity[seqCompactProj-1])/sumIntensities;
    PE[changeOneMode-1] = (intensity[changeOneMode-1])/sumIntensities;
    PE[changeTwoMode-1] = (intensity[changeTwoMode-1])/sumIntensities;
    PE[changeThreeMode-1] = (intensity[changeThreeMode-1])/sumIntensities;
    PE[changeFourMode-1] = (intensity[changeFourMode-1])/sumIntensities;
    PE[seqSwapJobFILS-1] = (intensity[seqSwapJobFILS-1])/sumIntensities;
    PE[seqInsertJobFILS-1] = (intensity[seqInsertJobFILS-1])/sumIntensities;
    PE[seqCompOnExtrem-1] = (intensity[seqCompOnExtrem-1])/sumIntensities;
    PE[seqMoveProj-1] = (intensity[seqMoveProj-1])/sumIntensities;

    PI[seqInvert-1] = (intensity[seqInvert-1])/sumIntensities;
    PI[seqShiftJob-1] = (intensity[seqShiftJob-1])/sumIntensities;
    PI[seqSwapJob-1] =(intensity[seqSwapJob-1])/sumIntensities;
    PI[seqShiftProj-1] =(intensity[seqShiftProj-1])/sumIntensities;
    PI[seqSwapProj-1] = (intensity[seqSwapProj-1])/sumIntensities;
    PI[seqCompactProj-1] = (intensity[seqCompactProj-1])/sumIntensities;
    PI[changeOneMode-1] = (intensity[changeOneMode-1])/sumIntensities;
    PI[changeTwoMode-1] = (intensity[changeTwoMode-1])/sumIntensities;
    PI[changeThreeMode-1] = (intensity[changeThreeMode-1])/sumIntensities;
    PI[changeFourMode-1] = (intensity[changeFourMode-1])/sumIntensities;
    PI[seqSwapJobFILS-1] = (intensity[seqSwapJobFILS-1])/sumIntensities;
    PI[seqInsertJobFILS-1] = (intensity[seqInsertJobFILS-1])/sumIntensities;
    PI[seqCompOnExtrem-1] = (intensity[seqCompOnExtrem-1])/sumIntensities;
    PI[seqMoveProj-1] = (intensity[seqMoveProj-1])/sumIntensities;


    /*    PE[seqInvert-1] = (intensity[seqInvert-1]+BETA)/sumIntensities;
        PE[seqShiftJob-1] = (intensity[seqShiftJob-1]+BETA)/sumIntensities;
        PE[seqSwapJob-1] = (intensity[seqSwapJob-1]+BETA)/sumIntensities;
        PE[seqShiftProj-1] =(intensity[seqShiftProj-1]+BETA)/sumIntensities;
        PE[seqSwapProj-1] = (intensity[seqSwapProj-1]+BETA)/sumIntensities;
        PE[seqCompactProj-1] = (intensity[seqCompactProj-1]+BETA)/sumIntensities;
        PE[changeOneMode-1] = (intensity[changeOneMode-1]+BETA)/sumIntensities;
        PE[changeTwoMode-1] = (intensity[changeTwoMode-1]+BETA)/sumIntensities;
        PE[changeThreeMode-1] = (intensity[changeThreeMode-1]+BETA)/sumIntensities;
        PE[changeFourMode-1] = (intensity[changeFourMode-1]+BETA)/sumIntensities;
        PE[seqSwapJobFILS-1] = (intensity[seqSwapJobFILS-1]+BETA)/sumIntensities;
        PE[seqInsertJobFILS-1] = (intensity[seqInsertJobFILS-1]+BETA)/sumIntensities;
        PE[seqCompOnExtrem-1] = (intensity[seqCompOnExtrem-1]+BETA)/sumIntensities;
        PE[seqMoveProj-1] = (intensity[seqMoveProj-1]+BETA)/sumIntensities;

        PI[seqInvert-1] = (intensity[seqInvert-1]+BETA)/sumIntensities;
        PI[seqShiftJob-1] = (intensity[seqShiftJob-1]+BETA)/sumIntensities;
        PI[seqSwapJob-1] =(intensity[seqSwapJob-1]+BETA)/sumIntensities;
        PI[seqShiftProj-1] =(intensity[seqShiftProj-1]+BETA)/sumIntensities;
        PI[seqSwapProj-1] = (intensity[seqSwapProj-1]+BETA)/sumIntensities;
        PI[seqCompactProj-1] = (intensity[seqCompactProj-1]+BETA)/sumIntensities;
        PI[changeOneMode-1] = (intensity[changeOneMode-1]+BETA)/sumIntensities;
        PI[changeTwoMode-1] = (intensity[changeTwoMode-1]+BETA)/sumIntensities;
        PI[changeThreeMode-1] = (intensity[changeThreeMode-1]+BETA)/sumIntensities;
        PI[changeFourMode-1] = (intensity[changeFourMode-1]+BETA)/sumIntensities;
        PI[seqSwapJobFILS-1] = (intensity[seqSwapJobFILS-1]+BETA)/sumIntensities;
        PI[seqInsertJobFILS-1] = (intensity[seqInsertJobFILS-1]+BETA)/sumIntensities;
        PI[seqCompOnExtrem-1] = (intensity[seqCompOnExtrem-1]+BETA)/sumIntensities;
        PI[seqMoveProj-1] = (intensity[seqMoveProj-1]+BETA)/sumIntensities;
    */
    neighborhood->PE = PE;
    neighborhood->PI = PI;


    return neighborhood;
}

void Neighbor_free( Neighborhood **_neighbor )
{
    Neighborhood *neighbor = *_neighbor;

    free( neighbor->intensity );
    free( neighbor->intensityAux );
    free( neighbor->lastJ );
    free( neighbor->assortment );
    free( neighbor->neighborMaxK );
    free( neighbor->neighborMinK );
    free( neighbor->nLastJModify);
    free( neighbor->posNLastJModify);

    free( neighbor->I );
    free( neighbor->EQ );
    free( neighbor->TI );
    free( neighbor->TE );
    free( neighbor->TIV );
    free( neighbor->FI );
    free( neighbor->FE );
    free( neighbor->PE );
    free( neighbor->PI );

    free( neighbor->normFI );
    free( neighbor->normFE );

    /*  for(int i = 0; i < Inst_nProjects(neighbor->inst) ; i++)
          MSM_free(&neighbor->mss[i]);
      free(neighbor->mss);
    */
    free( neighbor );
    *_neighbor = NULL;
}

int Neighbor_getMoveSW(Neighborhood *neighborhood)
{

    assert(neighborhood != NULL);

    return neighborhood->moveSW;

}

void Neighbor_checkMoveSW(Neighborhood *neighborhood, Solution * sol, Solution * oldSol)
{
    assert(neighborhood != NULL);

    neighborhood->moveSW = 0;

    if(Sol_getCost(sol) == Sol_getCost(oldSol)) {
        for(int i = 0 ; i < Inst_nJobs(Sol_inst(sol)); i++) {
            //int posSolJob1 = Sol_getPosJob(sol,i);
            //int posOldSoljob1 = Sol_getPosJob(oldsol,i);
            if(Sol_getStartTime(sol,i) != Sol_getStartTime(oldSol,i)) {
                neighborhood->moveSW = 1;
                break;
            }
        }
    }

}

void Neighbor_getUpdatesIntensity(Neighborhood *neighborhood)
{
    assert(neighborhood != NULL);

    double sumIntensities = 0;

    for(int i = 0 ; i < neighborhood->nNeighborhood ; i++) {
        double value = neighborhood->intensity[i];
        neighborhood->intensity[i] = neighborhood->intensityAux[i];
        neighborhood->intensityAux[i] = value;
        sumIntensities += neighborhood->intensity[i];//+BETA;
    }

    neighborhood->PE[seqInvert-1] = (neighborhood->intensity[seqInvert-1])/sumIntensities;
    neighborhood->PE[seqShiftJob-1] = (neighborhood->intensity[seqShiftJob-1])/sumIntensities;
    neighborhood->PE[seqSwapJob-1] = (neighborhood->intensity[seqSwapJob-1])/sumIntensities;
    neighborhood->PE[seqShiftProj-1] =(neighborhood->intensity[seqShiftProj-1])/sumIntensities;
    neighborhood->PE[seqSwapProj-1] = (neighborhood->intensity[seqSwapProj-1])/sumIntensities;
    neighborhood->PE[seqCompactProj-1] = (neighborhood->intensity[seqCompactProj-1])/sumIntensities;
    neighborhood->PE[changeOneMode-1] = (neighborhood->intensity[changeOneMode-1])/sumIntensities;
    neighborhood->PE[changeTwoMode-1] = (neighborhood->intensity[changeTwoMode-1])/sumIntensities;
    neighborhood->PE[changeThreeMode-1] = (neighborhood->intensity[changeThreeMode-1])/sumIntensities;
    neighborhood->PE[changeFourMode-1] = (neighborhood->intensity[changeFourMode-1])/sumIntensities;
    neighborhood->PE[seqSwapJobFILS-1] = (neighborhood->intensity[seqSwapJobFILS-1])/sumIntensities;
    neighborhood->PE[seqInsertJobFILS-1] = (neighborhood->intensity[seqInsertJobFILS-1])/sumIntensities;
    neighborhood->PE[seqCompOnExtrem-1] = (neighborhood->intensity[seqCompOnExtrem-1])/sumIntensities;
    neighborhood->PE[seqMoveProj-1] = (neighborhood->intensity[seqMoveProj-1])/sumIntensities;

    neighborhood->PI[seqInvert-1] = (neighborhood->intensity[seqInvert-1])/sumIntensities;
    neighborhood->PI[seqShiftJob-1] = (neighborhood->intensity[seqShiftJob-1])/sumIntensities;
    neighborhood->PI[seqSwapJob-1] =(neighborhood->intensity[seqSwapJob-1])/sumIntensities;
    neighborhood->PI[seqShiftProj-1] =(neighborhood->intensity[seqShiftProj-1])/sumIntensities;
    neighborhood->PI[seqSwapProj-1] = (neighborhood->intensity[seqSwapProj-1])/sumIntensities;
    neighborhood->PI[seqCompactProj-1] = (neighborhood->intensity[seqCompactProj-1])/sumIntensities;
    neighborhood->PI[changeOneMode-1] = (neighborhood->intensity[changeOneMode-1])/sumIntensities;
    neighborhood->PI[changeTwoMode-1] = (neighborhood->intensity[changeTwoMode-1])/sumIntensities;
    neighborhood->PI[changeThreeMode-1] = (neighborhood->intensity[changeThreeMode-1])/sumIntensities;
    neighborhood->PI[changeFourMode-1] = (neighborhood->intensity[changeFourMode-1])/sumIntensities;
    neighborhood->PI[seqSwapJobFILS-1] = (neighborhood->intensity[seqSwapJobFILS-1])/sumIntensities;
    neighborhood->PI[seqInsertJobFILS-1] = (neighborhood->intensity[seqInsertJobFILS-1])/sumIntensities;
    neighborhood->PI[seqCompOnExtrem-1] = (neighborhood->intensity[seqCompOnExtrem-1])/sumIntensities;
    neighborhood->PI[seqMoveProj-1] = (neighborhood->intensity[seqMoveProj-1])/sumIntensities;

    /*
        neighborhood->PE[seqInvert-1] = (neighborhood->intensity[seqInvert-1]+BETA)/sumIntensities;
        neighborhood->PE[seqShiftJob-1] = (neighborhood->intensity[seqShiftJob-1]+BETA)/sumIntensities;
        neighborhood->PE[seqSwapJob-1] = (neighborhood->intensity[seqSwapJob-1]+BETA)/sumIntensities;
        neighborhood->PE[seqShiftProj-1] =(neighborhood->intensity[seqShiftProj-1]+BETA)/sumIntensities;
        neighborhood->PE[seqSwapProj-1] = (neighborhood->intensity[seqSwapProj-1]+BETA)/sumIntensities;
        neighborhood->PE[seqCompactProj-1] = (neighborhood->intensity[seqCompactProj-1]+BETA)/sumIntensities;
        neighborhood->PE[changeOneMode-1] = (neighborhood->intensity[changeOneMode-1]+BETA)/sumIntensities;
        neighborhood->PE[changeTwoMode-1] = (neighborhood->intensity[changeTwoMode-1]+BETA)/sumIntensities;
        neighborhood->PE[changeThreeMode-1] = (neighborhood->intensity[changeThreeMode-1]+BETA)/sumIntensities;
        neighborhood->PE[changeFourMode-1] = (neighborhood->intensity[changeFourMode-1]+BETA)/sumIntensities;
        neighborhood->PE[seqSwapJobFILS-1] = (neighborhood->intensity[seqSwapJobFILS-1]+BETA)/sumIntensities;
        neighborhood->PE[seqInsertJobFILS-1] = (neighborhood->intensity[seqInsertJobFILS-1]+BETA)/sumIntensities;
        neighborhood->PE[seqCompOnExtrem-1] = (neighborhood->intensity[seqCompOnExtrem-1]+BETA)/sumIntensities;
        neighborhood->PE[seqMoveProj-1] = (neighborhood->intensity[seqMoveProj-1]+BETA)/sumIntensities;

        neighborhood->PI[seqInvert-1] = (neighborhood->intensity[seqInvert-1]+BETA)/sumIntensities;
        neighborhood->PI[seqShiftJob-1] = (neighborhood->intensity[seqShiftJob-1]+BETA)/sumIntensities;
        neighborhood->PI[seqSwapJob-1] =(neighborhood->intensity[seqSwapJob-1]+BETA)/sumIntensities;
        neighborhood->PI[seqShiftProj-1] =(neighborhood->intensity[seqShiftProj-1]+BETA)/sumIntensities;
        neighborhood->PI[seqSwapProj-1] = (neighborhood->intensity[seqSwapProj-1]+BETA)/sumIntensities;
        neighborhood->PI[seqCompactProj-1] = (neighborhood->intensity[seqCompactProj-1]+BETA)/sumIntensities;
        neighborhood->PI[changeOneMode-1] = (neighborhood->intensity[changeOneMode-1]+BETA)/sumIntensities;
        neighborhood->PI[changeTwoMode-1] = (neighborhood->intensity[changeTwoMode-1]+BETA)/sumIntensities;
        neighborhood->PI[changeThreeMode-1] = (neighborhood->intensity[changeThreeMode-1]+BETA)/sumIntensities;
        neighborhood->PI[changeFourMode-1] = (neighborhood->intensity[changeFourMode-1]+BETA)/sumIntensities;
        neighborhood->PI[seqSwapJobFILS-1] = (neighborhood->intensity[seqSwapJobFILS-1]+BETA)/sumIntensities;
        neighborhood->PI[seqInsertJobFILS-1] = (neighborhood->intensity[seqInsertJobFILS-1]+BETA)/sumIntensities;
        neighborhood->PI[seqCompOnExtrem-1] = (neighborhood->intensity[seqCompOnExtrem-1]+BETA)/sumIntensities;
        neighborhood->PI[seqMoveProj-1] = (neighborhood->intensity[seqMoveProj-1]+BETA)/sumIntensities;
    */
}

void Neighbor_restartIntensity(Neighborhood *neighborhood)
{
    assert(neighborhood != NULL);

    for(int i = 0 ; i < Neighbor_nNeighborhood(neighborhood); i++) {
        neighborhood->intensity[i] = 1.000;
        neighborhood->PE[i] = 1.000/(neighborhood->nNeighborhood*1);
        neighborhood->PI[i] = 1.000/(neighborhood->nNeighborhood*1);
#ifdef DEBUG
        printf("\nRestart Intensities...");
#endif
    }


}

void Neighbor_setTime( Neighborhood *neighborhood)
{
    assert(neighborhood!=NULL);
    neighborhood->timePrint += 10;
}


int Neighbor_getTimePrint( Neighborhood *neighborhood)
{
    assert(neighborhood!=NULL);
    return neighborhood->timePrint;
}

void Neighbor_printProbabilitiesTime( Neighborhood *neighborhood, int sw, int it, float psw, char *argv, double timePrint)
{

    char fil[256] = "probabilitiesOnlineInstance_";
    char bufferII[256];
    char side[256];
    char penalty[256];
    char instance[256];
    sprintf(instance, "%s", argv);
    strcat (fil,instance);
    strcat (fil,"_");
    sprintf(bufferII, "%d", it);
    strcat (fil,bufferII);
    strcat (fil,"_");
    sprintf(side, "%d", sw);
    strcat (fil,side);
    strcat (fil,"_");
    sprintf(penalty, "%f", psw);
    strcat (fil,penalty);
    strcat (fil,".txt");


    FILE* ap = fopen(fil, "a");

    for ( int i = 0; i < neighborhood->nNeighborhood; ++i ) {
        char* sigla = "";
        switch(i) {
            case seqInvert-1:
                sigla = "ISJ";
                break;
            case seqShiftJob-1:
                sigla = "OJ";
                break;
            case seqSwapJob-1:
                sigla = "STJ";
                break;
            case seqShiftProj-1:
                sigla = "OP";
                break;
            case seqSwapProj-1:
                sigla = "SCTP";
                break;
            case seqCompactProj-1:
                sigla = "CPP";
                break;
            case changeOneMode-1:
                sigla = "COM";
                break;
            case changeTwoMode-1:
                sigla = "CTM";
                break;
            case changeThreeMode-1:
                sigla = "CThM";
                break;
            case changeFourMode-1:
                sigla = "CFM";
                break;
            case seqSwapJobFILS-1:
                sigla = "SSJW";
                break;
            case seqInsertJobFILS-1:
                sigla = "SSIW";
                break;
            case seqCompOnExtrem-1:
                sigla = "SPE";
                break;
            case seqMoveProj-1:
                sigla = "CSP";
                break;
        }


        if(sw) {
#ifdef DEBUG
            //      printf("\n%d %s %.8Lf", neighborhood->timePrint, sigla, neighborhood->PE[i]  );
            printf("\n%f %s %.8Lf", timePrint, sigla, neighborhood->PE[i]  );
#endif
            //    fprintf(ap,"\n%d %s %.8Lf", neighborhood->timePrint, sigla, neighborhood->PE[i] );
            fprintf(ap,"\n%f %s %.8Lf", timePrint, sigla, neighborhood->PE[i] );
        } else {
#ifdef DEBUG
            // printf("\n%d %s %.8Lf", neighborhood->timePrint, sigla, neighborhood->PI[i]  );
            printf("\n%f %s %.8Lf", timePrint, sigla, neighborhood->PI[i]  );
#endif
            //fprintf(ap,"\n%d %s %.8Lf", neighborhood->timePrint, sigla, neighborhood->PI[i] );
            fprintf(ap,"\n%f %s %.8Lf", timePrint, sigla, neighborhood->PI[i] );
        }
    }
    fclose(ap);
}

void Neighbor_printProbabilities( Neighborhood *neighborhood, int sw, int it, float psw, char *argv)
{

    char fil[256] = "probabilitiesOnline_";
    char bufferII[256];
    char side[256];
    char penalty[256];
    char instance[256];
    sprintf(instance, "%s", argv);
    strcat (fil,instance);
    strcat (fil,"_");
    sprintf(bufferII, "%d", it);
    strcat (fil,bufferII);
    strcat (fil,"_");
    sprintf(side, "%d", sw);
    strcat (fil,side);
    strcat (fil,"_");
    sprintf(penalty, "%f", psw);
    strcat (fil,penalty);
    strcat (fil,".txt");


    FILE* ap = fopen(fil, "a");

    for ( int i = 0; i < neighborhood->nNeighborhood; ++i ) {
        char* sigla = "";
        switch(i) {
            case seqInvert-1:
                sigla = "ISJ";
                break;
            case seqShiftJob-1:
                sigla = "OJ";
                break;
            case seqSwapJob-1:
                sigla = "STJ";
                break;
            case seqShiftProj-1:
                sigla = "OP";
                break;
            case seqSwapProj-1:
                sigla = "SCTP";
                break;
            case seqCompactProj-1:
                sigla = "CPP";
                break;
            case changeOneMode-1:
                sigla = "COM";
                break;
            case changeTwoMode-1:
                sigla = "CTM";
                break;
            case changeThreeMode-1:
                sigla = "CThM";
                break;
            case changeFourMode-1:
                sigla = "CFM";
                break;
            case seqSwapJobFILS-1:
                sigla = "SSJW";
                break;
            case seqInsertJobFILS-1:
                sigla = "SSIW";
                break;
            case seqCompOnExtrem-1:
                sigla = "SPE";
                break;
            case seqMoveProj-1:
                sigla = "CSP";
                break;
        }


        if(sw) {
#ifdef DEBUG
            printf("\n%d %s %.8Lf", neighborhood->timePrint, sigla, neighborhood->PE[i]  );

#endif
            fprintf(ap,"\n%d %s %.8Lf", neighborhood->timePrint, sigla, neighborhood->PE[i] );

        } else {
#ifdef DEBUG
            printf("\n%d %s %.8Lf", neighborhood->timePrint, sigla, neighborhood->PI[i]  );

#endif
            fprintf(ap,"\n%d %s %.8Lf", neighborhood->timePrint, sigla, neighborhood->PI[i] );

        }
    }
    fclose(ap);
}


void Neighbor_updatesIntensities(Neighborhood *neighborhood, int sw)
{
    assert(neighborhood != NULL);
    for(int i = 0 ; i < neighborhood->nNeighborhood ; i++) {
        if(sw) {
            neighborhood->intensity[i] = neighborhood->normFE[i];
            // printf("\nIntensitie %d %.8f", i, neighborhood->intensity[i]);
        } else {
            neighborhood->intensity[i] = neighborhood->normFI[i];
            // printf("\nIntensitie %d %.8f", i, neighborhood->intensity[i]);
        }
    }
}

void Neighbor_incrementI(Neighborhood *neighbor, int i)
{
    assert( neighbor != NULL );
    assert( i >= 0 );
    assert( i < Neighbor_nNeighborhood(neighbor) );

    neighbor->I[i] += 1;

}

void Neighbor_incrementEQ(Neighborhood *neighbor, int i)
{
    assert( neighbor != NULL );
    assert( i >= 0 );
    assert( i < Neighbor_nNeighborhood(neighbor) );


    neighbor->EQ[i] += 0.1;

}

long double Neighbor_getI(Neighborhood *neighbor, int i)
{
    assert( neighbor != NULL );
    assert( i >= 0 );
    assert( i < Neighbor_nNeighborhood(neighbor) );


    return neighbor->I[i];

}

long double Neighbor_getEQ(Neighborhood *neighbor, int i)
{
    assert( neighbor != NULL );
    assert( i >= 0 );
    assert( i < Neighbor_nNeighborhood(neighbor) );

    return neighbor->EQ[i];

}

long double Neighbor_getTI(Neighborhood *neighbor, int i)
{
    assert( neighbor != NULL );
    assert( i >= 0 );
    assert( i < Neighbor_nNeighborhood(neighbor) );

    return neighbor->TI[i];
}

long double Neighbor_getTIV(Neighborhood *neighbor, int i)
{
    assert( neighbor != NULL );
    assert( i >= 0 );
    assert( i < Neighbor_nNeighborhood(neighbor) );

    return neighbor->TIV[i];
}

long double Neighbor_getTE(Neighborhood *neighbor, int i)
{
    assert( neighbor != NULL );
    assert( i >= 0 );
    assert( i < Neighbor_nNeighborhood(neighbor) );

    return neighbor->TE[i];
}

void Neighbor_setTI(Neighborhood *neighbor, int i, long double value)
{
    assert( neighbor != NULL );
    assert( i >= 0 );
    assert( i < Neighbor_nNeighborhood(neighbor) );

    neighbor->TI[i] += value;

}

void Neighbor_setFI(Neighborhood *neighbor, int i, long double value)
{
    assert( neighbor != NULL );
    assert( i >= 0 );
    assert( i < Neighbor_nNeighborhood(neighbor) );

    neighbor->FI[i] = value;

}

void Neighbor_setFE(Neighborhood *neighbor, int i, long  double value)
{
    assert( neighbor != NULL );
    assert( i >= 0 );
    assert( i < Neighbor_nNeighborhood(neighbor) );

    neighbor->FE[i] = value;

}

void Neighbor_setF(Neighborhood *neighborhood)
{

    for(int n=0; n< neighborhood->nNeighborhood; ++n) {
        long double value = neighborhood->TI[n] + neighborhood->TE[n]+ neighborhood->TIV[n];
        long double valueTime =0;
        long double valueI = 0;
        long double valueE = 0;

        if(value != 0) {
            if(neighborhood->log)
                valueTime = log2(value);
            else
                valueTime = value;

            assert (valueTime>0);

            if(valueTime != 0) {
                valueI = (long double) neighborhood->I[n]/(long double) (valueTime);
                valueE = (long double) ( neighborhood->I[n]+ neighborhood->EQ[n]) / (long double)(valueTime);
            }

        }
        //printf("I %f, TI%f,TIV%Lf,TE%Lf = FI %Lf \n", neighborhood->I[n], neighborhood->TI[n], neighborhood->TIV[n], neighborhood->TE[n], valueI);

        neighborhood->FI[n] = valueI;
        neighborhood->FE[n] = valueE;
        //printf("FE %Lf, FI %Lf\n", neighborhood->FE[n], neighborhood->FI[n]);
    }
}

void Neighbor_setTIV(Neighborhood *neighbor, int i, long double value)
{
    assert( neighbor != NULL );
    assert( i >= 0 );
    assert( i < Neighbor_nNeighborhood(neighbor) );

    neighbor->TIV[i] += value;

}

void Neighbor_setTE(Neighborhood *neighbor, int i, long double value)
{
    assert( neighbor != NULL );
    assert( i >= 0 );
    assert( i < Neighbor_nNeighborhood(neighbor) );

    neighbor->TE[i] += value;

}

void Neighbor_normF(Neighborhood *neighborhood)
{

    neighborhood->maxFI =0;
    neighborhood->maxFE =0;

    for(int n=0; n< neighborhood->nNeighborhood; ++n) {
        neighborhood->maxFI = MAX(neighborhood->maxFI, neighborhood->FI[n]);
        neighborhood->maxFE = MAX(neighborhood->maxFE, neighborhood->FE[n]);
    }

    for(int i=0; i< neighborhood->nNeighborhood; ++i) {
        if(neighborhood->maxFI == 0) neighborhood->normFI[i] = 1;
        else neighborhood->normFI[i] = (long double)(neighborhood->FI[i])/(long double)neighborhood->maxFI;

        if(neighborhood->maxFE == 0) neighborhood->normFE[i] = 1;
        else neighborhood->normFE[i] = (long double)(neighborhood->FE[i])/(long double)neighborhood->maxFE;
    }

    long double sumWeightsE = 0;
    long double sumWeightsI = 0;
    for(int i=0; i< neighborhood->nNeighborhood; ++i) {
        sumWeightsE += neighborhood->normFE[i] +BETA;
        sumWeightsI += neighborhood->normFI[i] +BETA;
    }

    for(int i=0; i< neighborhood->nNeighborhood; ++i) {
        neighborhood->PE[i] = (neighborhood->normFE[i]+BETA) / sumWeightsE;
        neighborhood->PI[i] = (neighborhood->normFI[i]+BETA) / sumWeightsI;
    }

}

void Neighbor_clearImpTime(Neighborhood *neighborhood)
{

    CLEAR_VECTOR( neighborhood->I, long double, neighborhood->nNeighborhood);
    CLEAR_VECTOR( neighborhood->EQ, long double, neighborhood->nNeighborhood);
    CLEAR_VECTOR( neighborhood->TI, long double, neighborhood->nNeighborhood);
    CLEAR_VECTOR( neighborhood->TE, long double, neighborhood->nNeighborhood);
    CLEAR_VECTOR( neighborhood->TIV, long double, neighborhood->nNeighborhood);
    CLEAR_VECTOR( neighborhood->FI, long double, neighborhood->nNeighborhood);
    CLEAR_VECTOR( neighborhood->FE, long double, neighborhood->nNeighborhood);
    CLEAR_VECTOR( neighborhood->normFI, long double, neighborhood->nNeighborhood);
    CLEAR_VECTOR( neighborhood->normFE, long double, neighborhood->nNeighborhood);

    //neighborhood->maxFI =0;
    //neighborhood->maxFE =0;
    //neighborhood->minFI = INT_MAX;
    //neighborhood->minFE =INT_MAX;
    //neighborhood->intervalI =0;
    //neighborhood->intervalE =0;

}

int Neighbor_getIdxAssortment(Neighborhood *neighborhood, int idx)
{
    assert(neighborhood != NULL);
    assert(idx >= 0);
    assert(idx < neighborhood->nNeighborhood);

    return neighborhood->assortment[idx];
}

int *Neighbor_getAssortment(Neighborhood *neighborhood)
{
    assert(neighborhood != NULL);

    return neighborhood->assortment;
}

int Neighbor_nNeighborhood(Neighborhood *neighborhood)
{
    assert(neighborhood != NULL);

    return neighborhood->nNeighborhood;
}



double Neighbor_getIdP(Neighborhood *neighborhood, int idxNeighbor, int sw)
{
    assert(neighborhood != NULL);
    assert(idxNeighbor >= 0 && idxNeighbor < neighborhood->nNeighborhood);

    if(sw)
        return neighborhood->PE[idxNeighbor];
    else return neighborhood->PI[idxNeighbor];

}




long double Neighbor_getIdIntensity(Neighborhood *neighborhood, int idxNeighbor)
{
    assert(neighborhood != NULL);
    assert(idxNeighbor >= 0 && idxNeighbor < neighborhood->nNeighborhood);

    return neighborhood->intensity[idxNeighbor];

}

int Neighbor_getMinK(Neighborhood *neighborhood, int idx)
{
    assert (neighborhood != NULL);
    assert(idx >= 0 && idx < neighborhood->nNeighborhood);

    return neighborhood->neighborMinK[idx];
}

int Neighbor_getMaxK(Neighborhood *neighborhood, int idx)
{
    assert (neighborhood != NULL);
    assert(idx >= 0 && idx < neighborhood->nNeighborhood);

    return neighborhood->neighborMaxK[idx];
}

void Neighbor_setIt(Neighborhood *neighborhood)
{
    assert (neighborhood != NULL);
    neighborhood->it++;
}

int Neighbor_getPenaltyChangeMode(Neighborhood *neighborhood)
{
    assert (neighborhood != NULL);
    return neighborhood->penaltyChangeMode;
}

int Neighbor_roulette(Neighborhood *neighborhood)
{
    assert (neighborhood != NULL);

    double sumAux = 0.0, perc = 0.0;

    //for(int i = 0; i < Neighbor_nNeighborhood(neighborhood); ++i)
    //    sum += Neighbor_getIdIntensity(neighborhood,i) + BETA;

    //if(sum <= 1e-5) Neighbor_restartIntensity(neighborhood);

    perc = DOUBLE_RANDOM(0.0, 1.0);

    for(int i = 0; i < Neighbor_nNeighborhood(neighborhood); ++i) {
        sumAux += Neighbor_getIdP(neighborhood,i, neighborhood->sw);
        if(sumAux >= perc)
            return i;
    }

    return -1;

}

/* Neighborhoods stochastic of modes */
int Neighbor_random_changeOneMode(Neighborhood *neighborhood, Solution *sol, Solution *current)
{
    assert (neighborhood != NULL);
    assert(sol != NULL);
    assert(current != NULL);

    int nJobs = Inst_nJobs(Sol_inst(current));

    ModeSet *ms = Sol_getModeSet(current);

    int j = rand() %nJobs;
    const Job *job = Inst_job(Sol_inst(current), j);
    int m = Sol_getMode(current,j), valid=0;
    // int ctrl1=0;

    const Mode * mode;
    if((Job_nModes(job) > 1) && (Job_nInfeasModes(job) < (Job_nModes(job)-1))) {
        do {
            m = INT_RANDOM(Job_nModes(job));
            mode = Job_mode(job,m);
        } while( Mode_index(mode) == Modes_job(ms,Job_index(job)) || !Mode_isFeasible(Sol_inst(current),mode));
        //    ctrl1=1;
    }

    //  if(ctrl1)
    valid = Neighbor_changeOneMode(neighborhood, current, j, m);

    if(!valid) {
        Sol_cpy(current,sol);
        return 0;
    }

    neighborhood->lastJ[0] = j;
    neighborhood->newModes[0] = m;

    return 1;
}

int Neighbor_random_changeNModes(Neighborhood *neighborhood, Solution *sol,  Solution *current, int nModes, int **residency)
{

    //    int proj = INT_RANDOM_LU(0, (Inst_nProjects(neighborhood->inst)-1));
    int valid=0;

    //    ModeSet *ms = Sol_getModeSet(current);

    //printf("Proj %d, nProj %d \n", proj, Inst_nProjects(neighborhood->inst));
    //    valid = MSM_changeModes( neighborhood->mss[proj], ms, nModes, nModes, (const int **) residency );

    return valid;
}

int Neighbor_random_changeTwoMode(Neighborhood *neighborhood, Solution *sol,  Solution *current)
{
    assert (neighborhood != NULL);
    assert(sol != NULL);
    assert(current != NULL);

    int nJobs = Inst_nJobs(Sol_inst(current));

    ModeSet *ms = Sol_getModeSet(current);

    int j1 = rand() %nJobs;
    int j2;
    do {
        j2 =  rand() %nJobs;
    } while(j1==j2);


    const Job *job1 = Inst_job(Sol_inst(current), j1);
    const Job *job2 = Inst_job(Sol_inst(current), j2);
    int m1= Sol_getMode(current,j1), m2=Sol_getMode(current,j2), valid=0;
    // int ctrl1=0, ctrl2=0;

    const Mode* mode1;
    if( (Job_nModes(job1) > 1) &&(Job_nInfeasModes(job1) < (Job_nModes(job1)-1))) {

        do {
            m1 = INT_RANDOM(Job_nModes(job1));
            mode1 = Job_mode(job1,m1);
        } while( Mode_index(mode1) == Modes_job(ms,Job_index(job1))  || !Mode_isFeasible(Sol_inst(current),mode1));
        // ctrl1=1;
    }

    const Mode* mode2;
    if((Job_nModes(job2) > 1) && (Job_nInfeasModes(job2) < (Job_nModes(job2)-1))) {

        do {
            m2 = INT_RANDOM(Job_nModes(job2));
            mode2 = Job_mode(job2,m2);
        } while( Mode_index(mode2) == Modes_job(ms,Job_index(job2)) || !Mode_isFeasible(Sol_inst(current),mode2));
        //  ctrl2=1;
    }

    //if(ctrl1&&ctrl2)
    valid =  Neighbor_changeTwoMode(neighborhood, current, j1, j2, m1, m2);
    if(!valid) {
        Sol_cpy(current,sol);
        return 0;
    }

    neighborhood->lastJ[0] = j1;
    neighborhood->lastJ[1] = j2;
    neighborhood->newModes[0] = m1;
    neighborhood->newModes[1] = m2;

    return 1;
}

int Neighbor_random_changeThreeMode(Neighborhood *neighborhood, Solution *sol, Solution *current)
{
    assert (neighborhood != NULL);
    assert(sol != NULL);
    assert(current != NULL);

    int nJobs = Inst_nJobs(Sol_inst(current));

    ModeSet *ms = Sol_getModeSet(current);

    int j1 = rand() %nJobs;
    int j2;
    do {
        j2 =  rand() %nJobs;
    } while(j1==j2);
    int j3;
    do {
        j3 =  rand() %nJobs;
    } while(j2==j3 || j3==j1);


    const Job *job1 = Inst_job(Sol_inst(current), j1);
    const Job *job2 = Inst_job(Sol_inst(current), j2);
    const Job *job3 = Inst_job(Sol_inst(current), j3);
    int m1=Sol_getMode(current,j1), m2=Sol_getMode(current,j2), m3=Sol_getMode(current,j3), valid=0;
    // int ctrl1=0,ctrl2=0,ctrl3=0;

    const Mode* mode1;
    if((Job_nModes(job1) > 1) && (Job_nInfeasModes(job1) < (Job_nModes(job1)-1)) ) {
        do {
            m1 = INT_RANDOM(Job_nModes(job1));
            mode1 = Job_mode(job1,m1);
        } while( Mode_index(mode1) == Modes_job(ms,Job_index(job1))|| !Mode_isFeasible(Sol_inst(current),mode1));
        //  ctrl1=1;
    }

    const Mode* mode2;
    if((Job_nModes(job1) > 1) && (Job_nInfeasModes(job2) < (Job_nModes(job2)-1))) {
        do {
            m2 = INT_RANDOM(Job_nModes(job2));
            mode2 = Job_mode(job2,m2);
        } while( Mode_index(mode2) == Modes_job(ms,Job_index(job2))|| !Mode_isFeasible(Sol_inst(current),mode2));
        //    ctrl2=1;
    }

    const Mode* mode3;
    if((Job_nModes(job3) > 1) && (Job_nInfeasModes(job3) < (Job_nModes(job3)-1))) {
        do {
            m3 = INT_RANDOM(Job_nModes(job3));
            mode3 = Job_mode(job3,m3);
        } while( Mode_index(mode3) == Modes_job(ms,Job_index(job3))|| !Mode_isFeasible(Sol_inst(current),mode3));
        //  ctrl3=1;
    }

    //if(ctrl1&&ctrl2&&ctrl3)
    valid =  Neighbor_changeThreeMode(neighborhood, current,Job_index(job1), Job_index(job2), Job_index(job3), m1, m2, m3);

    if(!valid) {
        Sol_cpy(current,sol);
        return 0;
    }

    neighborhood->lastJ[0] = Job_index(job1);
    neighborhood->lastJ[1] = Job_index(job2);
    neighborhood->lastJ[2] = Job_index(job3);

    neighborhood->newModes[0] = m1;
    neighborhood->newModes[1] = m2;
    neighborhood->newModes[2] = m3;

    return 1;
}

int Neighbor_random_changeFourMode(Neighborhood *neighborhood, Solution *sol, Solution *current)
{
    assert (neighborhood != NULL);
    assert(sol != NULL);
    assert(current != NULL);

    int nJobs = Inst_nJobs(Sol_inst(current));

    ModeSet *ms = Sol_getModeSet(current);


    int j1 = rand() %nJobs;
    int j2;
    do {
        j2 =  rand() %nJobs;
    } while(j1==j2);
    int j3;
    do {
        j3 =  rand() %nJobs;
    } while(j2==j3 || j3==j1);
    int j4;
    do {
        j4 =  rand() %nJobs;
    } while(j2==j4 || j3==j4 || j1==j4);


    const Job *job1 = Inst_job(Sol_inst(current), j1);
    const Job *job2 = Inst_job(Sol_inst(current), j2);
    const Job *job3 = Inst_job(Sol_inst(current), j3);
    const Job *job4 = Inst_job(Sol_inst(current), j4);
    int m1 = Sol_getMode(current,j1), m2=Sol_getMode(current,j2), m3=Sol_getMode(current,j3), m4 = Sol_getMode(current,j4);
    int valid=0;
    // int ctrl1=0, ctrl2=0,ctrl3=0,ctrl4=0;

    const Mode* mode1;
    if((Job_nModes(job1) > 1) && (Job_nInfeasModes(job1) < (Job_nModes(job1)-1))) {
        do {
            m1 = INT_RANDOM(Job_nModes(job1));
            mode1 = Job_mode(job1,m1);
        } while( Mode_index(mode1) == Modes_job(ms,Job_index(job1))|| !Mode_isFeasible(Sol_inst(current),mode1));
        //   ctrl1=1;
    }


    const Mode* mode2;
    if((Job_nModes(job2) > 1) && (Job_nInfeasModes(job2) < (Job_nModes(job2)-1))) {
        do {
            m2 = INT_RANDOM(Job_nModes(job2));
            mode2 = Job_mode(job2,m2);
        } while( Mode_index(mode2) == Modes_job(ms,Job_index(job2))|| !Mode_isFeasible(Sol_inst(current),mode2));
        //  ctrl2=1;
    }

    const Mode* mode3;
    if((Job_nModes(job3) > 1 ) && (Job_nInfeasModes(job3) < (Job_nModes(job3)-1))) {
        do {
            m3 = INT_RANDOM(Job_nModes(job3));
            mode3 = Job_mode(job3,m3);
        } while( Mode_index(mode3) == Modes_job(ms,Job_index(job3))|| !Mode_isFeasible(Sol_inst(current),mode3));
        // ctrl3=1;
    }

    const Mode* mode4;
    if((Job_nModes(job4) > 1) && (Job_nInfeasModes(job4) < (Job_nModes(job4)-1))) {
        do {
            m4 = INT_RANDOM(Job_nModes(job4));
            mode4 = Job_mode(job4,m4);
        } while( Mode_index(mode4) == Modes_job(ms,Job_index(job4))|| !Mode_isFeasible(Sol_inst(current),mode4));
        //  ctrl4=1;
    }

    // if(ctrl1&&ctrl2&&ctrl3&&ctrl4)
    valid =  Neighbor_changeFourMode(neighborhood, current, Job_index(job1), Job_index(job2), Job_index(job3), Job_index(job4), m1, m2, m3, m4);

    if(!valid) {
        Sol_cpy(current,sol);
        return 0;
    }

    neighborhood->lastJ[0] = Job_index(job1);
    neighborhood->lastJ[1] = Job_index(job2);
    neighborhood->lastJ[2] = Job_index(job3);
    neighborhood->lastJ[3] = Job_index(job4);

    neighborhood->newModes[0] = m1;
    neighborhood->newModes[1] = m2;
    neighborhood->newModes[2] = m3;
    neighborhood->newModes[3] = m4;

    return 1;

}

/* Neighborhoods stochastic of sequence */
int Neighbor_verifyPredInv( Solution *sol, int p1, int p2 )
{
    assert (sol != NULL);
    assert( p1>=0 && p1<Inst_nJobs(Sol_inst(sol)) );
    assert( p2>=0 && p2<Inst_nJobs(Sol_inst(sol)) );
    assert( p1 <= p2);

    for(int p = p2; p >= p1; p--) {
        int idxJob = Sol_getSequence(sol,p);
        const Job* job = Inst_job(Sol_inst(sol),idxJob);
        for( int pj = 0; pj<Job_nPred(job) ; pj++ ) {
            int idxJobPred = Job_pred(job,pj);
            const Job* jobPred = Inst_job(Sol_inst(sol),idxJobPred);
            if( (Sol_getPosJob(sol, Job_index(jobPred)) < p) && (Sol_getPosJob(sol,Job_index(jobPred)) >= p1))
                return 0;
        }
    }
    return 1;
}

int Neighbor_random_inv(Neighborhood *neighborhood,Solution *current)
{
    assert (neighborhood != NULL);
    assert( current!=NULL);

    int *sequence = Sol_sequence(current);
    int nJobs = Inst_nJobs(Sol_inst(current));

    int k = INT_RANDOM_LU(Neighbor_getMinK(neighborhood, seqInvert-1),Neighbor_getMaxK(neighborhood, seqInvert-1));
    while(k >= nJobs) k--;


#ifdef DEBUG
    assert(k>2);
    assert(k<nJobs);
#endif

    int pos=0, aux;

    pos = rand() % nJobs;

    while(pos > nJobs - k)
        pos = rand() % nJobs;

    int valid = Neighbor_verifyPredInv( current, pos, pos+k-1 );

    if(valid) {
        for(int i = 0; i < k/2; ++i) {
            neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[pos + i];
            neighborhood->posNLastJModify[neighborhood->contLastJ] = pos + i;
            neighborhood->contLastJ++;
            aux = sequence[pos + i];
            sequence[pos + i] = sequence[pos + k - i - 1];
            sequence[pos + k - i - 1] = aux;
        }

        //if(neighborhood->contLastJ > Inst_nJobs(Sol_inst(current))*Inst_nJobs(Sol_inst(current)))
        //          printf( "inv Cont %d \n", neighborhood->contLastJ);
        return 1;
    }
    return 0;
}

int Neighbor_verifyPredShiftJobAhead( Solution *sol, int p1, int p2 )
{
    assert (sol != NULL);
    assert( p1>=0 && p1<Inst_nJobs(Sol_inst(sol)) );
    assert( p2>=0 && p2<Inst_nJobs(Sol_inst(sol)) );
    assert( p1 <= p2);

    int idxJob = Sol_getSequence(sol,p1);
    const Job* job = Inst_job(Sol_inst(sol),idxJob);
    for( int sj = 0; sj<Job_nSucc(job) ; sj++ ) {
        int idxJobSucc = Job_succ(job,sj);
        const Job* jobSucc = Inst_job(Sol_inst(sol),idxJobSucc);
        if( (Sol_getPosJob(sol, Job_index(jobSucc)) > p1) && (Sol_getPosJob(sol,Job_index(jobSucc)) <= p2))
            return 0;
    }

    return 1;
}

int Neighbor_verifyPredShiftJobBack( Solution *sol, int p1, int p2 )
{
    assert (sol != NULL);
    assert( p1>=0 && p1<Inst_nJobs(Sol_inst(sol)) );
    assert( p2>=0 && p2<Inst_nJobs(Sol_inst(sol)) );
    assert( p1 >= p2);

    int idxJob = Sol_getSequence(sol,p1);
    const Job* job = Inst_job(Sol_inst(sol),idxJob);
    for( int pj = 0; pj<Job_nPred(job) ; pj++ ) {
        int idxJobPred = Job_pred(job,pj);
        const Job* jobPred = Inst_job(Sol_inst(sol),idxJobPred);
        if( (Sol_getPosJob(sol, Job_index(jobPred)) < p1) && (Sol_getPosJob(sol,Job_index(jobPred)) >= p2))
            return 0;
    }

    return 1;
}

int Neighbor_random_shiftJob(Neighborhood *neighborhood, Solution *current)
{
    assert (neighborhood != NULL);
    assert( current!=NULL);

    int *sequence = Sol_sequence(current);
    int nJobs = Inst_nJobs(Sol_inst(current));

    int k = INT_RANDOM_LU(Neighbor_getMinK(neighborhood,seqShiftJob-1),Neighbor_getMaxK(neighborhood, seqShiftJob-1));

    while(k >= nJobs) k--;

#ifdef DEBUG
    assert(k > 0);
    assert(k < nJobs);
#endif

    int pos=0, aux, dir;
    pos = rand() % nJobs;

    if(pos == (nJobs - 1))
        dir = -1;
    else if(pos == 0)
        dir = 1;
    else {
        dir = rand() % 100;
        if(dir > 50) dir = 1;
        else dir = -1;
    }

    int valid =0;

    aux = sequence[pos];

    if(dir == 1) {
        while(pos + k >= nJobs) --k;
        valid = Neighbor_verifyPredShiftJobAhead( current, pos, pos+k );
        if(valid) {
            for(int i = pos; i < pos + k; ++i) {

                neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[i];
                neighborhood->posNLastJModify[neighborhood->contLastJ] = i;
                neighborhood->contLastJ++;
                sequence[i] = sequence[i + 1];
            }

            neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[pos + k];
            neighborhood->posNLastJModify[neighborhood->contLastJ] = pos + k;
            neighborhood->contLastJ++;
            sequence[pos + k] = aux;
            //            if(neighborhood->contLastJ > Inst_nJobs(Sol_inst(current))*Inst_nJobs(Sol_inst(current)))
            //                  printf( "shift job Cont %d \n", neighborhood->contLastJ);

            return 1;
        }
    } else {
        while(pos - k < 0) --k;
        valid = Neighbor_verifyPredShiftJobBack( current, pos, pos-k );
        if(valid) {
            for(int i = pos; i > pos - k; --i) {
                neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[i];
                neighborhood->posNLastJModify[neighborhood->contLastJ] = i;
                neighborhood->contLastJ++;

                sequence[i] = sequence[i - 1];
            }

            neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[pos - k];
            neighborhood->posNLastJModify[neighborhood->contLastJ] = pos - k;
            neighborhood->contLastJ++;
            sequence[pos - k] = aux;

            //  if(neighborhood->contLastJ > Inst_nJobs(Sol_inst(current))*Inst_nJobs(Sol_inst(current)))
            //        printf( "shift jobCont %d \n", neighborhood->contLastJ);
            return 1;
        }
    }
    return 0;

}

int Neighbor_verifyPredSwapJob( Solution *sol, int p1, int p2 )
{
    assert (sol != NULL);
    assert( p1>=0 && p1<Inst_nJobs(Sol_inst(sol)) );
    assert( p2>=0 && p2<Inst_nJobs(Sol_inst(sol)) );
    assert( p1 <= p2);

    int idxJob1 = Sol_getSequence(sol,p1);
    const Job* job1 = Inst_job(Sol_inst(sol),idxJob1);
    for( int sj = 0; sj<Job_nSucc(job1) ; sj++ ) {
        int idxJobSucc = Job_succ(job1,sj);
        const Job* jobSucc = Inst_job(Sol_inst(sol),idxJobSucc);
        if( (Sol_getPosJob(sol, Job_index(jobSucc)) > p1) && (Sol_getPosJob(sol,Job_index(jobSucc)) <= p2))
            return 0;
    }

    int idxJob2 = Sol_getSequence(sol,p2);
    const Job* job2 = Inst_job(Sol_inst(sol),idxJob2);
    for( int pj = 0; pj<Job_nPred(job2) ; pj++ ) {
        int idxJobPred = Job_pred(job2,pj);
        const Job* jobPred = Inst_job(Sol_inst(sol),idxJobPred);
        if( (Sol_getPosJob(sol, Job_index(jobPred)) < p2) && (Sol_getPosJob(sol,Job_index(jobPred)) >= p1))
            return 0;
    }

    return 1;
}

int Neighbor_random_swapJob(Neighborhood *neighborhood, Solution *current)
{
    assert (neighborhood != NULL);
    assert( current!=NULL);

    int *sequence = Sol_sequence(current);
    int nJobs = Inst_nJobs(Sol_inst(current));


    int k = INT_RANDOM_LU(Neighbor_getMinK(neighborhood,seqSwapJob-1),Neighbor_getMaxK(neighborhood,seqSwapJob-1));
    while(k >= nJobs) k--;

#ifdef DEBUG
    assert(k > 0);
    assert(k < nJobs);
#endif

    int pos=0, aux, dir;
    pos = rand() % nJobs;


    if(pos == (nJobs - 1))
        dir = -1;
    else if(pos == 0)
        dir = 1;
    else {
        dir = rand() % 100;
        if(dir > 50) dir = 1;
        else dir = -1;
    }


    if(dir == 1) {
        while(pos + k >= nJobs) --k;
        int valid = Neighbor_verifyPredSwapJob( current, pos, pos+k );
        if(valid) {
            neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[pos];
            neighborhood->posNLastJModify[neighborhood->contLastJ] = pos;
            neighborhood->contLastJ++;
            neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[pos + k];
            neighborhood->posNLastJModify[neighborhood->contLastJ] = pos + k;
            neighborhood->contLastJ++;

            aux = sequence[pos];
            sequence[pos] = sequence[pos + k];
            sequence[pos + k] = aux;
            //  if(neighborhood->contLastJ > Inst_nJobs(Sol_inst(current))*Inst_nJobs(Sol_inst(current)))
            //        printf( "swap jobCont %d \n", neighborhood->contLastJ);

            return 1;
        }
    } else {
        while(pos - k < 0) --k;
        int valid = Neighbor_verifyPredSwapJob( current,  pos-k, pos );
        if(valid) {
            neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[pos];
            neighborhood->posNLastJModify[neighborhood->contLastJ] = pos;
            neighborhood->contLastJ++;
            neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[pos - k];
            neighborhood->posNLastJModify[neighborhood->contLastJ] = pos - k;
            neighborhood->contLastJ++;


            aux = sequence[pos];
            sequence[pos] = sequence[pos - k];
            sequence[pos - k] = aux;
            //  if(neighborhood->contLastJ > Inst_nJobs(Sol_inst(current))*Inst_nJobs(Sol_inst(current)))
            //        printf( "swap jobCont %d \n", neighborhood->contLastJ);

            return 1;
        }
    }

    return 0;

}

int Neighbor_random_shiftProj(Neighborhood *neighborhood, Solution *current)
{
    assert (neighborhood != NULL);
    assert(current!=NULL);

    int *sequence = Sol_sequence(current);
    int nJobs = Inst_nJobs(Sol_inst(current));

    int k = INT_RANDOM_LU(Neighbor_getMinK(neighborhood,seqShiftProj-1),Neighbor_getMaxK(neighborhood, seqShiftProj-1));
    int p = INT_RANDOM_LU(0,(Inst_nProjects(Sol_inst(current))-1));
    while(k >= nJobs) k--;

#ifdef DEBUG
    assert(k > 0);
    assert(k < nJobs);
#endif

    int pid = INT_MAX_M, aux, dir;

    dir = rand() % 100;
    if(dir > 50) dir = 1;
    else dir = -1;

    if(dir == -1) {
        for(int i = 0; i < nJobs; ++i) {
            const Job* job = Inst_job(Sol_inst(current),sequence[i]);
            if(Job_project(job) == p) {
                pid = i;
                break;
            }
        }
    } else {
        for(int i = nJobs - 1; i >= 0; --i) {
            const Job* job = Inst_job(Sol_inst(current),sequence[i]);
            if(Job_project(job) == p) {
                pid = i;
                break;
            }
        }
    }

    if(dir == -1) {
        while(pid - k < 0) --k;

        for(int i = pid; i < nJobs; ++i) {
            const Job* job = Inst_job(Sol_inst(current),sequence[i]);
            if(Job_project(job) == p) {
                aux = sequence[i];
                for(int j = i; j > i - k; --j) {
                    neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[j];
                    neighborhood->posNLastJModify[neighborhood->contLastJ] = j;
                    neighborhood->contLastJ++;
                    sequence[j] = sequence[j - 1];

                }
                neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[i-k];
                neighborhood->posNLastJModify[neighborhood->contLastJ] = i-k;
                neighborhood->contLastJ++;

                sequence[i - k] = aux;
            }
        }
    } else {
        while(pid + k >= nJobs) --k;
        for(int i = pid; i >= 0; --i) {
            const Job* job = Inst_job(Sol_inst(current),sequence[i]);
            if(Job_project(job) == p) {
                aux = sequence[i];
                for(int j = i; j < i + k; ++j) {

                    neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[j];
                    neighborhood->posNLastJModify[neighborhood->contLastJ] = j;
                    neighborhood->contLastJ++;

                    sequence[j] = sequence[j + 1];
                }
                neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[i+k];
                neighborhood->posNLastJModify[neighborhood->contLastJ] = i+k;
                neighborhood->contLastJ++;
                sequence[i + k] = aux;
            }
        }
    }

    //if(neighborhood->contLastJ > Inst_nJobs(Sol_inst(current))*Inst_nJobs(Sol_inst(current)))
    //              printf( "--------------shift proj Cont %d \n", neighborhood->contLastJ);

    return 1;

}

int Neighbor_random_swapJobFILS(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *current  )
{

    assert (neighborhood != NULL);
    assert (sol != NULL);
    assert (bestSol != NULL);
    assert (current != NULL);

    int j1, l = 0, tamMaxW=0, posWstart=0, posWend=0, iniW=0, aux;
    int nJobs = Inst_nJobs(Sol_inst(current));
    int *sequence = Sol_sequence(current);


    j1 = rand() % nJobs;
    int valid = 0;

    const Job* job1 = Inst_job(Sol_inst(current), sequence[j1]);
    const Job* job2;

    if(Job_nPred(job1)  == 0)
        posWstart = j1;
    else {

        for(int i = j1; i >= 0; --i) {
            job2 = Inst_job(Sol_inst(current), sequence[i]);
            if(Job_hasPred(Sol_inst(current), job1, Job_index(job2))) {
                posWstart = i;
                break;
            }
        }
    }

    if(Job_nSucc(job1) == 0)
        posWend = j1;
    else {

        for(int i = j1; i < nJobs; ++i) {
            job2 = Inst_job( Sol_inst(current), sequence[i]);
            if(Job_hasSucc(Sol_inst(current),job1,  Job_index(job2))) {
                posWend = i;
                break;
            }
        }
    }

    tamMaxW = posWend - posWstart - 1;

    if(tamMaxW < 3)
        return 0;

    do {
        l = rand() % tamMaxW;
    } while(l < 2);

    iniW = rand() % (tamMaxW - l);

    posWstart = posWstart + iniW + 1;
    posWend = posWstart + l - 1;

    for(int i = posWstart; i <= posWend; ++i) {

        if(i<j1)
            valid = Neighbor_verifyPredSwapJob( current, i, j1 );
        else
            valid = Neighbor_verifyPredSwapJob( current, j1, i );
        if(valid) {
            //            printf("cont %d\n", neighborhood->contLastJ);
            neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[i];
            neighborhood->posNLastJModify[neighborhood->contLastJ] = i;
            neighborhood->contLastJ++;
            neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[j1];
            neighborhood->posNLastJModify[neighborhood->contLastJ] = j1;
            neighborhood->contLastJ++;

            aux = sequence[i];
            sequence[i] = sequence[j1];
            sequence[j1] = aux;


            Sol_rebuild_opt(current,sol);

            if(Sol_getCost(current) < Sol_getCost(bestSol)) {

                Sol_cpy(sol,current);
                //printf("\nImprovement Inside swapJobFILS %ld to %ld \n ", Sol_getCost(bestSol), Sol_getCost(current) );
                Sol_cpy( bestSol, current );
                //if(neighborhood->contLastJ > Inst_nJobs(Sol_inst(current))*Inst_nJobs(Sol_inst(current)))
                //  printf( "swap jobs fill Cont %d \n", neighborhood->contLastJ);
                return 1;
            } else
                Sol_cpy(current,sol);
        }
    }
    //if(neighborhood->contLastJ > Inst_nJobs(Sol_inst(current))*Inst_nJobs(Sol_inst(current)))
    //              printf( "swap jobs fill Cont %d \n", neighborhood->contLastJ);

    return 0;

}

//Neighbor_random_ShiftJobFILS
int Neighbor_random_insertJobFILS(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *current  )
{

    assert (neighborhood != NULL);
    assert (sol != NULL);
    assert (bestSol != NULL);
    assert (current != NULL);

    int j1, l = 0, tamMaxW, posWstart=0, posWend=0, iniW=0, aux;
    int nJobs = Inst_nJobs(Sol_inst(current));

    int *sequence = Sol_sequence(current);

    j1 = rand() % nJobs;
    int valid = 0;

    const Job* job1 = Inst_job(Sol_inst(current), sequence[j1]);
    const Job* job2;

    if(Job_nPred(job1) == 0)
        posWstart = j1;
    else {

        for(int i = j1; i >= 0; --i) {
            job2 = Inst_job(Sol_inst(current), sequence[i]);
            if(Job_hasPred(Sol_inst(current), job1, Job_index(job2))) {
                posWstart = i;
                break;
            }
        }
    }

    if(Job_nSucc(job1) == 0)
        posWend = j1;
    else {

        for(int i = j1; i < nJobs; ++i) {
            job2 = Inst_job( Sol_inst(current), sequence[i]);
            if(Job_hasSucc(Sol_inst(current),job1,  Job_index(job2))) {
                posWend = i;
                break;
            }
        }
    }

    tamMaxW = posWend - posWstart - 1;

    if(tamMaxW < 3)
        return 0;

    do {
        l = rand() % tamMaxW;
    } while(l < 2);

    iniW = rand() % (tamMaxW - l);

    posWstart = posWstart + iniW + 1;
    posWend = posWstart + l - 1;

    for(int i = posWstart; i <= posWend; ++i) {

        if(i<j1)
            valid =  Neighbor_verifyPredShiftJobBack( current, j1, i );
        else
            valid = Neighbor_verifyPredShiftJobAhead( current, j1, i );
        if(valid) {
            aux = sequence[j1];

            if(i < j1)
                for(int j = j1; j > i; --j) {
                    neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[j];
                    neighborhood->posNLastJModify[neighborhood->contLastJ] = j;
                    neighborhood->contLastJ++;
                    sequence[j] = sequence[j - 1];
                } else if(i > j1)
                for(int j = j1; j < i; ++j) {
                    neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[j];
                    neighborhood->posNLastJModify[neighborhood->contLastJ] = j;
                    neighborhood->contLastJ++;
                    sequence[j] = sequence[j + 1];
                } else
                continue;
            neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[i];
            neighborhood->posNLastJModify[neighborhood->contLastJ] = i;
            neighborhood->contLastJ++;

            sequence[i] = aux;

            Sol_rebuild_opt(current,sol);

            if(Sol_getCost(current) < Sol_getCost(bestSol)) {


                Sol_cpy(sol,current);
                // printf("\nImprovement Inside insertJobFILS %ld to %ld \n ", Sol_getCost(bestSol), Sol_getCost(current) );
                Sol_cpy( bestSol, current );
                //if(neighborhood->contLastJ > Inst_nJobs(Sol_inst(current))*Inst_nJobs(Sol_inst(current)))
                //  printf( "insert jobs fill Cont %d \n", neighborhood->contLastJ);

                return valid;//this moviment was valid, but don't need run rebuild after, because it has already been done here.
            } else
                Sol_cpy(current,sol);

        }

    }
    //if(neighborhood->contLastJ > Inst_nJobs(Sol_inst(current))*Inst_nJobs(Sol_inst(current)))
    //              printf( "insert jobs fill Cont %d \n", neighborhood->contLastJ);

    return valid;

}

int Neighbor_random_compactOnExtreme(Neighborhood *neighborhood, Solution *current)
{
    assert (neighborhood != NULL);
    assert( current!=NULL);

    int *sequence = Sol_sequence(current);
    int nJobs = Inst_nJobs(Sol_inst(current));


    int j1, posWstart=0, posWend=0, aux, proj;

    j1 = rand() % nJobs;

    const Job* job1 = Inst_job(Sol_inst(current),sequence[j1]);
    const Job* job2;

    proj = Job_project(job1);

    for(int i = 0; i < nJobs; ++i) {
        job2 = Inst_job(Sol_inst(current),sequence[i]);
        if(Job_project(job2) == proj) {
            posWstart = i;
            break;
        }
    }

    for(int i = nJobs-1; i >= 0; --i) {
        job2 = Inst_job(Sol_inst(current),sequence[i]);
        if(Job_project(job2) == proj) {
            posWend = i;
            break;
        }
    }

    ++posWend;
    --posWstart;

    for(int i = j1-1; i > posWstart; --i) {
        job2 = Inst_job(Sol_inst(current),sequence[i]);
        if(Job_project(job2) == proj) {
            aux = sequence[i];
            for(int j = i; j < j1; ++j) {
                job2 = Inst_job(Sol_inst(current),sequence[j+1]);
                neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[j];
                neighborhood->posNLastJModify[neighborhood->contLastJ] = j;
                neighborhood->contLastJ++;
                if(Job_project(job2)== proj) {
                    sequence[j] = aux;
                    break;
                } else
                    sequence[j] = sequence[j + 1];

            }
        }

    }

    for(int i = j1+1; i < posWend; ++i) {
        job2 = Inst_job(Sol_inst(current),sequence[i]);
        if(Job_project(job2) == proj) {
            aux = sequence[i];
            for(int j = i; j > j1; --j) {
                job2 = Inst_job(Sol_inst(current),sequence[j-1]);
                neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[j];
                neighborhood->posNLastJModify[neighborhood->contLastJ] = j;
                neighborhood->contLastJ++;
                if(Job_project(job2) == proj) {
                    sequence[j] = aux;
                    break;
                } else
                    sequence[j] = sequence[j - 1];

            }
        }

    }

    //   if(neighborhood->contLastJ > Inst_nJobs(Sol_inst(current))*  Inst_nJobs(Sol_inst(current)))
    //                 printf( "-------------- compac on extrem Cont %d\n ", neighborhood->contLastJ);

    return 1;

}

int Neighbor_random_moveProj(Neighborhood *neighborhood, Solution *current)
{
    assert (neighborhood != NULL);
    assert(current!=NULL);

    int *sequence = Sol_sequence(current);
    int nJobs = Inst_nJobs(Sol_inst(current));
    int nProjects = Inst_nProjects(Sol_inst(current));

    int nProj, part, aux;
    int *projects;
    ALLOCATE_VECTOR_INI(projects, int, nProjects);

    const Job* job1;

    nProj = (rand() % nProjects) + 1;
    part = rand() % 2;
    int p = -1;
    if(part == 0) {
        job1 = Inst_job(Sol_inst(current),sequence[0]);
        p++;
        projects[p] = Job_project(job1);
        for(int i = 1; i < nJobs; ++i) {
            if(p == nProj)
                break;

            job1 = Inst_job(Sol_inst(current),sequence[i]);
            for(int j = 0; j < p; ++j) {
                if(Job_project(job1) == projects[j])
                    break;
                else if(j == p-1)
                    projects[p] = Job_project(job1);
            }
        }

        for(int pp = p-1; pp >= 0; --pp) {
            for(int i = 0; i < nJobs; ++i) {
                job1 = Inst_job(Sol_inst(current),sequence[i]);
                if(projects[pp] !=  Job_project(job1)) {
                    for(int j = i; j < nJobs; ++j) {
                        job1 = Inst_job(Sol_inst(current),sequence[j]);
                        if(projects[pp] == Job_project(job1)) {

                            aux = sequence[j];
                            for(int h = j; h > i; --h) {
                                neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[h];
                                neighborhood->posNLastJModify[neighborhood->contLastJ] = h;
                                neighborhood->contLastJ++;

                                sequence[h] = sequence[h-1];
                            }
                            neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[i];
                            neighborhood->posNLastJModify[neighborhood->contLastJ] = i;
                            neighborhood->contLastJ++;

                            sequence[i] = aux;

                            break;
                        }
                    }
                }
            }
        }

    } else {
        p++;
        job1 = Inst_job(Sol_inst(current),sequence[nJobs-1]);
        projects[p] = Job_project(job1);

        for(int i = nJobs-2; i >= 0; --i) {
            if(p == nProj)
                break;

            job1 = Inst_job(Sol_inst(current),sequence[i]);
            for(int j = 0; j < p; ++j) {
                if(Job_project(job1) == projects[j])
                    break;
                else if(j == p-1)
                    projects[p] = Job_project(job1);
            }
        }

        for(int pp = p-1; pp >= 0; --pp) {
            for(int i = nJobs-1; i >= 0; --i) {
                job1 = Inst_job(Sol_inst(current),sequence[i]);
                if(projects[p] != Job_project(job1)) {
                    for(int j = i; j >= 0; --j) {
                        job1 = Inst_job(Sol_inst(current),sequence[j]);
                        if(projects[p] == Job_project(job1)) {
                            aux = sequence[j];
                            for(int h = j; h < i; ++h) {

                                neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[h];
                                neighborhood->posNLastJModify[neighborhood->contLastJ] = h;
                                neighborhood->contLastJ++;

                                sequence[h] = sequence[h+1];
                            }
                            neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[i];
                            neighborhood->posNLastJModify[neighborhood->contLastJ] = i;
                            neighborhood->contLastJ++;

                            sequence[i] = aux;

                            break;
                        }
                    }
                }
            }
        }


    }

    //if(neighborhood->contLastJ > Inst_nJobs(Sol_inst(current))*Inst_nJobs(Sol_inst(current)))
    //              printf( "move pro Cont %d \n", neighborhood->contLastJ);

    free(projects);
    return 1;

}

int Neighbor_random_compactProj(Neighborhood *neighborhood, Solution *current)
{
    assert (neighborhood != NULL);
    assert(current!=NULL);

    double perc = DOUBLE_RANDOM(0,1);
    int p;
    if(Inst_nProjects(Sol_inst(current))==1)
        p = 0;
    else
        p = INT_RANDOM_LU(0,Inst_nProjects(Sol_inst(current))-1);

    int *sequence = Sol_sequence(current);

    int nJobs = Inst_nJobs(Sol_inst(current));
    int nJobsPerc = nJobs*perc;

    int *A;
    int *M1;
    int *M2;
    int *D;
    ALLOCATE_VECTOR_INI(A, int, nJobs);
    ALLOCATE_VECTOR_INI(M1, int, nJobs);
    ALLOCATE_VECTOR_INI(M2, int, nJobs);
    ALLOCATE_VECTOR_INI(D, int, nJobs);

    int controlMinJobProject = INT_MAX_M;
    int controlMaxJobProject = INT_MAX_M;

    int kJobs = nJobsPerc;
    int cM1 = 0;
    for(int i = nJobs-1; kJobs > 0 && i >= 0; i--) {
        const Job* job = Inst_job(Sol_inst(current),sequence[i]);
        if(Job_project(job) == p) {
            M1[cM1] = Job_index(job);
            cM1++;
            if(controlMaxJobProject == INT_MAX_M)
                controlMaxJobProject = i;
            controlMinJobProject = i;
            kJobs--;

        }
    }


    int cA = 0;
    int cD = 0;
    int cM2 = 0;

    for(int i = 0; i < nJobs; i++) {
        const Job* job = Inst_job(Sol_inst(current),sequence[i]);

        if(i < controlMinJobProject) {
            A[cA] = Job_index(job);
            cA++;
        }

        if(i > controlMaxJobProject) {
            D[cD] = Job_index(job);
            cD++;
        }

        if(i > controlMinJobProject &&  i < controlMaxJobProject && Job_project(job) != p) {
            M2[cM2] = Job_index(job);
            cM2++;
        }

    }

    int i = 0;
    for(int c = 0; c < cA; c++) {
        neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[i];
        neighborhood->posNLastJModify[neighborhood->contLastJ] = i;
        neighborhood->contLastJ++;

        sequence[i] = A[c];
        i++;
    }

    for(int c = cM1-1; c >=0 ; c--) {
        neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[i];
        neighborhood->posNLastJModify[neighborhood->contLastJ] = i;
        neighborhood->contLastJ++;

        sequence[i] = M1[c];
        i++;
    }

    for(int c = 0; c < cM2; c++) {
        neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[i];
        neighborhood->posNLastJModify[neighborhood->contLastJ] = i;
        neighborhood->contLastJ++;

        sequence[i] = M2[c];
        i++;
    }

    for(int c = 0; c < cD; c++) {
        neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[i];
        neighborhood->posNLastJModify[neighborhood->contLastJ] = i;
        neighborhood->contLastJ++;

        sequence[i] = D[c];
        i++;
    }

    free(A);
    free(M1);
    free(M2);
    free(D);


    //if(neighborhood->contLastJ > Inst_nJobs(Sol_inst(current))*Inst_nJobs(Sol_inst(current)))
    //              printf( "compact proj Cont %d \n", neighborhood->contLastJ);

    return 1;

}

int Neighbor_random_swapTwoProj(Neighborhood *neighborhood, Solution *current)
{
    assert (neighborhood != NULL);
    assert( current!=NULL);

    int p1 = INT_RANDOM_LU(0,(Inst_nProjects(Sol_inst(current))-1));
    int p2;
    do {
        p2 = INT_RANDOM_LU(0,(Inst_nProjects(Sol_inst(current))-1));
    } while(p1==p2);


    int *sequence = Sol_sequence(current);
    int nJobs = Inst_nJobs(Sol_inst(current));

    int controlFirstProject = INT_MAX_M;
    int controlLastProject = INT_MAX_M;
    int controlMinJobProject1 = INT_MAX_M;
    int controlMinJobProject2 = INT_MAX_M;
    int controlIndex = INT_MAX_M;

    int *R;
    int *M1;
    int *M2;
    int *D;
    ALLOCATE_VECTOR_INI(R, int, nJobs);
    ALLOCATE_VECTOR_INI(M1, int, nJobs);
    ALLOCATE_VECTOR_INI(M2, int, nJobs);
    ALLOCATE_VECTOR_INI(D, int, nJobs);

    int cM1 = 0, cM2 = 0, cR = 0, cD = 0;
    for(int i = nJobs-1; i >= 0; i--) {

        const Job* job =  Inst_job(Sol_inst(current),sequence[i]);
        if(Job_project(job) == p1) {
            if(controlLastProject == INT_MAX_M) {
                controlLastProject = p1;
                controlFirstProject = p2;
            }
            M1[cM1] = Job_index(job);
            cM1++;
            controlMinJobProject1 = i;
        }

        if(Job_project(job) == p2) {
            if(controlLastProject == INT_MAX_M) {
                controlLastProject = p1;
                controlFirstProject = p2;
            }
            M2[cM2] = Job_index(job);
            cM2++;
            controlMinJobProject2 = i;
        }
    }

    if(controlFirstProject == p1)  controlIndex = controlMinJobProject1;
    if(controlFirstProject == p2)  controlIndex = controlMinJobProject2;

    for(int i = 0; i < nJobs; i++) {
        const Job* job =  Inst_job(Sol_inst(current),sequence[i]);

        if(Job_project(job) != p2 && Job_project(job)!= p1 && i < controlIndex ) {
            R[cR] = Job_index(job);
            cR++;
        }

        if(Job_project(job) != p2 && Job_project(job) != p1 && i > controlIndex ) {
            D[cD] = Job_index(job);
            cD++;
        }

    }

    int i = 0;
    for(int c = 0; c < cR; c++) {
        neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[i];
        neighborhood->posNLastJModify[neighborhood->contLastJ] = i;
        neighborhood->contLastJ++;

        sequence[i] = R[c];
        i++;
    }

    if(controlLastProject == p2) {
        for(int c = cM2-1; c >= 0 ; c--) {
            neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[i];
            neighborhood->posNLastJModify[neighborhood->contLastJ] = i;
            neighborhood->contLastJ++;

            sequence[i] = M2[c];
            i++;
        }
        for(int c = cM1-1; c >= 0; c--) {
            neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[i];
            neighborhood->posNLastJModify[neighborhood->contLastJ] = i;
            neighborhood->contLastJ++;

            sequence[i] = M1[c];
            i++;
        }
    }

    if(controlLastProject == p1) {
        for(int c = cM1-1; c >=0; c--) {
            neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[i];
            neighborhood->posNLastJModify[neighborhood->contLastJ] = i;
            neighborhood->contLastJ++;

            sequence[i] = M1[c];
            i++;
        }
        for(int c = cM2-1; c >=0; c--) {
            neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[i];
            neighborhood->posNLastJModify[neighborhood->contLastJ] = i;
            neighborhood->contLastJ++;

            sequence[i] = M2[c];
            i++;
        }
    }


    for(int c = 0; c < cD; c++) {
        neighborhood->nLastJModify[neighborhood->contLastJ] = sequence[i];
        neighborhood->posNLastJModify[neighborhood->contLastJ] = i;
        neighborhood->contLastJ++;

        sequence[i] = D[c];
        i++;
    }

    free(R);
    free(M1);
    free(M2);
    free(D);

    //if(neighborhood->contLastJ > Inst_nJobs(Sol_inst(current))*Inst_nJobs(Sol_inst(current)))
    //              printf( "swap t proj Cont %d \n", neighborhood->contLastJ);


    return 1;

}

//Neighborhood parameter to call deterministic
int Neighbor_inv( Solution *sol, int i, int k)
{
    assert( sol!=NULL);

    int *sequence = Sol_sequence(sol);
    int nJobs = Inst_nJobs(Sol_inst(sol));
    while(k >= nJobs) k--;

#ifdef DEBUG
    assert(k > 2);
    assert(k < nJobs);
    assert (i >= 0 && i < nJobs);
#endif

    int aux, j;

    int valid = Neighbor_verifyPredInv( sol, i, i+k-1 );
    if(!valid) return 0;

    for (int l = i; (l <i+k) && (l-i < k/2); l++) {
        j=(i+k)-(l-i) -1;
        aux = sequence[l];
        sequence[l] = sequence[j];
        sequence[j] = aux;
    }

    return 1;

}

int Neighbor_swapJob(Neighborhood *neighborhood, Solution *sol, int j1, int j2)
{
    assert( sol!=NULL);
    assert( neighborhood!=NULL);

    assert (j1 >= 0 && j1 < Inst_nJobs(Sol_inst(sol)));
    assert (j2 >= 0 && j2 < Inst_nJobs(Sol_inst(sol)));

    int *sequence = Sol_sequence(sol);
    int aux;
    int valid1, valid2;

    if(j1<j2) {
        valid1 = Neighbor_verifyPredSwapJob( sol, j1, j2 );
        if(!valid1) return 0;
    } else {
        valid2 = Neighbor_verifyPredSwapJob( sol, j2, j1 );
        if(!valid2) return 0;
    }

    aux = sequence[j1];
    sequence[j1] = sequence[j2];
    sequence[j2] = aux;

    return 1;

}

int Neighbor_shiftJob(Neighborhood *neighborhood, Solution *sol,  int i, int k, int dir)
{
    assert( neighborhood!=NULL);
    assert( sol!=NULL);

    int *sequence = Sol_sequence(sol);
    int nJobs = Inst_nJobs(Sol_inst(sol));
    while(k >= nJobs) k--;

#ifdef DEBUG
    assert(k>0);
    assert(k < nJobs);
    assert (i >= 0 && i < nJobs);
#endif

    int aux, valid =0;

    aux = sequence[i];
    if(dir == 1) {
        while(i + k >= nJobs) --k;
        valid = Neighbor_verifyPredShiftJobAhead( sol, i, i+k );
        if(!valid) return 0;
        for(int p = i; p < i + k; ++p)
            sequence[p] = sequence[p + 1];

        sequence[i + k] = aux;
    } else {
        while(i - k < 0) --k;
        valid = Neighbor_verifyPredShiftJobBack( sol, i, i-k );
        if(!valid) return 0;
        for(int p = i; p > i - k; --p)
            sequence[p] = sequence[p - 1];

        sequence[i - k] = aux;
    }

    return 1;
}

int Neighbor_shiftProj(Neighborhood *neighborhood, Solution *sol, int k, int p, int dir)
{

    assert( neighborhood!=NULL);
    assert( sol!=NULL);

    int *sequence = Sol_sequence(sol);
    int nJobs = Inst_nJobs(Sol_inst(sol));
    while(k >= nJobs) k--;

#ifdef DEBUG
    assert(k > 0);
    assert(k < nJobs);
    assert (p >= 0 && p < Inst_nProjects(Sol_inst(sol)));
#endif

    int pid = INT_MAX_M, aux;

    if(dir == -1) {
        for(int i = 0; i < nJobs; ++i) {
            const Job* job = Inst_job(Sol_inst(sol),sequence[i]);
            if(Job_project(job) == p) {
                pid = i;
                break;
            }
        }
    } else {
        for(int i = nJobs - 1; i >= 0; --i) {
            const Job* job = Inst_job(Sol_inst(sol),sequence[i]);
            if(Job_project(job) == p) {
                pid = i;
                break;
            }
        }
    }

    if(dir == -1) {
        while(pid - k < 0) --k;

        for(int i = pid; i < nJobs; ++i) {
            const Job* job = Inst_job(Sol_inst(sol),sequence[i]);
            if(Job_project(job) == p) {
                aux = sequence[i];

                for(int j = i; j > i - k; --j)
                    sequence[j] = sequence[j - 1];
                sequence[i - k] = aux;
            }
        }
    } else {
        while(pid + k >= nJobs) --k;
        for(int i = pid; i >= 0; --i) {
            const Job* job = Inst_job(Sol_inst(sol),sequence[i]);
            if(Job_project(job) == p) {
                aux = sequence[i];
                for(int j = i; j < i + k; ++j)
                    sequence[j] = sequence[j + 1];
                sequence[i + k] = aux;
            }
        }
    }
    return 1;
}

int Neighbor_swapTwoProj(Neighborhood *neighborhood, Solution *sol, double p1, int p2)
{
    assert( neighborhood!=NULL);
    assert( sol!=NULL);
    assert (p1 >= 0 && p1 < Inst_nProjects(Sol_inst(sol)));
    assert (p2 >= 0 && p2 < Inst_nProjects(Sol_inst(sol)));

    int *sequence = Sol_sequence(sol);
    int nJobs = Inst_nJobs(Sol_inst(sol));

    int controlFirstProject = INT_MAX_M;
    int controlLastProject = INT_MAX_M;
    int controlMinJobProject1 = INT_MAX_M;
    int controlMinJobProject2 = INT_MAX_M;
    int controlIndex = INT_MAX_M;

    int *L;
    int *M1;
    int *M2;
    int *R;
    ALLOCATE_VECTOR_INI(L, int, nJobs);
    ALLOCATE_VECTOR_INI(M1, int, nJobs);
    ALLOCATE_VECTOR_INI(M2, int, nJobs);
    ALLOCATE_VECTOR_INI(R, int, nJobs);

    int cM1 = 0, cM2 = 0, cL = 0, cR = 0;
    for(int i = nJobs-1; i >= 0; i--) {

        const Job* job =  Inst_job(Sol_inst(sol),sequence[i]);
        if(Job_project(job) == p1) {
            if(controlLastProject == INT_MAX_M) {
                controlLastProject = p1;
                controlFirstProject = p2;
            }
            M1[cM1] = Job_index(job);
            cM1++;
            controlMinJobProject1 = i;
        }

        if(Job_project(job) == p2) {
            if(controlLastProject == INT_MAX_M) {
                controlLastProject = p1;
                controlFirstProject = p2;
            }
            M2[cM2] = Job_index(job);
            cM2++;
            controlMinJobProject2 = i;
        }
    }

    if(controlFirstProject == p1)  controlIndex = controlMinJobProject1;
    if(controlFirstProject == p2)  controlIndex = controlMinJobProject2;

    for(int i = 0; i < nJobs; i++) {
        const Job* job =  Inst_job(Sol_inst(sol),sequence[i]);

        if(Job_project(job) != p2 && Job_project(job)!= p1 && i < controlIndex ) {
            L[cL] = Job_index(job);
            cL++;
        }


        if(Job_project(job) != p2 && Job_project(job) != p1 && i > controlIndex ) {
            R[cR] = Job_index(job);
            cR++;
        }

    }

    int i = 0;
    for(int c = 0; c < cL; c++) {
        sequence[i] = L[c];
        i++;
    }

    if(controlLastProject == p2) {
        for(int c = cM2-1; c >= 0; c--) {
            sequence[i] = M2[c];
            i++;
        }
        for(int c = cM1-1; c >= 0; c--) {
            sequence[i] = M1[c];
            i++;
        }
    }

    if(controlLastProject == p1) {
        for(int c = cM1-1; c >= 0; c--) {
            sequence[i] = M1[c];
            i++;
        }
        for(int c = cM2-1; c >= 0; c--) {
            sequence[i] = M2[c];
            i++;
        }
    }


    for(int c = 0; c < cR; c++) {
        sequence[i] = R[c];
        i++;
    }

    free(L);
    free(M1);
    free(M2);
    free(R);
    return 1;
}

int Neighbor_compactProj(Neighborhood *neighborhood, Solution *sol, double perc, int p)
{

    assert( neighborhood!=NULL);
    assert( sol!=NULL);
    assert (p >= 0 && p < Inst_nProjects(Sol_inst(sol)));

    int *sequence = Sol_sequence(sol);

    int nJobs = Inst_nJobs(Sol_inst(sol));
    int nJobsPerc = nJobs*perc;

    //  do {

    int *A;
    int *M1;
    int *M2;
    int *D;
    ALLOCATE_VECTOR_INI(A, int, nJobs);
    ALLOCATE_VECTOR_INI(M1, int, nJobs);
    ALLOCATE_VECTOR_INI(M2, int, nJobs);
    ALLOCATE_VECTOR_INI(D, int, nJobs);



    int controlMinJobProject = INT_MAX_M;
    int controlMaxJobProject = INT_MAX_M;

    int kJobs = nJobsPerc; //nJobs;
    int cM1 = 0;
    for(int i = nJobs-1; kJobs > 0 && i >= 0; i--) {
        const Job* job = Inst_job(Sol_inst(sol),sequence[i]);
        if(Job_project(job) == p) {
            M1[cM1] = Job_index(job);
            cM1++;
            if(controlMaxJobProject == INT_MAX_M)
                controlMaxJobProject = i;
            controlMinJobProject = i;
            kJobs--;

        }
    }


    int cA = 0;
    int cD = 0;
    int cM2 = 0;

    for(int i = 0; i < nJobs; i++) {
        const Job* job = Inst_job(Sol_inst(sol),sequence[i]);

        if(i < controlMinJobProject) {
            A[cA] = Job_index(job);
            cA++;
        }

        if(i > controlMaxJobProject) {
            D[cD] = Job_index(job);
            cD++;
        }

        if(i > controlMinJobProject &&  i < controlMaxJobProject && Job_project(job) != p) {
            M2[cM2] = Job_index(job);
            cM2++;
        }

    }

    int i = 0;
    for(int c = 0; c < cA; c++) {
        sequence[i] = A[c];
        i++;
    }

    for(int c = cM1-1; c >= 0; c--) {
        sequence[i] = M1[c];
        i++;
    }


    for(int c = 0; c < cM2; c++) {
        sequence[i] = M2[c];
        i++;
    }

    for(int c = 0; c < cD; c++) {
        sequence[i] = D[c];
        i++;
    }

    free(A);
    free(M1);
    free(M2);
    free(D);
    return 1;
    //    nJobsPerc--;

    //  } while(nJobsPerc>1);
}

int Neighbor_swapJobFILS(Neighborhood *neighborhood, Solution *current, Solution *sol, Solution *bestSol, int j1 )
{

    assert( neighborhood!=NULL);
    assert( sol!=NULL);
    assert( current!=NULL);
    assert( bestSol!=NULL);

    int l = 0, tamMaxW=0, posWstart=0, posWend=0, iniW=0, aux;
    int nJobs = Inst_nJobs(Sol_inst(current));
    int *sequence = Sol_sequence(current);

    assert (j1 >= 0 &&  j1 < nJobs);

    const Job* job1 = Inst_job(Sol_inst(current), sequence[j1]);
    const Job* job2;

    if(Job_nPred(job1)  == 0)
        posWstart = j1;
    else {

        for(int i = j1; i >= 0; --i) {
            job2 = Inst_job(Sol_inst(current), sequence[i]);
            if(Job_hasPred(Sol_inst(current), job1, Job_index(job2))) {
                posWstart = i;
                break;
            }
        }
    }

    if(Job_nSucc(job1) == 0)
        posWend = j1;
    else {

        for(int i = j1; i < nJobs; ++i) {
            job2 = Inst_job( Sol_inst(current), sequence[i]);
            if(Job_hasSucc(Sol_inst(current),job1,  Job_index(job2))) {
                posWend = i;
                break;
            }
        }
    }

    tamMaxW = posWend - posWstart - 1;

    if(tamMaxW < 3)
        return 0;

    l = tamMaxW/2;

    //l = rand() % tamMaxW +2;

    iniW = tamMaxW-l;

    //iniW = rand() % (tamMaxW - l);

    posWstart = posWstart + iniW + 1;
    posWend = posWstart + l - 1;

    int valid = 0;

    for(int i = posWstart; i <= posWend; ++i) {

        if(i<j1)
            valid = Neighbor_verifyPredSwapJob( current, i, j1 );
        else
            valid = Neighbor_verifyPredSwapJob( current, j1, i );

        if(valid) {

            aux = sequence[i];
            sequence[i] = sequence[j1];
            sequence[j1] = aux;

            Sol_rebuild_opt(current, sol);
            Neighbor_setIt(neighborhood);

            if(Sol_getCost(current) < Sol_getCost(bestSol)) {
                Sol_cpy(sol, current);
                // printf("\nImprovement swapJobFILS %ld to %ld \n ", Sol_getCost(bestSol), Sol_getCost(current) );
                Sol_cpy( bestSol, current );
                return valid;
            } else
                Sol_cpy(current,sol);
        }
    }
    return valid;

}

int Neighbor_insertJobFILS(Neighborhood *neighborhood, Solution *current, Solution *sol, Solution *bestSol, int j1 )
{

    assert( neighborhood!=NULL);
    assert( sol!=NULL);
    assert( current!=NULL);
    assert( bestSol!=NULL);

    int l = 0, tamMaxW=0, posWstart=0, posWend=0, iniW=0, aux;
    int nJobs = Inst_nJobs(Sol_inst(current));

    assert (j1 >= 0 &&  j1 < nJobs);

    int *sequence = Sol_sequence(current);

    const Job* job1 = Inst_job(Sol_inst(current), sequence[j1]);
    const Job* job2;

    if(Job_nPred(job1) == 0)
        posWstart = j1;
    else {

        for(int i = j1; i >= 0; --i) {
            job2 = Inst_job(Sol_inst(current), sequence[i]);
            if(Job_hasPred(Sol_inst(current), job1, Job_index(job2))) {
                posWstart = i;
                break;
            }
        }
    }

    if(Job_nSucc(job1) == 0)
        posWend = j1;
    else {

        for(int i = j1; i < nJobs; ++i) {
            job2 = Inst_job( Sol_inst(current), sequence[i]);
            if(Job_hasSucc(Sol_inst(current),job1,  Job_index(job2))) {
                posWend = i;
                break;
            }
        }
    }

    tamMaxW = posWend - posWstart - 1;

    if(tamMaxW < 3)
        return 0;

    //do {
    //   l = rand() % tamMaxW;
    //} while(l < 2);

    //l = rand() % tamMaxW +2;

    l = tamMaxW/2;


    iniW = tamMaxW-l;

    //iniW = rand() % (tamMaxW - l);

    posWstart = posWstart + iniW + 1;
    posWend = posWstart + l - 1;

    int valid = 0;

    for(int i = posWstart; i <= posWend; ++i) {

        if(i<j1)
            valid =  Neighbor_verifyPredShiftJobBack( current, j1, i );
        else
            valid = Neighbor_verifyPredShiftJobAhead( current, j1, i );

        if(valid) {

            aux = sequence[j1];

            if(i < j1)
                for(int j = j1; j > i; --j)
                    sequence[j] = sequence[j - 1];
            else if(i > j1)
                for(int j = j1; j < i; ++j)
                    sequence[j] = sequence[j + 1];
            else
                continue;

            sequence[i] = aux;

            Sol_rebuild_opt(current, sol);
            Neighbor_setIt(neighborhood);

            if(Sol_getCost(current) < Sol_getCost(bestSol)) {
                Sol_cpy(sol, current);
                // printf("\nImprovement insertJobFILS %ld to %ld \n ", Sol_getCost(bestSol), Sol_getCost(current) );
                Sol_cpy( bestSol, current );
                return valid;
            } else
                Sol_cpy(current,sol);
        }
    }

    return valid;

}

int Neighbor_moveProj(Neighborhood *neighborhood, Solution *sol, int nP, int dir)
{

    assert( neighborhood!=NULL);
    assert( sol!=NULL);

    int *sequence = Sol_sequence(sol);
    int nJobs = Inst_nJobs(Sol_inst(sol));
    int nProjects = Inst_nProjects(Sol_inst(sol));

    assert (nP >= 0 &&  nP <= nProjects);

    int aux;
    int *projects;
    ALLOCATE_VECTOR_INI(projects, int, nProjects);

    const Job* job1;

    int p = -1;
    if(dir == 0) {
        job1 = Inst_job(Sol_inst(sol),sequence[0]);
        p++;
        projects[p] = Job_project(job1);
        for(int i = 1; i < nJobs; ++i) {
            if(p == nP)
                break;

            job1 = Inst_job(Sol_inst(sol),sequence[i]);
            for(int j = 0; j < p; ++j) {
                if(Job_project(job1) == projects[j])
                    break;
                else if(j == p-1)
                    projects[p] = Job_project(job1);
            }
        }

        for(int pp = p-1; pp >= 0; --pp) {
            for(int i = 0; i < nJobs; ++i) {
                job1 = Inst_job(Sol_inst(sol),sequence[i]);
                if(projects[pp] !=  Job_project(job1)) {
                    for(int j = i; j < nJobs; ++j) {
                        job1 = Inst_job(Sol_inst(sol),sequence[j]);
                        if(projects[pp] == Job_project(job1)) {
                            aux = sequence[j];
                            for(int h = j; h > i; --h)
                                sequence[h] = sequence[h-1];

                            sequence[i] = aux;

                            break;
                        }
                    }
                }
            }
        }

    } else {
        p++;
        job1 = Inst_job(Sol_inst(sol),sequence[nJobs-1]);
        projects[p] = Job_project(job1);

        for(int i = nJobs-2; i >= 0; --i) {
            if(p == nP)
                break;

            job1 = Inst_job(Sol_inst(sol),sequence[i]);
            for(int j = 0; j < p; ++j) {
                if(Job_project(job1) == projects[j])
                    break;
                else if(j == p-1)
                    projects[p] = Job_project(job1);
            }
        }

        for(int pp = p-1; pp >= 0; --pp) {
            for(int i = nJobs-1; i >= 0; --i) {
                job1 = Inst_job(Sol_inst(sol),sequence[i]);
                if(projects[p] != Job_project(job1)) {
                    for(int j = i; j >= 0; --j) {
                        job1 = Inst_job(Sol_inst(sol),sequence[j]);
                        if(projects[p] == Job_project(job1)) {
                            aux = sequence[j];
                            for(int h = j; h < i; ++h)
                                sequence[h] = sequence[h+1];

                            sequence[i] = aux;

                            break;
                        }
                    }
                }
            }
        }


    }
    free(projects);
    return 1;
}

int Neighbor_compactOnExtreme(Neighborhood *neighborhood, Solution *sol, int j1)
{
    assert( neighborhood!=NULL);
    assert( sol!=NULL);

    int *sequence = Sol_sequence(sol);
    int nJobs = Inst_nJobs(Sol_inst(sol));

    assert (j1 >= 0 &&  j1 < nJobs);

    int posWstart=0, posWend=0, aux, proj;

    const Job* job1 = Inst_job(Sol_inst(sol),sequence[j1]);
    const Job* job2;

    proj = Job_project(job1);

    for(int i = 0; i < nJobs; ++i) {
        job2 = Inst_job(Sol_inst(sol),sequence[i]);
        if(Job_project(job2) == proj) {
            posWstart = i;
            break;
        }
    }

    for(int i = nJobs-1; i >= 0; --i) {
        job2 = Inst_job(Sol_inst(sol),sequence[i]);
        if(Job_project(job2) == proj) {
            posWend = i;
            break;
        }
    }

    ++posWend;
    --posWstart;

    for(int i = j1-1; i > posWstart; --i) {
        job2 = Inst_job(Sol_inst(sol),sequence[i]);
        if(Job_project(job2) == proj) {
            aux = sequence[i];
            for(int j = i; j < j1; ++j) {
                job2 = Inst_job(Sol_inst(sol),sequence[j+1]);
                if(Job_project(job2)== proj) {
                    sequence[j] = aux;
                    break;
                } else
                    sequence[j] = sequence[j + 1];

            }
        }

    }

    for(int i = j1+1; i < posWend; ++i) {
        job2 = Inst_job(Sol_inst(sol),sequence[i]);
        if(Job_project(job2) == proj) {
            aux = sequence[i];
            for(int j = i; j > j1; --j) {
                job2 = Inst_job(Sol_inst(sol),sequence[j-1]);
                if(Job_project(job2) == proj) {
                    sequence[j] = aux;
                    break;
                } else
                    sequence[j] = sequence[j - 1];

            }
        }

    }

    return 1;

}

int Neighbor_changeOneMode(Neighborhood *neighborhood, Solution *sol, int j, int m)
{
    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert (j >= 0 &&  j < Inst_nJobs(Sol_inst(sol)));
    assert (m >= 0 &&  m < 3);

    ModeSet *ms = Sol_getModeSet(sol);

    int control = 0;

    if(Modes_modifyAndVerify( ms, j, m))
        control=1;

    return control;

}

int Neighbor_changeTwoMode(Neighborhood *neighborhood, Solution *sol, int j1, int j2, int m1, int m2)
{
    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert (j1 >= 0 &&  j1 < Inst_nJobs(Sol_inst(sol)));
    assert (j2 >= 0 &&  j2 < Inst_nJobs(Sol_inst(sol)));
    assert (m1 >= 0 &&  m1 < 3);
    assert (m2 >= 0 &&  m2 < 3);

    ModeSet *ms = Sol_getModeSet(sol);

    int ctr1 = 0, ctr2 =0;

    if(Modes_modifyAndVerify( ms, j1, m1))
        ctr1=1;

    if(Modes_modifyAndVerify( ms, j2, m2))
        ctr2=1;

    if(ctr1 && ctr2) return 1;

    return 0;

}

int Neighbor_changeThreeMode(Neighborhood *neighborhood, Solution *sol, int j1, int j2, int j3, int m1, int m2, int m3)
{
    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert (j1 >= 0 &&  j1 < Inst_nJobs(Sol_inst(sol)));
    assert (j2 >= 0 &&  j2 < Inst_nJobs(Sol_inst(sol)));
    assert (j3 >= 0 &&  j3 < Inst_nJobs(Sol_inst(sol)));
    assert (m1 >= 0 &&  m1 < 3);
    assert (m2 >= 0 &&  m2 < 3);
    assert (m3 >= 0 &&  m3 < 3);


    ModeSet *ms = Sol_getModeSet(sol);


    int ctr1 = 0, ctr2 =0, ctr3 = 0;

    if(Modes_modifyAndVerify( ms,j1, m1))
        ctr1=1;

    if(Modes_modifyAndVerify( ms, j2, m2))
        ctr2=1;

    if(Modes_modifyAndVerify( ms, j3, m3))
        ctr3=1;

    if(ctr1 && ctr2 && ctr3) return 1;

    return 0;

}

int Neighbor_changeFourMode(Neighborhood *neighborhood, Solution *sol,int j1, int j2, int j3, int j4, int m1, int m2, int m3, int m4)
{
    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert (j1 >= 0 &&  j1 < Inst_nJobs(Sol_inst(sol)));
    assert (j2 >= 0 &&  j2 < Inst_nJobs(Sol_inst(sol)));
    assert (j3 >= 0 &&  j3 < Inst_nJobs(Sol_inst(sol)));
    assert (j4 >= 0 &&  j4 < Inst_nJobs(Sol_inst(sol)));
    assert (m1 >= 0 &&  m1 < 3);
    assert (m2 >= 0 &&  m2 < 3);
    assert (m3 >= 0 &&  m3 < 3);
    assert (m4 >= 0 &&  m4 < 3);

    ModeSet *ms = Sol_getModeSet(sol);

    int ctr1 = 0, ctr2 =0, ctr3 = 0, ctr4 = 0;


    if(Modes_modifyAndVerify( ms, j1, m1))
        ctr1=1;

    if(Modes_modifyAndVerify( ms, j2, m2))
        ctr2=1;

    if(Modes_modifyAndVerify( ms, j3, m3))
        ctr3=1;

    if(Modes_modifyAndVerify( ms,j4, m4))
        ctr4=1;


    if(ctr1 && ctr2 && ctr3 && ctr4) return 1;

    return 0;
}

//Neighborhood search deterministic
int Neighbor_search_inv(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *current, int kN, int firstImp, double timeRem)
{

    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(current != NULL);
    assert(kN >= 0 && kN < neighborhood->nNeighborhood);

    double _time = 0;
    clock_t tStart = clock();
    int nJobs = Inst_nJobs( Sol_inst(sol));
    int improve = 0, valid =0;
    for(int k = Neighbor_getMinK(neighborhood, seqInvert-1); k < Neighbor_getMaxK(neighborhood, seqInvert-1); k++) {
        for(int i = 0; i < nJobs-k; i++) {

            valid = Neighbor_inv(current, i, k);
            if(!valid) continue;
            Sol_rebuild_opt(current, sol);

            if(Sol_getCost(current) < Sol_getCost(sol) ) {
                Sol_cpy(sol,current);
                if(Sol_getCost(current) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                    _time = ( (double)clock()/CLOCKS_PER_SEC );
                    printf("\n%ld %f", Sol_getCost(bestSol), _time);
                    fflush(stdout);
                    fflush(stdin);

#endif // DEBUG
                    Sol_cpy( bestSol, current );
                }
                improve = 1;
                if(firstImp && kN)  return improve;
            } else
                Sol_cpy(current,sol);

            _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            if( _time > timeRem) return improve;

        }
    }

    return improve;

}

int Neighbor_search_shiftJob(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *current, int kN, int firstImp, double timeRem)
{
    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(current != NULL);
    assert(kN >= 0 && kN < neighborhood->nNeighborhood);

    double _time = 0;
    clock_t tStart = clock();
    int nJobs = Inst_nJobs( Sol_inst(sol));
    int improve = 0, valid = 0;
    for(int k = Neighbor_getMinK(neighborhood, seqShiftJob-1); k < Neighbor_getMaxK(neighborhood, seqShiftJob-1); k++) {
        for(int dir = -1; dir < 1; dir++ ) {
            for(int i = 0; i < nJobs-k; i++) {

                valid = Neighbor_shiftJob( neighborhood, current, i, k, dir);
                if(!valid) continue;
                Sol_rebuild_opt(current, sol);

                if(Sol_getCost(current) < Sol_getCost(sol) ) {
                    Sol_cpy(sol,current);
                    if(Sol_getCost(current) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                        _time = ( (double)clock()/CLOCKS_PER_SEC );
                        printf("\n%ld %f", Sol_getCost(bestSol), _time);
                        fflush(stdout);
                        fflush(stdin);

#endif // DEBUG
                        Sol_cpy( bestSol, sol );
                    }
                    improve = 1;

                    if(firstImp && kN) return improve;
                } else
                    Sol_cpy(current,sol);

                _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
                if( _time > timeRem) return improve;

            }
        }
    }
    return improve;

}

int Neighbor_search_swapJob(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *current, int kN, int firstImp, double timeRem)
{
    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(current != NULL);
    assert(kN >= 0 && kN < neighborhood->nNeighborhood);

    double _time = 0;
    clock_t tStart = clock();
    int nJobs = Inst_nJobs( Sol_inst(sol));
    int improve = 0, valid = 0;
    for(int j1 = 0; j1 < nJobs; j1++) {
        for(int j2 = j1+1; j2 < nJobs; j2++) {

            valid = Neighbor_swapJob(neighborhood, current, j1, j2);
            if(!valid) continue;
            Sol_rebuild_opt(current, sol);

            if(Sol_getCost(current) < Sol_getCost(sol) ) {
                Sol_cpy(sol,current);
                if(Sol_getCost(current) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                    _time = ( (double)clock()/CLOCKS_PER_SEC );
                    printf("\n%ld %f", Sol_getCost(bestSol), _time);
                    fflush(stdout);
                    fflush(stdin);
#endif
                    Sol_cpy( bestSol, sol );
                }
                improve = 1;

                if(firstImp && kN) return improve;
            } else
                Sol_cpy(current,sol);

            _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            if( _time > timeRem) return improve;
        }
    }

    return improve;
}

int Neighbor_search_shiftProj(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *current, int kN, int firstImp, double timeRem)
{
    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(current != NULL);
    assert(kN >= 0 && kN < neighborhood->nNeighborhood);

    double _time = 0;
    clock_t tStart = clock();
    int improve = 0, valid = 0;
    for(int k = Neighbor_getMinK(neighborhood, seqShiftProj-1); k<Neighbor_getMaxK(neighborhood, seqShiftProj-1); k++) {
        for(int dir = -1; dir < 1; dir++ ) {
            for(int p = 0; p < Inst_nProjects(Sol_inst(sol)); p++) {

                valid = Neighbor_shiftProj(neighborhood, current,k,p,dir);
                if(!valid) continue;
                Sol_rebuild_opt(current, sol);

                if(Sol_getCost(current) < Sol_getCost(sol) ) {
                    Sol_cpy(sol,current);
                    if(Sol_getCost(current) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                        _time = ( (double)clock()/CLOCKS_PER_SEC );
                        printf("\n%ld %f", Sol_getCost(bestSol), _time);
                        fflush(stdout);
                        fflush(stdin);
#endif // DEBUG
                        Sol_cpy( bestSol, sol );
                    }
                    improve = 1;

                    if(firstImp && kN) return improve;
                } else
                    Sol_cpy(current,sol);

                _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
                if( _time > timeRem) return improve;

            }
        }
    }
    return improve;
}

int Neighbor_search_swapProj(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *current, int kN, int firstImp, double timeRem)
{
    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(current != NULL);
    assert(kN >= 0 && kN < neighborhood->nNeighborhood);

    double _time = 0;
    clock_t tStart = clock();

    int improve = 0,valid = 0;
    for(int p1 = 0; p1 < Inst_nProjects(Sol_inst(sol)); p1++) {
        for(int p2 = p1+1; p2 < Inst_nProjects(Sol_inst(sol)); p2++) {

            valid = Neighbor_swapTwoProj(neighborhood, current,p1,p2);
            if(!valid) continue;
            Sol_rebuild_opt(current, sol);

            if(Sol_getCost(current) < Sol_getCost(sol) ) {
                Sol_cpy(sol,current);
                if(Sol_getCost(current) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                    _time = ( (double)clock()/CLOCKS_PER_SEC );
                    printf("\n%ld %f", Sol_getCost(bestSol), _time);
                    fflush(stdout);
                    fflush(stdin);
#endif // DEBUG
                    Sol_cpy( bestSol, sol );
                }
                improve = 1;
                if(firstImp && kN) return improve;
            } else
                Sol_cpy(current,sol);
            _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            if( _time > timeRem) return improve;

        }
    }

    return improve;

}

int Neighbor_search_compactProj(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *current, int kN, int firstImp, double timeRem)
{
    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(current != NULL);
    assert(kN >= 0 && kN < neighborhood->nNeighborhood);

    double _time = 0;
    clock_t tStart = clock();
    int improve = 0, valid=0;
    for(int p = 0; p < Inst_nProjects(Sol_inst(sol)); p++) {
        for(double perc = 0; perc < 1; perc = perc+0.25) {
            valid = Neighbor_compactProj(neighborhood, current,perc,p);
            if(!valid) continue;
            Sol_rebuild_opt(current, sol);

            if(Sol_getCost(current) < Sol_getCost(sol) ) {
                Sol_cpy(sol,current);
                if(Sol_getCost(current) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                    _time = ( (double)clock()/CLOCKS_PER_SEC );
                    printf("\n%ld %f", Sol_getCost(bestSol), _time);
                    fflush(stdout);
                    fflush(stdin);
#endif
                    Sol_cpy( bestSol, sol );
                }
                improve = 1;

                if(firstImp && kN) return improve;
            } else
                Sol_cpy(current,sol);

            _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            if( _time > timeRem) return improve;


        }
    }
    return improve;

}

int Neighbor_search_changeOneMode(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *current, int kN, int firstImp, double timeRem)
//int Neighbor_search_changeOneMode(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *current, int kN, int *nChanges,  int **nTimesJobOnModes, int firstImp, double timeRem)
{
    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(current != NULL);
    assert(kN >= 0 && kN < neighborhood->nNeighborhood);

    double _time = 0;
    clock_t tStart = clock();
    int nJobs = Inst_nJobs( Sol_inst(sol));

    int valid = 0;
    int  improve = 0;

    for(int j = 0; j < nJobs; j++) {
        const Job* job = Inst_job(Sol_inst(sol),j);
        for(int m = 0; m < Job_nModes(job); m++) {
            const Mode * mode = Job_mode(job,m);
            if(!Mode_isFeasible(Sol_inst(sol),mode)) continue;

            /*int mm = Sol_getMode(sol,Job_index(job));
            nTimesJobOnModes[Job_index(job)][mm]++;
            */
            valid = Neighbor_changeOneMode(neighborhood, current,j,m);

            if(!valid) {
                Sol_cpy(current,sol);
                valid = 0;
                continue;
            }

            //nChanges[j]++;

            Sol_rebuild_opt(current, sol);

            //Cost newFO =  Sol_getCost(current) + nChanges[j] * Neighbor_getPenaltyChangeMode(neighborhood);

            //if(  newFO < Sol_getCost(sol) ) {
            if(  Sol_getCost(current) < Sol_getCost(sol) ) {
                Sol_cpy(sol,current);
                if(Sol_getCost(current) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                    _time = ( (double)clock()/CLOCKS_PER_SEC );
                    printf("\n%ld %f", Sol_getCost(bestSol), _time);
                    fflush(stdout);
                    fflush(stdin);
#endif
                    Sol_cpy( bestSol, sol );
                }
                improve = 1;

                if(firstImp && kN) return improve;
            } else
                Sol_cpy(current,sol);

            _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            if( _time > timeRem) return improve;

        }
    }
    return improve;
}

int Neighbor_search_changeTwoMode(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *current, int kN, int firstImp, double timeRem)
//int Neighbor_search_changeTwoMode(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *current, int kN, int *nChanges,  int **nTimesJobOnModes, int firstImp, double timeRem)
{
    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(current != NULL);
    assert(kN >= 0 && kN < neighborhood->nNeighborhood);

    double _time = 0;
    clock_t tStart = clock();
    int nJobs = Inst_nJobs( Sol_inst(sol));

    int valid = 0;

    int  improve = 0;

    for(int j1 = 0; j1 < nJobs; j1++) {
        const Job* job1 = Inst_job(Sol_inst(sol),j1);
        for(int m1 = 0; m1 < Job_nModes(job1); m1++) {
            const Mode * mode1 = Job_mode(job1,m1);
            if(!Mode_isFeasible(Sol_inst(sol),mode1)) continue;
            for(int j2 = j1+1; j2 < nJobs; j2++) {
                const Job* job2 = Inst_job(Sol_inst(sol),j2);
                for(int m2 = 0; m2 < Job_nModes(job2); m2++) {
                    const Mode * mode2 = Job_mode(job2,m2);
                    if(!Mode_isFeasible(Sol_inst(sol),mode2)) continue;

                    /*  int mm1 = Sol_getMode(sol,Job_index(job1));
                      nTimesJobOnModes[Job_index(job1)][mm1]++;
                      int mm2 = Sol_getMode(sol,Job_index(job2));
                      nTimesJobOnModes[Job_index(job2)][mm2]++;
                      */
                    valid =  Neighbor_changeTwoMode(neighborhood, current, j1, j2, m1, m2);
                    if(!valid) {
                        Sol_cpy(current,sol);
                        valid = 0;
                        continue;
                    }

                    //nChanges[j1]++;
                    //nChanges[j2]++;

                    Sol_rebuild_opt(current, sol);

                    //Cost newFO =  Sol_getCost(current) + (nChanges[j1]+ nChanges[j2]) * Neighbor_getPenaltyChangeMode(neighborhood);

                    //        if( newFO < Sol_getCost(sol) ) {
                    if( Sol_getCost(current) < Sol_getCost(sol) ) {
                        Sol_cpy(sol,current);
                        if(Sol_getCost(current) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                            _time = ( (double)clock()/CLOCKS_PER_SEC );
                            printf("\n%ld %f", Sol_getCost(bestSol), _time);
                            fflush(stdout);
                            fflush(stdin);
#endif // DEBUG
                            Sol_cpy( bestSol, sol );
                        }
                        improve = 1;
                        if(firstImp && kN) return improve;
                    } else
                        Sol_cpy(current,sol);

                    _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
                    if( _time > timeRem)  return improve;


                }
            }
        }
    }
    return improve;
}

//int Neighbor_search_changeThreeMode(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *current, int kN, int *nChanges, int **nTimesJobOnModes, int firstImp, double timeRem)
int Neighbor_search_changeThreeMode(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *current, int kN, int firstImp, double timeRem)
{

    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(current != NULL);
    assert(kN >= 0 && kN < neighborhood->nNeighborhood);

    double _time = 0;
    clock_t tStart = clock();
    int nJobs = Inst_nJobs( Sol_inst(sol));

    int valid = 0;
    int  improve = 0;

    for(int j1 = 0; j1 < nJobs; j1++) {
        const Job* job1 = Inst_job(Sol_inst(sol),j1);
        for(int m1 = 0; m1 < Job_nModes(job1); m1++) {
            const Mode * mode1 = Job_mode(job1,m1);
            if(!Mode_isFeasible(Sol_inst(sol),mode1)) continue;
            for(int j2 =  0; j2 < Job_nSucc(job1); j2++) {
                int idxjob2 = Job_succ(job1,j2);
                const Job* job2 = Inst_job(Sol_inst(sol),idxjob2);
                for(int m2 = 0; m2 < Job_nModes(job2); m2++) {
                    const Mode * mode2 = Job_mode(job2,m2);
                    if(!Mode_isFeasible(Sol_inst(sol),mode2)) continue;
                    for(int j3 = 0; j3 < Job_nSucc(job2); j3++) {
                        int idxjob3 = Job_succ(job2,j3);
                        const Job* job3 = Inst_job(Sol_inst(sol),idxjob3);
                        for(int m3 = 0; m3 < Job_nModes(job3); m3++) {
                            const Mode * mode3 = Job_mode(job3,m3);
                            if(!Mode_isFeasible(Sol_inst(sol),mode3)) continue;

                            /*  int mm1 = Sol_getMode(sol,Job_index(job1));
                              nTimesJobOnModes[Job_index(job1)][mm1]++;
                              int mm2 = Sol_getMode(sol,Job_index(job2));
                              nTimesJobOnModes[Job_index(job2)][mm2]++;
                              int mm3 = Sol_getMode(sol,Job_index(job3));
                              nTimesJobOnModes[Job_index(job3)][mm3]++;
                              */
                            valid =  Neighbor_changeThreeMode(neighborhood, current,Job_index(job1), Job_index(job2), Job_index(job3), m1, m2, m3);

                            if(!valid) {
                                Sol_cpy(current,sol);
                                continue;
                            }

                            //   nChanges[Job_index(job1)]++;
                            //nChanges[Job_index(job2)]++;
                            //nChanges[Job_index(job3)]++;

                            Sol_rebuild_opt(current, sol);

                            //Cost newFO =  Sol_getCost(current) + (nChanges[Job_index(job1)]+ nChanges[Job_index(job2)]+ nChanges[Job_index(job3)]) * Neighbor_getPenaltyChangeMode(neighborhood);
                            //if( newFO < Sol_getCost(sol) ) {
                            if( Sol_getCost(current) < Sol_getCost(sol) ) {
                                Sol_cpy(sol,current);
                                if(Sol_getCost(current) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                                    _time = ( (double)clock()/CLOCKS_PER_SEC );
                                    printf("\n%ld %f", Sol_getCost(bestSol), _time);
                                    fflush(stdout);
                                    fflush(stdin);
#endif // DEBUG
                                    Sol_cpy( bestSol, sol );
                                }
                                improve = 1;

                                if(firstImp && kN)  return improve;
                            } else
                                Sol_cpy(current,sol);

                            _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
                            if( _time > timeRem) return improve;

                        }
                    }
                }
            }
        }
    }
    return improve;
}

int Neighbor_search_changeFourMode(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *current, int kN, int firstImp, double timeRem)
//int Neighbor_search_changeFourMode(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *current, int kN, int *nChanges,  int **nTimesJobOnModes, int firstImp, double timeRem)
{

    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(current != NULL);
    assert(kN >= 0 && kN < neighborhood->nNeighborhood);

    double _time = 0;
    clock_t tStart = clock();
    int nJobs = Inst_nJobs( Sol_inst(sol));

    int valid = 0;
    int  improve = 0;
    for(int j1 = 0; j1 < nJobs; j1++) {
        const Job* job1 = Inst_job(Sol_inst(sol),j1);
        for(int m1 = 0; m1 < Job_nModes(job1); m1++) {
            const Mode * mode1 = Job_mode(job1,m1);
            if(!Mode_isFeasible(Sol_inst(sol),mode1)) continue;
            for(int j2 = 0; j2 < Job_nSucc(job1); j2++) {
                int idxjob2 = Job_succ(job1,j2);
                const Job* job2 = Inst_job(Sol_inst(sol),idxjob2);
                for(int m2 = 0; m2 < Job_nModes(job2); m2++) {
                    const Mode * mode2 = Job_mode(job2,m2);
                    if(!Mode_isFeasible(Sol_inst(sol),mode2)) continue;
                    for(int j3 = 0; j3 < Job_nSucc(job2); j3++) {
                        int idxjob3 = Job_succ(job2,j3);
                        const Job* job3 = Inst_job(Sol_inst(sol),idxjob3);
                        for(int m3 = 0; m3 <Job_nModes(job3); m3++) {
                            const Mode * mode3 = Job_mode(job3,m3);
                            if(!Mode_isFeasible(Sol_inst(sol),mode3)) continue;
                            for(int j4 = 0; j4 < Job_nSucc(job3); j4++) {
                                int idxjob4 = Job_succ(job3,j4);
                                const Job* job4 = Inst_job(Sol_inst(sol),idxjob4);
                                for(int m4 = 0; m4 < Job_nModes(job4); m4++) {
                                    const Mode * mode4 = Job_mode(job4,m4);
                                    if(!Mode_isFeasible(Sol_inst(sol),mode4)) continue;

                                    /*int mm1 = Sol_getMode(sol,Job_index(job1));
                                    nTimesJobOnModes[Job_index(job1)][mm1]++;
                                    int mm2 = Sol_getMode(sol,Job_index(job2));
                                    nTimesJobOnModes[Job_index(job2)][mm2]++;
                                    int mm3 = Sol_getMode(sol,Job_index(job3));
                                    nTimesJobOnModes[Job_index(job3)][mm3]++;
                                    int mm4 = Sol_getMode(sol,Job_index(job4));
                                    nTimesJobOnModes[Job_index(job4)][mm4]++;
                                    */

                                    valid =  Neighbor_changeFourMode(neighborhood, current, Job_index(job1), Job_index(job2), Job_index(job3), Job_index(job4), m1, m2, m3, m4);

                                    if(!valid) {
                                        Sol_cpy(current,sol);
                                        continue;
                                    }
                                    //nChanges[Job_index(job1)]++;
                                    //nChanges[Job_index(job2)]++;
                                    //nChanges[Job_index(job3)]++;
                                    //nChanges[Job_index(job4)]++;

                                    Sol_rebuild_opt(current, sol);

                                    //Cost newFO =  Sol_getCost(current) + (nChanges[Job_index(job1)]+ nChanges[Job_index(job2)]+ nChanges[Job_index(job3)] +nChanges[Job_index(job4)]) * Neighbor_getPenaltyChangeMode(neighborhood);
                                    //if( newFO < Sol_getCost(sol) ) {
                                    if( Sol_getCost(current) < Sol_getCost(sol) ) {
                                        Sol_cpy(sol,current);
                                        if(Sol_getCost(current) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                                            _time = ( (double)clock()/CLOCKS_PER_SEC );
                                            printf("\n%ld %f", Sol_getCost(bestSol), _time);
                                            fflush(stdout);
                                            fflush(stdin);
#endif // DEBUG
                                            Sol_cpy( bestSol, sol );
                                        }
                                        improve = 1;
                                        if(firstImp && kN) return improve;
                                    } else
                                        Sol_cpy(current,sol);
                                    _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
                                    if( _time > timeRem) return improve;

                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return improve;
}

int Neighbor_search_swapJobFILS(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *current, int kN, int firstImp, double timeRem)
{
    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(current != NULL);
    assert(kN >= 0 && kN < neighborhood->nNeighborhood);

    double _time = 0;
    clock_t tStart = clock();
    int nJobs = Inst_nJobs( Sol_inst(sol));
    int  valid = 0;
    int  improved = 0;

    for(int j1 = 0; j1 < nJobs; j1++) {

        valid = Neighbor_swapJobFILS(neighborhood, current, sol, bestSol, j1);

        if(valid) {
            improved = 1;
#ifdef DEBUG
            _time = ( (double)clock()/CLOCKS_PER_SEC );
            printf("\n%ld %f", Sol_getCost(bestSol), _time);
            fflush(stdout);
            fflush(stdin);
#endif // DEBUG

            if(firstImp && kN) return improved;
        } else
            Sol_cpy(current,sol);

        _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
        if( _time > timeRem) return improved;
    }

    return improved;

}

int Neighbor_search_insertJobFILS(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *current, int kN, int firstImp, double timeRem)
{

    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(current != NULL);
    assert(kN >= 0 && kN < neighborhood->nNeighborhood);

    double _time = 0;
    clock_t tStart = clock();
    int nJobs = Inst_nJobs( Sol_inst(sol));
    int  valid = 0;
    int  improved = 0;

    for(int j1 = 0; j1 < nJobs; j1++) {

        valid = Neighbor_insertJobFILS(neighborhood, current, sol, bestSol, j1);

        if(valid) {
            improved = 1;
#ifdef DEBUG
            _time = ( (double)clock()/CLOCKS_PER_SEC );
            printf("\n%ld %f", Sol_getCost(bestSol), _time);
            fflush(stdout);
            fflush(stdin);
#endif // DEBUG
            if(firstImp && kN) return improved;
        } else
            Sol_cpy(current,sol);


        _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
        if( _time > timeRem) return improved;
    }

    return improved;
}

int Neighbor_search_compactOnExtreme(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *current, int kN, int firstImp, double timeRem)
{

    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(current != NULL);
    assert(kN >= 0 && kN < neighborhood->nNeighborhood);

    double _time = 0;
    clock_t tStart = clock();
    int nJobs = Inst_nJobs( Sol_inst(sol));
    int  improve = 0, valid=0;
    for(int j = 0; j < nJobs; j++) {


        valid = Neighbor_compactOnExtreme(neighborhood, current, j);
        if(!valid) continue;
        Sol_rebuild_opt(current, sol);

        if(Sol_getCost(current) < Sol_getCost(sol) ) {
            Sol_cpy(sol,current);
            if(Sol_getCost(current) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                _time = ( (double)clock()/CLOCKS_PER_SEC );
                printf("\n%ld %f", Sol_getCost(bestSol), _time);
                fflush(stdout);
                fflush(stdin);
#endif // DEBUG
                Sol_cpy( bestSol, sol );
            }
            improve = 1;
            if(firstImp && kN) return improve;
        } else
            Sol_cpy(current,sol);
        _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
        if( _time > timeRem)  return improve;

    }
    return improve;

}

int Neighbor_search_moveProj(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *current, int kN, int firstImp, double timeRem)
{
    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(current != NULL);
    assert(kN >= 0 && kN < neighborhood->nNeighborhood);

    double _time = 0;
    clock_t tStart = clock();
    int  improve = 0, valid = 0;
    for(int p = 1; p <= Inst_nProjects(Sol_inst(sol)); p++) {
        for(int dir = 0; dir <2; dir++) {

            valid = Neighbor_moveProj(neighborhood, current,p,dir);
            if(!valid) continue;
            Sol_rebuild_opt(current, sol);

            if(Sol_getCost(current) < Sol_getCost(sol) ) {
                Sol_cpy(sol,current);
                if(Sol_getCost(current) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                    _time = ( (double)clock()/CLOCKS_PER_SEC );
                    printf("\n%ld %f", Sol_getCost(bestSol), _time);
                    fflush(stdout);
                    fflush(stdin);
#endif // DEBUG
                    Sol_cpy( bestSol, sol );
                }
                improve = 1;

                if(firstImp && kN) return improve;
            } else
                Sol_cpy(current,sol);

            _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            if( _time > timeRem)  return improve;
        }
    }
    return improve;

}


int Neighbor_callDetLS( Neighborhood *neighborhood,
                        Solution *sol,
                        Solution *bestSol,
                        Solution *current,
                        // int *nChanges,
                        // int **nTimesJobOnModes,
                        int kN,
                        int firstImp,
                        double timeRem,
                        Test *test)
{

    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(neighborhood != NULL);
    assert(current != NULL);

    int improve = 0;
    int enumNeighbor = Neighbor_getIdxAssortment(neighborhood, kN);
    clock_t tStart = clock();
    double _time = 0;

    _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
    if( _time > timeRem) goto TERMINATE;

    switch(enumNeighbor) {
        case seqInvert:
#ifdef DEBUG
            printf("\n[cDLS] seqInvert : %d", seqInvert);
            fflush(stdout);
            fflush(stdin);
#endif
            improve = Neighbor_search_inv( neighborhood, sol, bestSol, current, kN, firstImp, timeRem);

            break;
        case seqShiftJob:
#ifdef DEBUG
            printf("\n[cDLS] seqShiftJob : %d", seqShiftJob);
            fflush(stdout);
            fflush(stdin);
#endif
            improve = Neighbor_search_shiftJob( neighborhood,  sol, bestSol,current, kN, firstImp, timeRem);
            break;
        case seqSwapJob:
#ifdef DEBUG
            printf("\n[cDLS] seqSwapJob : %d", seqSwapJob);
            fflush(stdout);
            fflush(stdin);
#endif
            improve = Neighbor_search_swapJob( neighborhood, sol, bestSol, current, kN, firstImp, timeRem);
            break;
        case seqShiftProj:
#ifdef DEBUG
            printf("\n[cDLS] seqShiftProj : %d", seqShiftProj);
            fflush(stdout);
            fflush(stdin);
#endif
            improve = Neighbor_search_shiftProj( neighborhood, sol, bestSol, current, kN, firstImp, timeRem);
            break;
        case seqSwapProj:
#ifdef DEBUG
            printf("\n[cDLS] seqSwapProj : %d", seqSwapProj);
            fflush(stdout);
            fflush(stdin);
#endif
            improve = Neighbor_search_swapProj( neighborhood, sol, bestSol,  current, kN, firstImp, timeRem);
            break;
        case seqCompactProj:
#ifdef DEBUG
            printf("\n[cDLS] seqCompactProj : %d", seqCompactProj);
            fflush(stdout);
            fflush(stdin);
#endif
            improve = Neighbor_search_compactProj( neighborhood, sol, bestSol, current, kN, firstImp, timeRem);
            break;
        case changeOneMode:
#ifdef DEBUG
            printf("\n[cDLS] changeOneMode : %d", changeOneMode);
            fflush(stdout);
            fflush(stdin);
#endif
            //improve = Neighbor_search_changeOneMode( neighborhood, sol, bestSol, current, kN, nChanges, nTimesJobOnModes, firstImp, timeRem);
            improve = Neighbor_search_changeOneMode( neighborhood, sol, bestSol, current, kN, firstImp, timeRem);
            break;
        case changeTwoMode:
#ifdef DEBUG
            printf("\n[cDLS] changeTwoMode : %d", changeTwoMode);
            fflush(stdout);
            fflush(stdin);
#endif
            //improve = Neighbor_search_changeTwoMode( neighborhood, sol, bestSol,current, kN, nChanges,  nTimesJobOnModes, firstImp, timeRem);
            improve = Neighbor_search_changeTwoMode( neighborhood, sol, bestSol,current, kN, firstImp, timeRem);
            break;
        case changeThreeMode:
#ifdef DEBUG
            printf("\n[cDLS] changeThreeMode : %d", changeThreeMode);
            fflush(stdout);
            fflush(stdin);
#endif
            //improve = Neighbor_search_changeThreeMode( neighborhood, sol, bestSol,current, kN, nChanges,  nTimesJobOnModes, firstImp, timeRem);
            improve = Neighbor_search_changeThreeMode( neighborhood, sol, bestSol,current, kN, firstImp, timeRem);
            break;
        case changeFourMode:
#ifdef DEBUG
            printf("\n[cDLS] changeFourMode : %d", changeFourMode);
            fflush(stdout);
            fflush(stdin);
#endif
            //improve = Neighbor_search_changeFourMode( neighborhood, sol, bestSol,current, kN, nChanges,  nTimesJobOnModes, firstImp, timeRem);
            improve = Neighbor_search_changeFourMode( neighborhood, sol, bestSol,current, kN, firstImp, timeRem);
            break;
        case seqSwapJobFILS:
#ifdef DEBUG
            printf("\n[cDLS] seqSwapJobFILS : %d", seqSwapJobFILS);
            fflush(stdout);
            fflush(stdin);
#endif
            improve = Neighbor_search_swapJobFILS( neighborhood, sol, bestSol, current, kN, firstImp, timeRem);
            break;
        case seqInsertJobFILS:
#ifdef DEBUG
            printf("\n[cDLS] seqInsertJobFILS : %d", seqInsertJobFILS);
            fflush(stdout);
            fflush(stdin);
#endif
            improve = Neighbor_search_insertJobFILS( neighborhood, sol, bestSol, current, kN, firstImp, timeRem);
            break;
        case seqCompOnExtrem:
#ifdef DEBUG
            printf("\n[cDLS] seqCompOnExtrem : %d", seqCompOnExtrem);
            fflush(stdout);
            fflush(stdin);
#endif
            improve = Neighbor_search_compactOnExtreme( neighborhood, sol, bestSol, current,kN, firstImp,timeRem);
            break;
        case seqMoveProj:
#ifdef DEBUG
            printf("\n[cDLS] seqMoveProj : %d", seqMoveProj);
            fflush(stdout);
            fflush(stdin);
#endif
            improve = Neighbor_search_moveProj( neighborhood, sol, bestSol, current, kN, firstImp, timeRem);
            break;

    }

TERMINATE:

    Sol_cpy(sol,bestSol);

    if(improve && kN) return improve;
    else return 0;


}

int Neighbor_callDetLS_Parallel( Neighborhood *neigh,
                                 Solution *sol,
                                 Solution *bestSol,
                                 Solution *currentT[],
                                 int kN,
                                 int firstImp,
                                 double timeRem,
                                 int nThreads,
                                 Test *test)
{


    assert(neigh != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(currentT != NULL);


    int improve = 0;

    int enumNeighbor = Neighbor_getIdxAssortment( neigh, kN);
    clock_t tStart = clock();
    double _time = 0, _timeNeigh = 0;

    _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
    if( _time > timeRem) goto TERMINATE;
    switch(enumNeighbor) {
        case seqInvert:
#ifdef DEBUG
            printf("\n[cDLS] seqInvert : %d", seqInvert);
            fflush(stdout);
            fflush(stdin);
#endif
            Test_setCurrentFO(test, Sol_getCost(bestSol));
            improve = Neighbor_search_inv_parallel( neigh, sol, bestSol, currentT, kN, firstImp, timeRem,nThreads);
            _timeNeigh = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            Test_callTest(test, seqInvert, _timeNeigh, Sol_getCost(bestSol));
            break;
        case seqShiftJob:
#ifdef DEBUG
            printf("\n[cDLS] seqShiftJob : %d", seqShiftJob);
            fflush(stdout);
            fflush(stdin);
#endif
            Test_setCurrentFO(test, Sol_getCost(bestSol));
            improve = Neighbor_search_shiftJob_parallel( neigh,  sol, bestSol,currentT, kN, firstImp, timeRem,nThreads);
            _timeNeigh = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            Test_callTest(test, seqShiftJob, _timeNeigh, Sol_getCost(bestSol));
            break;
        case seqSwapJob:
#ifdef DEBUG
            printf("\n[cDLS] seqSwapJob : %d", seqSwapJob);
            fflush(stdout);
            fflush(stdin);
#endif
            Test_setCurrentFO(test, Sol_getCost(bestSol));
            improve = Neighbor_search_swapJob_parallel( neigh, sol, bestSol, currentT, kN, firstImp, timeRem,nThreads);
            _timeNeigh = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            Test_callTest(test, seqSwapJob, _timeNeigh, Sol_getCost(bestSol));
            break;
        case seqShiftProj:
#ifdef DEBUG
            printf("\n[cDLS] seqShiftProj : %d", seqShiftProj);
            fflush(stdout);
            fflush(stdin);
#endif
            Test_setCurrentFO(test, Sol_getCost(bestSol));
            improve = Neighbor_search_shiftProj_parallel( neigh, sol, bestSol, currentT, kN, firstImp, timeRem,nThreads);
            _timeNeigh = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            Test_callTest(test, seqShiftProj, _timeNeigh, Sol_getCost(bestSol));
            break;
        case seqSwapProj:
#ifdef DEBUG
            printf("\n[cDLS] seqSwapProj : %d", seqSwapProj);
            fflush(stdout);
            fflush(stdin);
#endif
            Test_setCurrentFO(test, Sol_getCost(bestSol));
            improve = Neighbor_search_swapProj_parallel( neigh, sol, bestSol,  currentT, kN, firstImp, timeRem,nThreads);
            _timeNeigh = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            Test_callTest(test, seqSwapProj, _timeNeigh, Sol_getCost(bestSol));
            break;
        case seqCompactProj:
#ifdef DEBUG
            printf("\n[cDLS] seqCompactProj : %d", seqCompactProj);
            fflush(stdout);
            fflush(stdin);
#endif
            Test_setCurrentFO(test, Sol_getCost(bestSol));
            improve = Neighbor_search_compactProj_parallel( neigh, sol, bestSol, currentT, kN, firstImp, timeRem,nThreads);
            _timeNeigh = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            Test_callTest(test, seqCompactProj, _timeNeigh, Sol_getCost(bestSol));
            break;
        case changeOneMode:
#ifdef DEBUG
            printf("\n[cDLS] changeOneMode : %d", changeOneMode);
            fflush(stdout);
            fflush(stdin);
#endif
            Test_setCurrentFO(test, Sol_getCost(bestSol));
            improve = Neighbor_search_changeOneMode_parallel( neigh, sol, bestSol, currentT, kN, firstImp, timeRem,nThreads);
            _timeNeigh = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            Test_callTest(test, changeOneMode, _timeNeigh, Sol_getCost(bestSol));
            break;
        case changeTwoMode:
#ifdef DEBUG
            printf("\n[cDLS] changeTwoMode : %d", changeTwoMode);
            fflush(stdout);
            fflush(stdin);
#endif
            fflush(stdin);
            fflush(stdout);
            Test_setCurrentFO(test, Sol_getCost(bestSol));
            improve = Neighbor_search_changeTwoMode_parallel( neigh, sol, bestSol,currentT, kN, firstImp, timeRem,nThreads);
            _timeNeigh = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            Test_callTest(test, changeTwoMode, _timeNeigh, Sol_getCost(bestSol));
            break;
        case changeThreeMode:
#ifdef DEBUG
            printf("\n[cDLS] changeThreeMode : %d", changeThreeMode);
            fflush(stdout);
            fflush(stdin);
#endif
            Test_setCurrentFO(test, Sol_getCost(bestSol));
            improve = Neighbor_search_changeThreeMode_parallel( neigh, sol, bestSol,currentT, kN, firstImp, timeRem,nThreads);
            _timeNeigh = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            Test_callTest(test, changeThreeMode, _timeNeigh, Sol_getCost(bestSol));
            break;
        case changeFourMode:
#ifdef DEBUG
            printf("\n[cDLS] changeFourMode : %d", changeFourMode);
            fflush(stdout);
            fflush(stdin);
#endif
            Test_setCurrentFO(test, Sol_getCost(bestSol));
            improve = Neighbor_search_changeFourMode_parallel( neigh, sol, bestSol,currentT, kN, firstImp, timeRem,nThreads);
            _timeNeigh = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            Test_callTest(test, changeFourMode, _timeNeigh, Sol_getCost(bestSol));
            break;

        case seqSwapJobFILS:
#ifdef DEBUG
            printf("\n[cDLS] seqSwapJobFILS : %d", seqSwapJobFILS);
            fflush(stdout);
            fflush(stdin);
#endif
            Test_setCurrentFO(test, Sol_getCost(bestSol));
            improve = Neighbor_search_swapJobFILS_parallel( neigh, sol, bestSol, currentT, kN, firstImp, timeRem,nThreads);
            _timeNeigh = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            Test_callTest(test, seqSwapJobFILS, _timeNeigh, Sol_getCost(bestSol));
            break;
        case seqInsertJobFILS:
#ifdef DEBUG
            printf("\n[cDLS] seqInsertJobFILS : %d", seqInsertJobFILS);
            fflush(stdout);
            fflush(stdin);
#endif
            Test_setCurrentFO(test, Sol_getCost(bestSol));
            improve = Neighbor_search_insertJobFILS_parallel( neigh, sol, bestSol, currentT, kN, firstImp, timeRem, nThreads);
            _timeNeigh = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            Test_callTest(test, seqInsertJobFILS, _timeNeigh, Sol_getCost(bestSol));
            break;
        case seqCompOnExtrem:
#ifdef DEBUG
            printf("\n[cDLS] seqCompOnExtrem : %d", seqCompOnExtrem);
            fflush(stdout);
            fflush(stdin);
#endif
            Test_setCurrentFO(test, Sol_getCost(bestSol));
            improve = Neighbor_search_compactOnExtreme_parallel( neigh, sol, bestSol, currentT,kN, firstImp,timeRem,nThreads);
            _timeNeigh = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            Test_callTest(test, seqCompOnExtrem, _timeNeigh, Sol_getCost(bestSol));
            break;
        case seqMoveProj:
#ifdef DEBUG
            printf("\n[cDLS] seqMoveProj : %d", seqMoveProj);
            fflush(stdout);
            fflush(stdin);
#endif
            Test_setCurrentFO(test, Sol_getCost(bestSol));
            improve = Neighbor_search_moveProj_parallel( neigh, sol, bestSol, currentT, kN, firstImp, timeRem,nThreads);
            _timeNeigh = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            Test_callTest(test, seqMoveProj, _timeNeigh, Sol_getCost(bestSol));
            break;

    }

TERMINATE:

    Sol_cpy(sol,bestSol);

    if(improve && kN) return improve;
    else return 0;


}

//call deterministic in parallel
int Neighbor_search_inv_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads)
{

    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(currentT != NULL);

    const double tStart = omp_get_wtime();
    int nJobs = Inst_nJobs( Sol_inst(sol));
    int improve = 0, abort = 0;

    int vec[ nJobs*(nJobs-1) ][2];
    int l=0;
    for(int k = Neighbor_getMinK(neighborhood,seqInvert-1); k < Neighbor_getMaxK( neighborhood, seqInvert-1); k++) {
        for(int  i = 0; i < nJobs-k; i++) {
            vec[l][0] = i;
            vec[l][1] = k;
            l++;
        }
    }

    omp_set_num_threads(nTreads);


    #pragma omp parallel for shared (bestSol,sol)
    for(int j =0 ; j < l; j++) {

        #pragma omp flush (abort)
        if (!abort) {

            int th_id, valid =0;

            th_id = omp_get_thread_num();

            valid = Neighbor_inv(currentT[th_id], vec[j][0], vec[j][1]);
            if(!valid) continue;
            //Sol_rebuild_opt(currentT[th_id], sol);
            Sol_rebuild(currentT[th_id]);

            #pragma omp critical
            {

                if(Sol_getCost(currentT[th_id]) < Sol_getCost(sol) ) {
                    Sol_cpy(sol,currentT[th_id]);
                    if(Sol_getCost(currentT[th_id]) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                        th_id = omp_get_thread_num();
                        printf("\n%ld thread %d", Sol_getCost(bestSol), th_id);
#endif
                        Sol_cpy( bestSol, currentT[th_id] );
                    }
                    improve = 1;
                    if(firstImp && kN) {
                        abort = 1;
                        #pragma omp flush (abort)
                    }
                } else
                    Sol_cpy(currentT[th_id], sol);

                double _time  = omp_get_wtime()-tStart;

                if( _time > timeRem) {
                    abort = 1;
                    #pragma omp flush (abort)
                }

            }
        }
    }

    return improve;

}

int Neighbor_search_shiftJob_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *currentT[],int kN, int firstImp, double timeRem, int nTreads)
{
    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(currentT != NULL);

    double _time = 0;
    clock_t tStart = clock();
    int nJobs = Inst_nJobs( Sol_inst(sol));
    int improve = 0,  abort = 0;



    int vec[ nJobs*(nJobs-1) ][3];

    int l=0;
    for(int k = Neighbor_getMinK(neighborhood,seqShiftJob-1); k < Neighbor_getMaxK(neighborhood,seqShiftJob-1); k++) {
        for(int dir = -1; dir < 1; dir++ ) {
            for(int i = 0; i < nJobs-k; i++) {
                vec[l][0] = i;
                vec[l][1] = k;
                vec[l][2] = dir;
                l++;
            }
        }
    }

    omp_set_num_threads(nTreads);

    #pragma omp parallel for shared (bestSol,sol)
    for(int j =0 ; j < l; j++) {

        #pragma omp flush (abort)
        if (!abort) {

            int th_id, valid =0;

            th_id = omp_get_thread_num();

            valid= Neighbor_shiftJob( neighborhood, currentT[th_id], vec[j][0], vec[j][1], vec[j][2]);
            if(!valid) continue;
            //Sol_rebuild_opt(currentT[th_id], sol);
            Sol_rebuild(currentT[th_id]);

            #pragma omp critical
            {
                if(Sol_getCost(currentT[th_id]) < Sol_getCost(sol) ) {
                    Sol_cpy(sol,currentT[th_id]);
                    if(Sol_getCost(currentT[th_id]) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                        th_id = omp_get_thread_num();
                        printf("\n%ld thread %d", Sol_getCost(bestSol), th_id);
#endif
                        Sol_cpy( bestSol, sol );
                    }
                    improve = 1;

                    if(firstImp && kN) {
                        abort = 1;
                        #pragma omp flush (abort)
                    }
                } else
                    Sol_cpy(currentT[th_id], sol);
                _time = omp_get_wtime()-tStart;
                if( _time > timeRem) {
                    abort = 1;
                    #pragma omp flush (abort)
                }

            }
        }
    }
    return improve;

}

int Neighbor_search_swapJob_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads)
{

    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(currentT != NULL);

    double _time = 0;
    clock_t tStart = clock();
    int nJobs = Inst_nJobs( Sol_inst(sol));
    int  improve = 0, abort = 0;

    int vecJobs[nJobs*(nJobs-1)][2];
    int l=0;
    for(int j1 = 0; j1 < nJobs; j1++) {
        for(int j2 = j1+1; j2 < nJobs; j2++) {
            vecJobs[l][0] = j1;
            vecJobs[l][1] = j2;
            l++;

        }
    }

    omp_set_num_threads(nTreads);

    #pragma omp parallel for shared (bestSol,sol)
    for(int v = 0; v < l; v++) {
        #pragma omp flush (abort)
        if (!abort) {

            int th_id, valid =0;

            th_id = omp_get_thread_num();

            valid = Neighbor_swapJob(neighborhood, currentT[th_id], vecJobs[v][0], vecJobs[v][1]);
            if(!valid) continue;
            //Sol_rebuild_opt(currentT[th_id], sol);
            Sol_rebuild(currentT[th_id]);

            #pragma omp critical
            {

                if(Sol_getCost(currentT[th_id]) < Sol_getCost(sol) ) {
                    Sol_cpy(sol,currentT[th_id]);
                    if(Sol_getCost(currentT[th_id]) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                        th_id = omp_get_thread_num();
                        printf("\n%ld thread %d", Sol_getCost(bestSol), th_id);
#endif
                        Sol_cpy( bestSol, sol );
                    }
                    improve = 1;

                    if(firstImp && kN) {
                        abort = 1;
                        #pragma omp flush (abort)
                    }
                } else
                    Sol_cpy(currentT[th_id], sol);

                _time = omp_get_wtime()-tStart;
                if( _time > timeRem) {
                    abort = 1;
                    #pragma omp flush (abort)
                }
            }
        }
    }

    return improve;
}

int Neighbor_search_shiftProj_parallel(Neighborhood *neighborhood,Solution *sol, Solution *bestSol, Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads)
{

    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(currentT != NULL);

    double _time = 0;
    clock_t tStart = clock();
    int improve = 0, abort = 0;
    int nJobs = Inst_nJobs( Sol_inst(sol));

    int vec[ nJobs*(nJobs-1) ][3];

    int l=0;
    for(int k = Neighbor_getMinK(neighborhood, seqShiftProj-1); k<Neighbor_getMaxK(neighborhood,seqShiftProj-1); k++) {
        for(int dir = -1; dir < 1; dir++ ) {
            for(int p = 0; p < Inst_nProjects(Sol_inst(sol)); p++) {
                vec[l][0] = k;
                vec[l][1] = p;
                vec[l][2] = dir;
                l++;
            }
        }
    }

    omp_set_num_threads(nTreads);

    #pragma omp parallel for shared (bestSol,sol)
    for(int j =0 ; j < l; j++) {
        #pragma omp flush (abort)
        if (!abort) {

            int th_id, valid =0;

            th_id = omp_get_thread_num();

            valid = Neighbor_shiftProj(neighborhood,  currentT[th_id],vec[j][0],vec[j][1],vec[j][2]);
            if(!valid) continue;
            //Sol_rebuild_opt(currentT[th_id], sol);
            Sol_rebuild(currentT[th_id]);

            #pragma omp critical
            {
                if(Sol_getCost(currentT[th_id]) < Sol_getCost(sol) ) {
                    Sol_cpy(sol,currentT[th_id]);
                    if(Sol_getCost(currentT[th_id]) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                        th_id = omp_get_thread_num();
                        printf("\n%ld thread %d", Sol_getCost(bestSol), th_id);
#endif
                        Sol_cpy( bestSol, sol );
                    }
                    improve = 1;

                    if(firstImp && kN) {
                        abort = 1;
                        #pragma omp flush (abort)
                    }
                } else
                    Sol_cpy(currentT[th_id], sol);

                _time = omp_get_wtime()-tStart;
                if( _time > timeRem) {
                    abort = 1;
                    #pragma omp flush (abort)
                }
            }
        }
    }

    return improve;
}

int Neighbor_search_swapProj_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads)
{

    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(currentT != NULL);
    assert(kN >= 0 && kN < neighborhood->nNeighborhood);

    double _time = 0;
    clock_t tStart = clock();
    int nJobs = Inst_nJobs( Sol_inst(sol));
    int improve = 0, abort = 0;

    int vec[ nJobs*(nJobs-1) ][2];
    int l=0;
    for(int p1 = 0; p1 < Inst_nProjects(Sol_inst(sol)); p1++) {
        for(int p2 = p1+1; p2 < Inst_nProjects(Sol_inst(sol)); p2++) {
            vec[l][0] = p1;
            vec[l][1] = p2;
            l++;
        }
    }

    omp_set_num_threads(nTreads);
    #pragma omp parallel for shared (bestSol,sol)
    for(int j =0 ; j < l; j++) {
        #pragma omp flush (abort)
        if (!abort) {

            int th_id, valid =0;

            th_id = omp_get_thread_num();

            valid = Neighbor_swapTwoProj(neighborhood, currentT[th_id],vec[j][0],vec[j][1]);
            if(!valid) continue;
            //Sol_rebuild_opt(currentT[th_id], sol);
            Sol_rebuild(currentT[th_id]);

            #pragma omp critical
            {

                if(Sol_getCost(currentT[th_id]) < Sol_getCost(sol) ) {
                    Sol_cpy(sol,currentT[th_id]);
                    if(Sol_getCost(currentT[th_id]) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                        th_id = omp_get_thread_num();
                        printf("\n%ld thread %d", Sol_getCost(bestSol), th_id);
#endif
                        Sol_cpy( bestSol, sol );
                    }
                    improve = 1;
                    if(firstImp && kN) {
                        abort = 1;
                        #pragma omp flush (abort)
                    }
                } else
                    Sol_cpy(currentT[th_id], sol);
                _time = omp_get_wtime()-tStart;
                if( _time > timeRem) {
                    abort = 1;
                    #pragma omp flush (abort)
                }
            }
        }
    }

    return improve;

}

int Neighbor_search_compactProj_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads)
{
    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(currentT != NULL);


    double _time = 0;
    clock_t tStart = clock();
    int nJobs = Inst_nJobs( Sol_inst(sol));

    int improve = 0, abort = 0;

    double vec[ nJobs*(nJobs-1) ][2];
    int l=0;
    for(int p = 0; p < Inst_nProjects(Sol_inst(sol)); p++) {
        for(double perc = 0; perc < 1; perc = perc+0.25) {
            vec[l][0] = perc;
            vec[l][1] = p;
            l++;
        }
    }

    omp_set_num_threads(nTreads);
    #pragma omp parallel for shared (bestSol,sol)
    for(int j =0 ; j < l; j++) {
        #pragma omp flush (abort)
        if (!abort) {

            int th_id, valid =0;

            th_id = omp_get_thread_num();

            valid = Neighbor_compactProj(neighborhood, currentT[th_id], vec[j][0],  (int)vec[j][1] );
            if(!valid) continue;
            //Sol_rebuild_opt(currentT[th_id], sol);
            Sol_rebuild(currentT[th_id]);

            #pragma omp critical
            {
                if(Sol_getCost(currentT[th_id]) < Sol_getCost(sol) ) {
                    Sol_cpy(sol,currentT[th_id]);
                    if(Sol_getCost(currentT[th_id]) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                        th_id = omp_get_thread_num();
                        printf("\n%ld thread %d", Sol_getCost(bestSol), th_id);
#endif
                        Sol_cpy( bestSol, sol );
                    }
                    improve = 1;

                    if(firstImp && kN) {
                        abort = 1;
                        #pragma omp flush (abort)
                    }
                } else
                    Sol_cpy(currentT[th_id], sol);

                _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
                if( _time > timeRem) {
                    abort = 1;
                    #pragma omp flush (abort)
                }
            }

        }
    }
    return improve;

}

int Neighbor_search_changeOneMode_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads)
{

    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(currentT != NULL);

    double _time = 0;
    clock_t tStart = clock();
    int nJobs = Inst_nJobs( Sol_inst(sol));
    int improve = 0, abort = 0;
    int nComb = nJobs*3;
    int **vec = (int **) calloc (nComb, sizeof(int *));
    for ( int i = 0; i < nComb; i++ )
        vec[i] = (int*) calloc (2, sizeof(int));



    int l=0;
    for(int j = 0; j < nJobs; j++) {
        const Job* job = Inst_job(Sol_inst(sol),j);
        for(int m = 0; m < Job_nModes(job); m++) {
            const Mode * mode = Job_mode(job,m);
            if(!Mode_isFeasible(Sol_inst(sol),mode)) continue;
            vec[l][0] = Job_index(job);
            vec[l][1] = m;
            l++;
        }
    }

    omp_set_num_threads(nTreads);
    #pragma omp parallel for shared (bestSol,sol)
    for(int j =0 ; j < l; j++) {
        #pragma omp flush (abort)
        if (!abort) {

            int th_id, valid =0;

            th_id = omp_get_thread_num();

            valid = Neighbor_changeOneMode(neighborhood, currentT[th_id],vec[j][0],vec[j][1]);
            if(!valid) continue;
            //Sol_rebuild_opt(currentT[th_id], sol);
            Sol_rebuild(currentT[th_id]);

            #pragma omp critical
            {
                if(  Sol_getCost(currentT[th_id])  < Sol_getCost(sol) ) {
                    Sol_cpy(sol,currentT[th_id]);
                    if(Sol_getCost(currentT[th_id]) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                        th_id = omp_get_thread_num();
                        printf("\n%ld thread %d", Sol_getCost(bestSol), th_id);
#endif
                        Sol_cpy( bestSol, sol );
                    }
                    improve = 1;
                    if(firstImp && kN) {
                        abort = 1;
                        #pragma omp flush (abort)
                    }
                } else
                    Sol_cpy(currentT[th_id], sol);
                _time = omp_get_wtime()-tStart;
                if( _time > timeRem) {
                    abort = 1;
                    #pragma omp flush (abort)
                }
            }
        }

    }
    return improve;
}

int Neighbor_search_changeTwoMode_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *currentT[], int kN,  int firstImp, double timeRem, int nTreads)
{
    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(currentT != NULL);

    double _time = 0;
    clock_t tStart = clock();
    int nJobs = Inst_nJobs( Sol_inst(sol));

    int  improve = 0, abort = 0;

    int nComb = nJobs*3*(nJobs-1)*3;
    int l=0;

    int **vec = (int **) calloc (nComb, sizeof(int *));
    for ( int i = 0; i < nComb; i++ )
        vec[i] = (int*) calloc (4, sizeof(int));


    for(int j1 = 0; j1 < nJobs; j1++) {

        const Job* job1 = Inst_job(Sol_inst(sol),j1);

        for(int m1 = 0; m1 < Job_nModes(job1); m1++) {
            const Mode * mode1 = Job_mode(job1,m1);
            if(!Mode_isFeasible(Sol_inst(sol),mode1)) continue;
            for(int j2 = j1+1; j2 < nJobs; j2++) {
                const Job* job2 = Inst_job(Sol_inst(sol),j2);
                for(int m2 = 0; m2 < Job_nModes(job2); m2++) {
                    const Mode * mode2 = Job_mode(job2,m2);
                    if(!Mode_isFeasible(Sol_inst(sol),mode2)) continue;
                    vec[l][0] = j1;
                    vec[l][1] = j2;
                    vec[l][2] = m1;
                    vec[l][3] = m2;
                    l++;

                }
            }
        }
    }

    omp_set_num_threads(nTreads);
    #pragma omp parallel for shared (bestSol,sol)
    for(int j =0 ; j < l; j++) {
        #pragma omp flush (abort)
        if (!abort) {

            int th_id, valid =0;

            th_id = omp_get_thread_num();

            valid =  Neighbor_changeTwoMode(neighborhood, currentT[th_id], vec[j][0], vec[j][1], vec[j][2],vec[j][3]);

            if(!valid)   continue;
            //Sol_rebuild_opt(currentT[th_id], sol);
            Sol_rebuild(currentT[th_id]);

            #pragma omp critical
            {


                if( Sol_getCost(currentT[th_id])  < Sol_getCost(sol) ) {
                    Sol_cpy(sol,currentT[th_id]);
                    if(Sol_getCost(currentT[th_id]) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                        th_id = omp_get_thread_num();
                        printf("\n%ld thread %d", Sol_getCost(bestSol), th_id);
#endif
                        Sol_cpy( bestSol, sol );
                    }
                    improve = 1;
                    if(firstImp && kN) {
                        abort = 1;
                        #pragma omp flush (abort)
                    }
                } else
                    Sol_cpy(currentT[th_id], sol);
                _time = omp_get_wtime()-tStart;
                if( _time > timeRem) {
                    abort = 1;
                    #pragma omp flush (abort)
                }
            }
        }
    }

    return improve;
}

int Neighbor_search_changeThreeMode_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *currentT[],int kN,  int firstImp, double timeRem, int nTreads)
{

    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(currentT != NULL);

    double _time = 0;
    clock_t tStart = clock();
    int nJobs = Inst_nJobs( Sol_inst(sol));
    int  improve = 0, abort = 0;
    int nComb = nJobs*3*3*3*3*3;
    int **vec = (int **) calloc (nComb, sizeof(int *));
    for ( int i = 0; i < nComb; i++ )
        vec[i] = (int*) calloc (6, sizeof(int));

    int l=0;

    for(int j1 = 0; j1 < nJobs; j1++) {
        const Job* job1 = Inst_job(Sol_inst(sol),j1);
        for(int m1 = 0; m1 < Job_nModes(job1); m1++) {
            const Mode * mode1 = Job_mode(job1,m1);
            if(!Mode_isFeasible(Sol_inst(sol),mode1)) continue;
            for(int j2 =  0; j2 < Job_nSucc(job1); j2++) {
                int idxjob2 = Job_succ(job1,j2);
                const Job* job2 = Inst_job(Sol_inst(sol),idxjob2);
                for(int m2 = 0; m2 < Job_nModes(job2); m2++) {
                    const Mode * mode2 = Job_mode(job2,m2);
                    if(!Mode_isFeasible(Sol_inst(sol),mode2)) continue;
                    for(int j3 = 0; j3 < Job_nSucc(job2); j3++) {
                        int idxjob3 = Job_succ(job2,j3);
                        const Job* job3 = Inst_job(Sol_inst(sol),idxjob3);
                        for(int m3 = 0; m3 < Job_nModes(job3); m3++) {
                            const Mode * mode3 = Job_mode(job3,m3);
                            if(!Mode_isFeasible(Sol_inst(sol),mode3)) continue;
                            vec[l][0] = j1;
                            vec[l][1] = idxjob2;
                            vec[l][2] = idxjob3;
                            vec[l][3] = m1;
                            vec[l][4] = m2;
                            vec[l][5] = m3;
                            l++;
                        }
                    }
                }
            }
        }
    }

    omp_set_num_threads(nTreads);
    #pragma omp parallel for shared (bestSol,sol)
    for(int j =0 ; j < l; j++) {
        #pragma omp flush (abort)
        if (!abort) {

            int th_id, valid =0;

            th_id = omp_get_thread_num();

            valid =  Neighbor_changeThreeMode(neighborhood, currentT[th_id], vec[j][0], vec[j][1], vec[j][2], vec[j][3], vec[j][4], vec[j][5]);
            if(!valid)  continue;
            //Sol_rebuild_opt(currentT[th_id], sol);
            Sol_rebuild(currentT[th_id]);

            #pragma omp critical
            {

                if( Sol_getCost(currentT[th_id])  < Sol_getCost(sol) ) {
                    Sol_cpy(sol,currentT[th_id]);
                    if(Sol_getCost(currentT[th_id]) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                        th_id = omp_get_thread_num();
                        printf("\n%ld thread %d", Sol_getCost(bestSol), th_id);
#endif
                        Sol_cpy( bestSol, sol );
                    }
                    improve = 1;

                    if(firstImp && kN) {
                        abort = 1;
                        #pragma omp flush (abort)
                    }
                } else
                    Sol_cpy(currentT[th_id], sol);
                _time = omp_get_wtime()-tStart;
                if( _time > timeRem) {
                    abort = 1;
                    #pragma omp flush (abort)
                }
            }
        }
    }

    return improve;
}

int Neighbor_search_changeFourMode_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads)
{
    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(currentT != NULL);


    double _time = 0;
    clock_t tStart = clock();
    int nJobs = Inst_nJobs( Sol_inst(sol));
    int improve = 0, abort = 0;
    int nComb = nJobs*3*3*3*3*3*3*3;

    int **vec = (int **) calloc (nComb, sizeof(int *));
    for ( int i = 0; i < nComb; i++ )
        vec[i] = (int*) calloc (8, sizeof(int));

    int l=0;

    for(int j1 = 0; j1 < nJobs; j1++) {
        const Job* job1 = Inst_job(Sol_inst(sol),j1);
        for(int m1 = 0; m1 < Job_nModes(job1); m1++) {
            const Mode * mode1 = Job_mode(job1,m1);
            if(!Mode_isFeasible(Sol_inst(sol),mode1)) continue;
            for(int j2 = 0; j2 < Job_nSucc(job1); j2++) {
                int idxjob2 = Job_succ(job1,j2);
                const Job* job2 = Inst_job(Sol_inst(sol),idxjob2);
                for(int m2 = 0; m2 < Job_nModes(job2); m2++) {
                    const Mode * mode2 = Job_mode(job2,m2);
                    if(!Mode_isFeasible(Sol_inst(sol),mode2)) continue;
                    for(int j3 = 0; j3 < Job_nSucc(job2); j3++) {
                        int idxjob3 = Job_succ(job2,j3);
                        const Job* job3 = Inst_job(Sol_inst(sol),idxjob3);
                        for(int m3 = 0; m3 <Job_nModes(job3); m3++) {
                            const Mode * mode3 = Job_mode(job3,m3);
                            if(!Mode_isFeasible(Sol_inst(sol),mode3)) continue;
                            for(int j4 = 0; j4 < Job_nSucc(job3); j4++) {
                                int idxjob4 = Job_succ(job3,j4);
                                const Job* job4 = Inst_job(Sol_inst(sol),idxjob4);
                                for(int m4 = 0; m4 < Job_nModes(job4); m4++) {
                                    const Mode * mode4 = Job_mode(job4,m4);
                                    if(!Mode_isFeasible(Sol_inst(sol),mode4)) continue;
                                    vec[l][0] = j1;
                                    vec[l][1] = idxjob2;
                                    vec[l][2] = idxjob3;
                                    vec[l][3] = idxjob4;
                                    vec[l][4] = m1;
                                    vec[l][5] = m2;
                                    vec[l][6] = m3;
                                    vec[l][7] = m4;
                                    l++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    omp_set_num_threads(nTreads);
    #pragma omp parallel for shared (bestSol,sol)
    for(int j =0 ; j < l; j++) {
        #pragma omp flush (abort)
        if (!abort) {

            int th_id, valid =0;

            th_id = omp_get_thread_num();

            valid =  Neighbor_changeFourMode(neighborhood, currentT[th_id], vec[j][0], vec[j][1], vec[j][2], vec[j][3], vec[j][4], vec[j][5], vec[j][6], vec[j][7]);
            if(!valid) continue;
            //Sol_rebuild_opt(currentT[th_id], sol);
            Sol_rebuild(currentT[th_id]);

            #pragma omp critical
            {

                if( Sol_getCost(currentT[th_id])  < Sol_getCost(sol) ) {
                    Sol_cpy(sol,currentT[th_id]);
                    if(Sol_getCost(currentT[th_id]) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                        th_id = omp_get_thread_num();
                        printf("\n%ld thread %d", Sol_getCost(bestSol), th_id);
#endif
                        Sol_cpy( bestSol, sol );
                    }
                    improve = 1;
                    if(firstImp && kN) {
                        abort = 1;
                        #pragma omp flush (abort)
                    }
                } else
                    Sol_cpy(currentT[th_id], sol);
                _time = omp_get_wtime()-tStart;
                if( _time > timeRem) {
                    abort = 1;
                    #pragma omp flush (abort)
                }
            }
        }

    }
    return improve;
}

int Neighbor_swapJobFILS_parallel(Neighborhood *neighborhood, Solution *current, Solution *sol, Solution *bestSol, int j1 )
{

    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(current != NULL);

    int l = 0, tamMaxW=0, posWstart=0, posWend=0, iniW=0, aux;
    int nJobs = Inst_nJobs(Sol_inst(current));
    int *sequence = Sol_sequence(current);

    assert(j1 >=0 && j1 < nJobs);

    const Job* job1 = Inst_job(Sol_inst(current), sequence[j1]);
    const Job* job2;

    if(Job_nPred(job1)  == 0)
        posWstart = j1;
    else {

        for(int i = j1; i >= 0; --i) {
            job2 = Inst_job(Sol_inst(current), sequence[i]);
            if(Job_hasPred(Sol_inst(current), job1, Job_index(job2))) {
                posWstart = i;
                break;
            }
        }
    }

    if(Job_nSucc(job1) == 0)
        posWend = j1;
    else {

        for(int i = j1; i < nJobs; ++i) {
            job2 = Inst_job( Sol_inst(current), sequence[i]);
            if(Job_hasSucc(Sol_inst(current),job1,  Job_index(job2))) {
                posWend = i;
                break;
            }
        }
    }

    tamMaxW = posWend - posWstart - 1;

    if(tamMaxW < 3)
        return 0;

    l = tamMaxW/2;

    iniW = tamMaxW-l;

    posWstart = posWstart + iniW + 1;
    posWend = posWstart + l - 1;

    int valid = 0;

    for(int i = posWstart; i <= posWend; ++i) {

        if(i<j1)
            valid = Neighbor_verifyPredSwapJob( current, i, j1 );
        else
            valid = Neighbor_verifyPredSwapJob( current, j1, i );

        if(valid) {

            aux = sequence[i];
            sequence[i] = sequence[j1];
            sequence[j1] = aux;

            //Sol_rebuild_opt(current, sol);
            Sol_rebuild(current);
            #pragma omp critical
            {
                if(Sol_getCost(current) < Sol_getCost(bestSol)) {
                    Sol_cpy(sol, current);
                    // printf("\nImprovement swapJobFILS %ld to %ld \n ", Sol_getCost(bestSol), Sol_getCost(current) );
                    Sol_cpy( bestSol, current );
                    i=posWend+1;
                } else
                    Sol_cpy(current,sol);
            }
        }
    }
    return 0;

}

int Neighbor_insertJobFILS_parallel(Neighborhood *neighborhood, Solution *current, Solution *sol, Solution *bestSol, int j1 )
{

    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(current != NULL);

    int l = 0, tamMaxW=0, posWstart=0, posWend=0, iniW=0, aux;
    int nJobs = Inst_nJobs(Sol_inst(current));

    assert(j1 >=0 && j1 < nJobs);

    int *sequence = Sol_sequence(current);

    const Job* job1 = Inst_job(Sol_inst(current), sequence[j1]);
    const Job* job2;

    if(Job_nPred(job1) == 0)
        posWstart = j1;
    else {

        for(int i = j1; i >= 0; --i) {
            job2 = Inst_job(Sol_inst(current), sequence[i]);
            if(Job_hasPred(Sol_inst(current), job1, Job_index(job2))) {
                posWstart = i;
                break;
            }
        }
    }

    if(Job_nSucc(job1) == 0)
        posWend = j1;
    else {

        for(int i = j1; i < nJobs; ++i) {
            job2 = Inst_job( Sol_inst(current), sequence[i]);
            if(Job_hasSucc(Sol_inst(current),job1,  Job_index(job2))) {
                posWend = i;
                break;
            }
        }
    }

    tamMaxW = posWend - posWstart - 1;

    if(tamMaxW < 3)
        return 0;

    l = tamMaxW/2;


    iniW = tamMaxW-l;

    posWstart = posWstart + iniW + 1;
    posWend = posWstart + l - 1;

    int valid = 0;

    for(int i = posWstart; i <= posWend; ++i) {

        if(i<j1)
            valid =  Neighbor_verifyPredShiftJobBack( current, j1, i );
        else
            valid = Neighbor_verifyPredShiftJobAhead( current, j1, i );

        if(valid) {

            aux = sequence[j1];

            if(i < j1)
                for(int j = j1; j > i; --j)
                    sequence[j] = sequence[j - 1];
            else if(i > j1)
                for(int j = j1; j < i; ++j)
                    sequence[j] = sequence[j + 1];
            else
                continue;

            sequence[i] = aux;

            //Sol_rebuild_opt(current, sol);
            Sol_rebuild(current);

            #pragma omp critical
            {
                if(Sol_getCost(current) < Sol_getCost(bestSol)) {
                    Sol_cpy(sol, current);
                    // printf("\nImprovement insertJobFILS %ld to %ld \n ", Sol_getCost(bestSol), Sol_getCost(current) );
                    Sol_cpy( bestSol, current );
                    i = posWend+ 1;
                } else
                    Sol_cpy(current,sol);
            }
        }
    }

    return 0;

}

int Neighbor_search_swapJobFILS_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads)
{

    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(currentT != NULL);

    double _time = 0;
    double tStart = omp_get_wtime();
    int nJobs = Inst_nJobs( Sol_inst(sol));
    int improved = 0,abort = 0;


    omp_set_num_threads(nTreads);
    #pragma omp parallel for shared (bestSol,sol)
    for(int j1 = 0; j1 < nJobs ; j1++) {
        #pragma omp flush (abort)
        if (!abort) {

            int th_id, valid =0;

            th_id = omp_get_thread_num();

            Sol_cpy(currentT[th_id],sol);

            valid = Neighbor_swapJobFILS_parallel(neighborhood, currentT[th_id], sol, bestSol, j1);
            if(valid) {
                improved = 1;
#ifdef DEBUG
                th_id = omp_get_thread_num();
                printf("\n%ld thread %d", Sol_getCost(bestSol), th_id);
#endif
                if(firstImp && kN) {
                    abort = 1;
                    #pragma omp flush (abort)
                }
            }
            _time = omp_get_wtime()-tStart;
            if( _time > timeRem) {
                abort = 1;
                #pragma omp flush (abort)
            }
        }

    }

    return improved;

}

int Neighbor_search_insertJobFILS_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads)
{


    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(currentT != NULL);

    double _time = 0;
    double tStart = omp_get_wtime();
    int nJobs = Inst_nJobs( Sol_inst(sol));
    int  improved = 0, abort = 0;

    omp_set_num_threads(nTreads);
    #pragma omp parallel for shared (bestSol,sol)
    for(int j1 = 0; j1 < nJobs ; j1++) {
        #pragma omp flush (abort)
        if (!abort) {

            int th_id, valid =0;

            th_id = omp_get_thread_num();

            Sol_cpy(currentT[th_id],sol);

            valid = Neighbor_insertJobFILS_parallel(neighborhood, currentT[th_id], sol, bestSol, j1);

            if(valid) {
                improved = 1;
#ifdef DEBUG
                th_id = omp_get_thread_num();
                printf("\n%ld thread %d", Sol_getCost(bestSol), th_id);
#endif
                if(firstImp && kN) {
                    abort = 1;
                    #pragma omp flush (abort)
                }
            }

            _time = omp_get_wtime()-tStart;
            if( _time > timeRem) {
                abort = 1;
                #pragma omp flush (abort)
            }
        }

    }

    return improved;
}

int Neighbor_search_compactOnExtreme_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,Solution *currentT[],  int kN, int firstImp, double timeRem, int nTreads)
{


    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(currentT != NULL);

    double _time = 0;
    double tStart = omp_get_wtime();
    int nJobs = Inst_nJobs( Sol_inst(sol));
    int  improve = 0, abort = 0;;

    omp_set_num_threads(nTreads);
    #pragma omp parallel for shared (bestSol,sol)
    for(int j = 0; j < nJobs; j++) {
        #pragma omp flush (abort)
        if (!abort) {

            int th_id, valid =0;


            th_id = omp_get_thread_num();

            valid = Neighbor_compactOnExtreme(neighborhood, currentT[th_id], j);
            if(!valid) continue;
            //Sol_rebuild_opt(currentT[th_id], sol);
            Sol_rebuild(currentT[th_id]);

            #pragma omp critical
            {

                if( Sol_getCost(currentT[th_id])  < Sol_getCost(sol) ) {
                    Sol_cpy(sol,currentT[th_id]);
                    if(Sol_getCost(currentT[th_id]) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                        th_id = omp_get_thread_num();
                        printf("\n%ld thread %d", Sol_getCost(bestSol), th_id);
#endif
                        Sol_cpy( bestSol, sol );
                    }
                    improve = 1;
                    if(firstImp && kN) {
                        abort = 1;
                        #pragma omp flush (abort)
                    }
                } else
                    Sol_cpy(currentT[th_id], sol);
                _time = omp_get_wtime()-tStart;
                if( _time > timeRem) {
                    abort = 1;
                    #pragma omp flush (abort)
                }
            }
        }

    }
    return improve;

}

int Neighbor_search_moveProj_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads)
{


    assert(neighborhood != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);
    assert(currentT != NULL);

    double _time = 0;
    double tStart = omp_get_wtime();
    int  improve = 0, abort = 0;
    int nJobs = Inst_nJobs( Sol_inst(sol));

    int vec[ nJobs*(nJobs-1) ][2];
    int l=0;
    for(int p = 1; p <= Inst_nProjects(Sol_inst(sol)); p++) {
        for(int dir = 0; dir <2; dir++) {
            vec[l][0] = p;
            vec[l][1] = dir;
            l++;
        }
    }

    omp_set_num_threads(nTreads);
    #pragma omp parallel for shared (bestSol,sol)
    for(int j =0 ; j < l; j++) {
        #pragma omp flush (abort)
        if (!abort) {

            int th_id, valid =0;

            th_id = omp_get_thread_num();


            valid= Neighbor_moveProj(neighborhood, currentT[th_id], vec[j][0], vec[j][1] );
            if(!valid) continue;
            //Sol_rebuild_opt(currentT[th_id], sol);
            Sol_rebuild(currentT[th_id]);

            #pragma omp critical
            {

                if(Sol_getCost(currentT[th_id]) < Sol_getCost(sol) ) {
                    Sol_cpy(sol,currentT[th_id]);
                    if(Sol_getCost(currentT[th_id]) < Sol_getCost(bestSol)) {
#ifdef DEBUG
                        th_id = omp_get_thread_num();
                        printf("\n%ld thread %d", Sol_getCost(bestSol), th_id);
#endif
                        Sol_cpy( bestSol, sol );
                    }
                    improve = 1;
                    if(firstImp && kN) {
                        abort = 1;
                        #pragma omp flush (abort)
                    }
                } else
                    Sol_cpy(currentT[th_id], sol);
                _time = omp_get_wtime()-tStart;
                if( _time > timeRem) {
                    abort = 1;
                    #pragma omp flush (abort)
                }
            }
        }
    }
    return improve;

}

/*Call neighborhood stochastic local search to VNS by idx assortment*/
int Neighbor_callStocChosen( Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *current, int kN,  Test *test)
{

    assert(neighborhood != NULL);
    assert(current != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);

    //int pos = kN;
    int pos = Neighbor_getIdxAssortment( neighborhood, kN-1);
    int valid = 0;



    assert( pos > 0 );
    assert( pos <= Neighbor_nNeighborhood(neighborhood) );

    switch(pos) {
        case seqInvert:
#ifdef DEBUG
            //  printf("\n[cSC] seqInvert : %d", seqInvert);
#endif
            valid = Neighbor_random_inv( neighborhood,  current );
            Neighbor_setLastNeigh(neighborhood,seqInvert);
            break;
        case seqShiftJob:
#ifdef DEBUG
            //printf("\n[cSC] seqShiftJob : %d", seqShiftJob);
#endif
            valid = Neighbor_random_shiftJob( neighborhood, current);
            Neighbor_setLastNeigh(neighborhood,seqShiftJob);
            break;
        case seqSwapJob:
#ifdef DEBUG
            //printf("\n[cSC] seqSwapJob : %d", seqSwapJob);
#endif
            valid = Neighbor_random_swapJob( neighborhood, current);
            Neighbor_setLastNeigh(neighborhood,seqSwapJob);
            break;
        case seqShiftProj:
#ifdef DEBUG
            //printf("\n[cSC] seqShiftProj : %d", seqShiftProj);
#endif
            valid = Neighbor_random_shiftProj( neighborhood, current);
            Neighbor_setLastNeigh(neighborhood,seqShiftProj);
            break;
        case seqSwapProj:
#ifdef DEBUG
            //printf("\n[cSC] seqSwapProj : %d", seqSwapProj);
#endif
            valid = Neighbor_random_swapTwoProj( neighborhood, current);
            Neighbor_setLastNeigh(neighborhood,seqSwapProj);
            break;
        case seqCompactProj:
#ifdef DEBUG
            // printf("\n[cSC] seqCompactProj : %d", seqCompactProj);
#endif
            valid = Neighbor_random_compactProj( neighborhood,  current);
            Neighbor_setLastNeigh(neighborhood,seqCompactProj);
            break;
        case changeOneMode:
#ifdef DEBUG
            // printf("\n[cSC] changeOneMode : %d", changeOneMode);
#endif
            valid = Neighbor_random_changeOneMode( neighborhood, sol, current);
            Neighbor_setLastNeigh(neighborhood,changeOneMode);
            break;
        case changeTwoMode:
#ifdef DEBUG
            //       printf("\n[cSC] changeTwoMode : %d", changeTwoMode);
#endif
            valid = Neighbor_random_changeTwoMode( neighborhood, sol, current);
            Neighbor_setLastNeigh(neighborhood,changeTwoMode);
            break;
        case changeThreeMode:
#ifdef DEBUG
            //     printf("\n[cSC] changeThreeMode : %d", changeThreeMode);
#endif
            valid = Neighbor_random_changeThreeMode( neighborhood, sol, current);
            Neighbor_setLastNeigh(neighborhood,changeThreeMode);
            break;
        case changeFourMode:
#ifdef DEBUG
            //   printf("\n[cSC] changeFourMode : %d", changeFourMode);
#endif
            valid = Neighbor_random_changeFourMode( neighborhood, sol, current);
            Neighbor_setLastNeigh(neighborhood,changeFourMode);
            break;
        case seqSwapJobFILS:
#ifdef DEBUG
            //printf("\n[cSC] seqSwapJobFILS : %d", seqSwapJobFILS);
#endif
            valid = Neighbor_random_swapJobFILS( neighborhood, sol, bestSol, current);
            Neighbor_setLastNeigh(neighborhood,seqSwapJobFILS);
            break;
        case seqInsertJobFILS:
#ifdef DEBUG
            //printf("\n[cSC] seqInsertJobFILS : %d", seqInsertJobFILS);
#endif
            valid = Neighbor_random_insertJobFILS( neighborhood, sol, bestSol, current);
            Neighbor_setLastNeigh(neighborhood,seqInsertJobFILS);
            break;
        case seqCompOnExtrem:
#ifdef DEBUG
            //  printf("\n[cSC] seqCompOnExtrem : %d", seqCompOnExtrem);
#endif
            valid = Neighbor_random_compactOnExtreme( neighborhood, current);
            Neighbor_setLastNeigh(neighborhood,seqCompOnExtrem);
            break;
        case seqMoveProj:
#ifdef DEBUG
            // printf("\n[cSC] seqMoveProj : %d", seqMoveProj);
#endif
            valid = Neighbor_random_moveProj( neighborhood, current);
            Neighbor_setLastNeigh(neighborhood,seqMoveProj);
            break;
    }
    return valid;
}

void Neighbor_Shake( VNS *vns, Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *current, int nMoves,  Test *test)
{
    assert(neighborhood != NULL);
    assert(current != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);




    int pos;
    int valid = 0;
    //100-200 aceita o melhor com o criterio de diversificacao.
    Cost oldFo;
    Solution *bestSolShake = Sol_create(Sol_inst(current));
    Solution *solInitial = Sol_create(Sol_inst(current));
    Sol_cpy(solInitial,sol);



    for(int i = 0 ; i < nMoves ; i++) {


        if(VNS_getDivRM(vns))
            VNS_increasingResidencyJobInMode(vns, current);
        if(VNS_getDivRJ(vns))
            VNS_increasingResidencyJobInSequence(vns, current);

        pos = Neighbor_roulette( neighborhood );

        assert( pos >= 0 );
        assert( pos < Neighbor_nNeighborhood(neighborhood) );

        Neighbor_setLastNeigh(neighborhood,pos+1);

        switch(pos) {
            case seqInvert-1:
#ifdef DEBUG
                // printf("\n[S] seqInvert : %d", seqInvert);
                // fflush(stdout);
                // fflush(stdin);
#endif
                valid = Neighbor_random_inv( neighborhood,  current );
                break;
            case seqShiftJob-1:
#ifdef DEBUG
                //printf("\n[S] seqShiftJob : %d", seqShiftJob);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_shiftJob( neighborhood, current);
                break;
            case seqSwapJob-1:
#ifdef DEBUG
                // printf("\n[S] seqSwapJob : %d", seqSwapJob);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_swapJob( neighborhood, current);
                break;
            case seqShiftProj-1:
#ifdef DEBUG
                //printf("\n[S] seqShiftProj : %d", seqShiftProj);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_shiftProj( neighborhood, current);
                break;
            case seqSwapProj-1:
#ifdef DEBUG
                //printf("\n[S] seqSwapProj : %d", seqSwapProj);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_swapTwoProj( neighborhood, current);
                break;
            case seqCompactProj-1:
#ifdef DEBUG
                //printf("\n[S] seqCompactProj : %d", seqCompactProj);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_compactProj( neighborhood,  current);
                break;
            case changeOneMode-1:
#ifdef DEBUG
                //printf("\n[S] changeOneMode : %d", changeOneMode);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_changeOneMode( neighborhood, solInitial, current);
                break;
            case changeTwoMode-1:
#ifdef DEBUG
                //printf("\n[S] changeTwoMode : %d", changeTwoMode);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_changeTwoMode( neighborhood, solInitial, current);
                break;
            case changeThreeMode-1:
#ifdef DEBUG
                //printf("\n[S] changeThreeMode : %d", changeThreeMode);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_changeThreeMode( neighborhood, solInitial, current);
                break;
            case changeFourMode-1:
#ifdef DEBUG
                //printf("\n[S] changeFourMode : %d", changeFourMode);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_changeFourMode( neighborhood, solInitial, current);
                break;
            case seqSwapJobFILS-1:
#ifdef DEBUG
                // printf("\n[S] seqSwapJobFILS : %d", seqSwapJobFILS);
                // fflush(stdout);
                // fflush(stdin);
#endif
                valid = Neighbor_random_swapJobFILS( neighborhood, solInitial, bestSol, current);
                break;
            case seqInsertJobFILS-1:
#ifdef DEBUG
                //printf("\n[S] seqInsertJobFILS : %d", seqInsertJobFILS);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_insertJobFILS( neighborhood, solInitial, bestSol, current);
                break;
            case seqCompOnExtrem-1:
#ifdef DEBUG
                //printf("\n[S] seqCompOnExtrem : %d", seqCompOnExtrem);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_compactOnExtreme( neighborhood, current);
                break;
            case seqMoveProj-1:
#ifdef DEBUG
                //printf("\n[S] seqMoveProj : %d", seqMoveProj);
                //fflush(stdout);
                //fflush(stdin);

#endif
                valid = Neighbor_random_moveProj( neighborhood, current);
                break;
        }
        if(valid) {
            Sol_rebuild_opt(current, solInitial);
            Sol_cpy(solInitial,current);
        }

        if(i==0) {
            oldFo = Sol_getCost(current);
            Sol_cpy(bestSolShake,current);
        }

        Cost fo = Sol_getCost(current);


        /*penalizing residency*/

        //   if( nItDiversification == lahc->nDiversification || (nItStay < lahc->nStayDiversification && nItStay >= 0)) {

        int applyDiversitication = 15;//NT_RANDOM( 14 );
        //printf("applyDiversitication %d pos %d\n", applyDiversitication,pos);

        //if(pos ==13)  getchar();
        if(applyDiversitication == pos) {

            if(VNS_getDivRM(vns)) {
#ifdef DEBUG
                //          printf("\nPenalizing residency of modes...\n");
                fflush(stdout);
                fflush(stdin);
#endif
                fo += VNS_penaltyModes(vns, neighborhood, Sol_getCost(current) );
            }
            if(VNS_getDivRJ(vns)) {
#ifdef DEBUG
                //            printf("\nPenalizing residency of jobs...\n");
                fflush(stdout);
                fflush(stdin);
#endif
                fo += VNS_penaltyJobs(vns, neighborhood, Sol_getCost(current) );
            }
            if(VNS_getDivTM(vns)) {
#ifdef DEBUG
                //        printf("\nPenalizing transitivity of modes...\n");
                fflush(stdout);
                fflush(stdin);
#endif
                fo += VNS_penaltyTransModes(vns, neighborhood, Sol_getCost(current) );
            }
            if(VNS_getDivTJ(vns)) {
                fo += VNS_penaltyTransJobs(vns, neighborhood, Sol_getCost(current) );
#ifdef DEBUG
                //      printf("\nPenalizing transitivity of jobs...\n");
                fflush(stdout);
                fflush(stdin);
#endif

            }
        }


        if( fo < oldFo) {

            Sol_cpy( bestSolShake, current );
            oldFo = fo;
#ifdef DEBUG
            //    printf("\nShake %ld I %d ", Sol_getCost(current), i );
#endif
            /*increasing the transitivity*/
            if(VNS_getDivTM(vns))
                VNS_increasingTransitivityOfModes(vns, neighborhood);
            if(VNS_getDivTJ(vns))
                VNS_increasingTransitivityOfSequence(vns, neighborhood);
            if(Sol_getCost(current)<Sol_getCost(bestSol))
                Sol_cpy( bestSol, current );

        }

        Neighbor_setNullLastJobModify(neighborhood);


    }

    Sol_cpy( current, bestSolShake );

    Sol_free(&solInitial);
    Sol_free(&bestSolShake);


}


/*Call neighborhood stochastic local search to VNS randomly*/
int Neighbor_callStocRandom( Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *current)
{

    assert(neighborhood != NULL);
    assert(current != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);

    int valid = 0;

    int pos = Neighbor_roulette( neighborhood );

    assert( pos >= 0 );
    assert( pos < Neighbor_nNeighborhood(neighborhood) );

    Neighbor_setLastNeigh(neighborhood,pos+1);

    switch(pos) {
        case seqInvert-1:
#ifdef DEBUG
            //  printf("\n[cSC] seqInvert : %d", seqInvert);
#endif
            valid = Neighbor_random_inv( neighborhood,  current );
            break;
        case seqShiftJob-1:
#ifdef DEBUG
            //printf("\n[cSC] seqShiftJob : %d", seqShiftJob);
#endif
            valid = Neighbor_random_shiftJob( neighborhood, current);
            break;
        case seqSwapJob-1:
#ifdef DEBUG
            //printf("\n[cSC] seqSwapJob : %d", seqSwapJob);
#endif
            valid = Neighbor_random_swapJob( neighborhood, current);
            break;
        case seqShiftProj-1:
#ifdef DEBUG
            //printf("\n[cSC] seqShiftProj : %d", seqShiftProj);
#endif
            valid = Neighbor_random_shiftProj( neighborhood, current);
            break;
        case seqSwapProj-1:
#ifdef DEBUG
            //printf("\n[cSC] seqSwapProj : %d", seqSwapProj);
#endif
            valid = Neighbor_random_swapTwoProj( neighborhood, current);
            break;
        case seqCompactProj-1:
#ifdef DEBUG
            // printf("\n[cSC] seqCompactProj : %d", seqCompactProj);
#endif
            valid = Neighbor_random_compactProj( neighborhood,  current);
            break;
        case changeOneMode-1:
#ifdef DEBUG
            // printf("\n[cSC] changeOneMode : %d", changeOneMode);
#endif
            valid = Neighbor_random_changeOneMode( neighborhood, sol, current);
            break;
        case changeTwoMode-1:
#ifdef DEBUG
            //       printf("\n[cSC] changeTwoMode : %d", changeTwoMode);
#endif
            valid = Neighbor_random_changeTwoMode( neighborhood, sol, current);
            break;
        case changeThreeMode-1:
#ifdef DEBUG
            //     printf("\n[cSC] changeThreeMode : %d", changeThreeMode);
#endif
            valid = Neighbor_random_changeThreeMode( neighborhood, sol, current);
            break;
        case changeFourMode-1:
#ifdef DEBUG
            //   printf("\n[cSC] changeFourMode : %d", changeFourMode);
#endif
            valid = Neighbor_random_changeFourMode( neighborhood, sol, current);
            break;
        case seqSwapJobFILS-1:
#ifdef DEBUG
            //printf("\n[cSC] seqSwapJobFILS : %d", seqSwapJobFILS);
#endif
            valid = Neighbor_random_swapJobFILS( neighborhood, sol, bestSol, current);
            break;
        case seqInsertJobFILS-1:
#ifdef DEBUG
            //printf("\n[cSC] seqInsertJobFILS : %d", seqInsertJobFILS);
#endif
            valid = Neighbor_random_insertJobFILS( neighborhood, sol, bestSol, current);
            break;
        case seqCompOnExtrem-1:
#ifdef DEBUG
            //  printf("\n[cSC] seqCompOnExtrem : %d", seqCompOnExtrem);
#endif
            valid = Neighbor_random_compactOnExtreme( neighborhood, current);
            break;
        case seqMoveProj-1:
#ifdef DEBUG
            // printf("\n[cSC] seqMoveProj : %d", seqMoveProj);
#endif
            valid = Neighbor_random_moveProj( neighborhood, current);
            break;
    }
    return valid;
}


void Neighbor_SmartShake( VNS *vns, Neighborhood *neighborhood, Solution *sol, Solution *bestSol, int nMovesShake, int sizeSamplingShake, double timeRem)
{
    assert(neighborhood != NULL);
    assert(vns != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);

    int valid = 0;
    double _time = 0;
    clock_t tStart = clock();

    Solution *current = Sol_create(Sol_inst(sol));


    for(int m = 1 ; m <= nMovesShake ; m++) {

        Solution *bestSolShake = Sol_create(Sol_inst(sol));
        Sol_cpy(bestSolShake,sol);
        Cost bestCost = INT_MAX;

        for(int i = 1 ; i <= sizeSamplingShake ; i++) {

            Sol_cpy(current,sol);

            valid = Neighbor_callStocRandom(neighborhood, sol, bestSol, current);
            if(valid) {
                Sol_rebuild_opt(current, sol);

                Cost fo = Sol_getCost(current);

                //                printf("\n FO %ld\n", fo );
                /*penalizing residency and transitivity*/
                if(VNS_getDivRM(vns)) {
#ifdef DEBUG
                    printf("\nPenalizing residency of modes...\n");
#endif
                    fo += VNS_penaltyModes(vns, neighborhood, Sol_getCost(current) );
                }
                if(VNS_getDivRJ(vns)) {
#ifdef DEBUG
                    printf("\nPenalizing residency of jobs...\n");
#endif
                    fo += VNS_penaltyJobs(vns, neighborhood, Sol_getCost(current) );
                }
                if(VNS_getDivTM(vns)) {
#ifdef DEBUG
                    printf("\nPenalizing transitivity of modes...\n");
#endif
                    fo += VNS_penaltyTransModes(vns, neighborhood, Sol_getCost(current) );
                }
                if(VNS_getDivTJ(vns)) {
#ifdef DEBUG
                    printf("\nPenalizing transitivity of jobs...\n");
#endif
                    fo += VNS_penaltyTransJobs(vns, neighborhood, Sol_getCost(current) );
                }


                //                printf("\n FO_Penalty %ld < bestCost %ld \n", fo,  bestCost);
                if( fo < bestCost) {

                    Sol_cpy( bestSolShake, current );
                    bestCost = Sol_getCost(bestSolShake);
#ifdef DEBUG
                    printf("\n BestShakeCurrent %ld SolInit %ld \n nMoveShake %d sizeSamplingShake %d ", Sol_getCost(current), Sol_getCost(sol), m, i );
                    //           getchar();
#endif

                    /* garantindo que se for melhor atualiza o best do VNS*/
                    //  if(Sol_getCost(current)<Sol_getCost(bestSol))
                    //    Sol_cpy( bestSol, current );

                }
            }

            /*set null the last change made for the last neighborhood*/
            Neighbor_setNullLastJobModify(neighborhood);

            _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            if(_time >= timeRem ) {
                Sol_cpy(sol,bestSolShake);
                Sol_free(&bestSolShake);
                Sol_free(&current);
                return 0;
            }
        }
        Sol_cpy(sol,bestSolShake);
        Sol_free(&bestSolShake);

    }

    Sol_free(&current);

}



void Neighbor_Shake2( VNS *vns, Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *current, int kN, int nMoves,  Test *test)
{
    assert(neighborhood != NULL);
    assert(current != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);

    int pos = Neighbor_getIdxAssortment( neighborhood, kN-1);

    assert( pos > 0 );
    assert( pos <= Neighbor_nNeighborhood(neighborhood) );

    int valid = 0;

    Cost oldFo;
    Solution *bestSolShake = Sol_create(Sol_inst(current));
    Solution *solInitial = Sol_create(Sol_inst(current));
    Sol_cpy(solInitial,sol);

    for(int i = 0 ; i < nMoves ; i++) {


        if(VNS_getDivRM(vns))
            VNS_increasingResidencyJobInMode(vns, current);
        if(VNS_getDivRJ(vns))
            VNS_increasingResidencyJobInSequence(vns, current);


        Neighbor_setLastNeigh(neighborhood,pos);

        switch(pos) {
            case seqInvert:
#ifdef DEBUG
                // printf("\n[S] seqInvert : %d", seqInvert);
                // fflush(stdout);
                // fflush(stdin);
#endif
                valid = Neighbor_random_inv( neighborhood,  current );
                break;
            case seqShiftJob:
#ifdef DEBUG
                //printf("\n[S] seqShiftJob : %d", seqShiftJob);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_shiftJob( neighborhood, current);
                break;
            case seqSwapJob:
#ifdef DEBUG
                // printf("\n[S] seqSwapJob : %d", seqSwapJob);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_swapJob( neighborhood, current);
                break;
            case seqShiftProj:
#ifdef DEBUG
                //printf("\n[S] seqShiftProj : %d", seqShiftProj);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_shiftProj( neighborhood, current);
                break;
            case seqSwapProj:
#ifdef DEBUG
                //printf("\n[S] seqSwapProj : %d", seqSwapProj);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_swapTwoProj( neighborhood, current);
                break;
            case seqCompactProj:
#ifdef DEBUG
                //printf("\n[S] seqCompactProj : %d", seqCompactProj);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_compactProj( neighborhood,  current);
                break;
            case changeOneMode:
#ifdef DEBUG
                //printf("\n[S] changeOneMode : %d", changeOneMode);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_changeOneMode( neighborhood, solInitial, current);
                break;
            case changeTwoMode:
#ifdef DEBUG
                //printf("\n[S] changeTwoMode : %d", changeTwoMode);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_changeTwoMode( neighborhood, solInitial, current);
                break;
            case changeThreeMode:
#ifdef DEBUG
                //printf("\n[S] changeThreeMode : %d", changeThreeMode);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_changeThreeMode( neighborhood, solInitial, current);
                break;
            case changeFourMode:
#ifdef DEBUG
                //printf("\n[S] changeFourMode : %d", changeFourMode);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_changeFourMode( neighborhood, solInitial, current);
                break;
            case seqSwapJobFILS:
#ifdef DEBUG
                // printf("\n[S] seqSwapJobFILS : %d", seqSwapJobFILS);
                // fflush(stdout);
                // fflush(stdin);
#endif
                valid = Neighbor_random_swapJobFILS( neighborhood, solInitial, bestSol, current);
                break;
            case seqInsertJobFILS:
#ifdef DEBUG
                //printf("\n[S] seqInsertJobFILS : %d", seqInsertJobFILS);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_insertJobFILS( neighborhood, solInitial, bestSol, current);
                break;
            case seqCompOnExtrem:
#ifdef DEBUG
                //printf("\n[S] seqCompOnExtrem : %d", seqCompOnExtrem);
                //fflush(stdout);
                //fflush(stdin);
#endif
                valid = Neighbor_random_compactOnExtreme( neighborhood, current);
                break;
            case seqMoveProj:
#ifdef DEBUG
                //printf("\n[S] seqMoveProj : %d", seqMoveProj);
                //fflush(stdout);
                //fflush(stdin);

#endif
                valid = Neighbor_random_moveProj( neighborhood, current);
                break;
        }
        if(valid) {
            Sol_rebuild_opt(current, solInitial);
            Sol_cpy(solInitial,current);
        }

        if(i==0) {
            oldFo = Sol_getCost(current);
            Sol_cpy(bestSolShake,current);
        }

        Cost fo = Sol_getCost(current);


        /*penalizing residency*/


        //int applyDiversitication = INT_RANDOM( 14 );

        //  if(applyDiversitication == pos){

        if(VNS_getDivRM(vns)) {
#ifdef DEBUG
            //          printf("\nPenalizing residency of modes...\n");
            fflush(stdout);
            fflush(stdin);
#endif
            fo += VNS_penaltyModes(vns, neighborhood, Sol_getCost(current) );
        }
        if(VNS_getDivRJ(vns)) {
#ifdef DEBUG
            //            printf("\nPenalizing residency of jobs...\n");
            fflush(stdout);
            fflush(stdin);
#endif
            fo += VNS_penaltyJobs(vns, neighborhood, Sol_getCost(current) );
        }
        if(VNS_getDivTM(vns)) {
#ifdef DEBUG
            //        printf("\nPenalizing transitivity of modes...\n");
            fflush(stdout);
            fflush(stdin);
#endif
            fo += VNS_penaltyTransModes(vns, neighborhood, Sol_getCost(current) );
        }
        if(VNS_getDivTJ(vns)) {
            fo += VNS_penaltyTransJobs(vns, neighborhood, Sol_getCost(current) );
#ifdef DEBUG
            //      printf("\nPenalizing transitivity of jobs...\n");
            fflush(stdout);
            fflush(stdin);
#endif

        }
        //   }

        if( fo < oldFo) {

            Sol_cpy( bestSolShake, current );
            oldFo = fo;
#ifdef DEBUG
            printf("\nShake same neigh %ld I:nMoves %d ", Sol_getCost(current), i );
#endif
            /*increasing the transitivity*/
            if(VNS_getDivTM(vns))
                VNS_increasingTransitivityOfModes(vns, neighborhood);
            if(VNS_getDivTJ(vns))
                VNS_increasingTransitivityOfSequence(vns, neighborhood);
            if(Sol_getCost(current)<Sol_getCost(bestSol))
                Sol_cpy( bestSol, current );

        }

        Neighbor_setNullLastJobModify(neighborhood);


    }

    Sol_cpy( current, bestSolShake );

    Sol_free(&solInitial);
    Sol_free(&bestSolShake);


}


void Neighbor_setLastNeigh(Neighborhood *neigh, int value)
{
    assert( neigh != NULL );
    assert( value >= 0 );

    neigh->lastN = value;

}

int Neighbor_getLastNeigh(Neighborhood *neigh)
{
    assert( neigh != NULL );

    return neigh->lastN;

}

int Neighbor_getLastJob(Neighborhood *neigh, int idx)
{
    assert( neigh != NULL );

    return neigh->lastJ[idx];

}

int Neighbor_getNewModes(Neighborhood *neigh, int idx)
{
    assert( neigh != NULL );

    return neigh->newModes[idx];

}

int Neighbor_getLastJobModify(Neighborhood *neigh, int idx)
{

    assert( neigh != NULL );

    return neigh->nLastJModify[idx];

}


int Neighbor_getPosLastJobModify(Neighborhood *neigh, int idx)
{

    assert( neigh != NULL );

    return neigh->posNLastJModify[idx];

}

void Neighbor_setNullLastJobModify(Neighborhood *neigh)
{

    assert( neigh != NULL );


    CLEAR_VECTOR( neigh->lastJ, int, 4);
    CLEAR_VECTOR( neigh->newModes, int, 4);
    CLEAR_VECTOR( neigh->nLastJModify, int, neigh->contLastJ);
    CLEAR_VECTOR( neigh->posNLastJModify, int, neigh->contLastJ);

    neigh->contLastJ = 0;

}

int Neighbor_getContLastJ(Neighborhood *neigh)
{

    assert( neigh != NULL );

    return neigh->contLastJ;

}

/*Call neighborhood stochastic local search to LAHC and SA*/
int Neighbor_callStocLS( Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *current, LearningAutomata *la, int learning)
{

    assert(neighborhood != NULL);
    assert(current != NULL);
    assert(sol != NULL);
    assert(bestSol != NULL);

    int pos;
    if(learning) {
        pos = LA_next( la );
        //  printf("learnig %d", pos);
    } else {
        pos = Neighbor_roulette( neighborhood );
        //  printf("roulette %d", pos);
    }

    int valid = 0;
    assert( pos >= 0 );
    assert( pos < Neighbor_nNeighborhood(neighborhood) );


    switch(pos) {
        case seqInvert-1:
            //printf("\n[S] seqInvert : %d", seqInvert);
            valid = Neighbor_random_inv( neighborhood,  current );
            Neighbor_setLastNeigh(neighborhood,seqInvert-1);
            break;
        case seqShiftJob-1:
            //   printf("\n[S] seqShiftJob : %d", seqShiftJob);
            valid = Neighbor_random_shiftJob( neighborhood, current);
            Neighbor_setLastNeigh(neighborhood,seqShiftJob-1);
            break;
        case seqSwapJob-1:
            //   printf("\n[S] seqSwapJob : %d", seqSwapJob);
            valid = Neighbor_random_swapJob( neighborhood, current);
            Neighbor_setLastNeigh(neighborhood,seqSwapJob-1);
            break;
        case seqShiftProj-1:
            if(neighborhood->typeinstance==2) {
                //   printf("\n[S] seqShiftProj : %d", seqShiftProj);
                valid = Neighbor_random_shiftProj( neighborhood,  current);
                Neighbor_setLastNeigh(neighborhood,seqShiftProj-1);
                break;
            }
        case seqSwapProj-1:
            if(neighborhood->typeinstance==2) {
                //    printf("\n[S] seqSwapProj : %d", seqSwapProj);
                valid =  Neighbor_random_swapTwoProj( neighborhood,  current);
                Neighbor_setLastNeigh(neighborhood,seqSwapProj-1);
                break;
            }
        case seqCompactProj-1:
            if(neighborhood->typeinstance==2) {
                //   printf("\n[S] seqCompactProj : %d", seqCompactProj);
                valid = Neighbor_random_compactProj( neighborhood,  current);
                Neighbor_setLastNeigh(neighborhood,seqCompactProj-1);
                break;
            }
        case changeOneMode-1:
            if(neighborhood->typeinstance>=1) {
                //printf("\n[S] changeOneMode : %d", changeOneMode);
                valid = Neighbor_random_changeOneMode( neighborhood, sol, current);
                Neighbor_setLastNeigh(neighborhood,changeOneMode-1);
                break;
            }
        case changeTwoMode-1:
            if(neighborhood->typeinstance>=1) {
                //   printf("\n[S] changeTwoMode : %d", changeTwoMode);
                valid = Neighbor_random_changeTwoMode( neighborhood, sol, current);
                Neighbor_setLastNeigh(neighborhood,changeTwoMode-1);
                break;
            }
        case changeThreeMode-1:
            //printf("\n[S] changeThreeMode : %d", changeThreeMode);
            valid = Neighbor_random_changeThreeMode( neighborhood, sol, current);
            Neighbor_setLastNeigh(neighborhood,changeThreeMode-1);
            break;
        case changeFourMode-1:
            if(neighborhood->typeinstance>=1) {
                //printf("\n[S] changeFourMode : %d", changeFourMode);
                valid = Neighbor_random_changeFourMode( neighborhood,sol, current);
                Neighbor_setLastNeigh(neighborhood,changeFourMode-1);
                break;
            }
        case seqSwapJobFILS-1:
            //    printf("\n[S] seqSwapJobFILS : %d", seqSwapJobFILS);
            valid = Neighbor_random_swapJobFILS( neighborhood, sol, bestSol, current);
            Neighbor_setLastNeigh(neighborhood,seqSwapJobFILS-1);
            break;
        case seqInsertJobFILS-1:
            //   printf("\n[S] seqInsertJobFILS : %d", seqInsertJobFILS);
            valid = Neighbor_random_insertJobFILS( neighborhood, sol, bestSol, current);
            Neighbor_setLastNeigh(neighborhood,seqInsertJobFILS-1);
            break;
        case seqCompOnExtrem-1:
            //   printf("\n[S] seqCompOnExtrem : %d", seqCompOnExtrem);
            valid = Neighbor_random_compactOnExtreme( neighborhood, current);
            Neighbor_setLastNeigh(neighborhood,seqCompOnExtrem-1);
            break;
        case seqMoveProj-1:
            if(neighborhood->typeinstance==2) {
                //    printf("\n[S] seqMoveProj : %d", seqMoveProj);
                valid = Neighbor_random_moveProj( neighborhood, current);
                Neighbor_setLastNeigh(neighborhood,seqMoveProj-1);
                break;
            }
    }

    return valid;

}

long double *Neighbor_getIntensity(Neighborhood *neighborhood)
{
    assert(neighborhood != NULL);
    return neighborhood->intensity;
}








/*mode set*/
//#define NEIGHBOR_CHANGEONEMODE_STSMS 0
//#define NEIGHBOR_CHANGETWOMODE_STSMS 1
// #define N_NEIGHBORHOOD_MS 2

/*
Neighborhood *Neighbor_createMS( const Instance* inst )
{

    assert( inst!=NULL );

    Neighborhood *neighborhood;
    ALLOCATE_INI(neighborhood, Neighborhood);

    neighborhood->nNeighborhood = N_NEIGHBORHOOD_MS;

    double *intensity;
    ALLOCATE_VECTOR_INI(intensity, double, N_NEIGHBORHOOD_MS);

    intensity[0] = 1.0;
    intensity[1] = 1.0;

    neighborhood->intensity = intensity;
    neighborhood->inst = inst;
    neighborhood->penaltyChangeMode = 0;

    return neighborhood;
}
*/





/*void Neighbor_setIntensity(Neighborhood *neighborhood, int idxNeighbor, double intens)
{
    assert(neighborhood != NULL);
    assert(idxNeighbor >= 0 && idxNeighbor < neighborhood->nNeighborhood);

    neighborhood->intensity[idxNeighbor] = intens;
}


*/
/*
void Neighbor_setIdxAssortment(Neighborhood *neighborhood, int idx, int idxNeighbor)
{
    assert(neighborhood != NULL);
    assert(idx >= 0 && idx < neighborhood->nNeighborhood);
    assert(idxNeighbor > 0 && idxNeighbor <= neighborhood->nNeighborhood);

    neighborhood->assortment[idx] = idxNeighbor;
}
*/

/*
int Neighbor_getIt(Neighborhood *neig)
{
    return neig->it;
}
*/


/*neighborhood stochastic mode set min*/
/*
void Neighbor_callStocMS( Neighborhood *neighborhood, ModeSet *modeSet)
{

    int pos = Neighbor_roulette(neighborhood);

    assert( pos>=0 );
    assert( pos<Neighbor_nNeighborhood(neighborhood) );

    switch(pos) {
        case NEIGHBOR_CHANGEONEMODE_STSMS:
            Neighbor_changeOneModeStocMS( neighborhood, modeSet );
            break;
        case NEIGHBOR_CHANGETWOMODE_STSMS:
            Neighbor_changeTwoModeStocMS( neighborhood, modeSet );
            break;
    }
}

void Neighbor_changeOneModeStocMS(Neighborhood *neighborhood, ModeSet *modeSet)
{
    assert(modeSet != NULL);

    int j = INT_RANDOM_LU(Modes_firstJob(modeSet), Modes_lastJob(modeSet));

    const Job *job = Inst_job( Modes_inst(modeSet), j);
    int m, cont=0;
    const Mode* mode;
    if(Job_nModes(job) > 1) {
        do {
            m = INT_RANDOM(Job_nModes(job));
            mode = Job_mode(job,m);
            cont++;
        } while( ( Mode_index(mode) == Modes_job(modeSet,j) ) || ( !Mode_isFeasible( Modes_inst(modeSet), mode) && cont <3) );
        if(cont< 3) Modes_modifyCount( modeSet, j, m);
    }
}

void Neighbor_changeTwoModeStocMS(Neighborhood *neighborhood, ModeSet *modeSet)
{

    assert(modeSet != NULL);

    int j1 = INT_RANDOM_LU(Modes_firstJob(modeSet), Modes_lastJob(modeSet));
    int j2;
    do {
        j2 = INT_RANDOM_LU(Modes_firstJob(modeSet), Modes_lastJob(modeSet));
    } while(j1==j2);

    const Job *job1 = Inst_job( Modes_inst(modeSet), j1);
    const Job *job2 = Inst_job( Modes_inst(modeSet), j2);
    int m1,m2, cont = 0;
    const Mode* mode1;
    if(Job_nModes(job1) > 1) {
        do {
            m1 = INT_RANDOM(Job_nModes(job1));
            mode1 = Job_mode(job1,m1);
            cont++;
        } while( ( Mode_index(mode1) == Modes_job(modeSet,j1) ) || ( !Mode_isFeasible( Modes_inst(modeSet), mode1)&& cont <3 ) );
        if(cont< 3)Modes_modifyCount( modeSet, j1, m1);
    }

    cont = 0;
    const Mode* mode2;
    if(Job_nModes(job2) > 1) {
        do {
            m2 = INT_RANDOM(Job_nModes(job2));
            mode2 = Job_mode(job2,m2);
            cont++;
        } while( ( Mode_index(mode2) == Modes_job(modeSet,j2) ) || ( !Mode_isFeasible( Modes_inst(modeSet), mode2 )&& cont <3) );
        if(cont< 3) Modes_modifyCount( modeSet, j2,m2);
    }
}

*/
