#ifndef NEIGHBORHOOD_H
#define NEIGHBORHOOD_H

#include "solution.h"
#include "learning.h"
#include "test.h"

typedef struct _Neighborhood Neighborhood;

enum NeighborhoodType { seqInvert=1, seqShiftJob=2, seqSwapJob=3, seqShiftProj=4, seqSwapProj=5, seqCompactProj=6,
                        changeOneMode=7, changeTwoMode=8, changeThreeMode=9, changeFourMode=10,
                        seqSwapJobFILS=11, seqInsertJobFILS=12, seqCompOnExtrem=13, seqMoveProj=14
                      };

/* returns an initialized object Neighborhood by default to Solution*/
Neighborhood *Neighbor_create( const Instance* inst, char **argv, int argc );
void Neighbor_printProbabilitiesTime( Neighborhood *neighborhood, int sw, int it, float psw, char *argv, double timePrint);
void Neighbor_printProbabilities( Neighborhood *neighborhood, int sw, int it, float psw, char *argv);
int Neighbor_getTimePrint( Neighborhood *neighborhood);

double Neighbor_getIdP(Neighborhood *neighborhood, int idxNeighbor, int sw);

int Neighbor_getLastNeigh(Neighborhood *neigh);
int Neighbor_getLastJob(Neighborhood *neigh, int idx);
int Neighbor_getPosLastJobModify(Neighborhood *neigh, int idx);
int Neighbor_getNewModes(Neighborhood *neigh, int idx);
/* updates the parameters through the arguments .
  -intensity idx value
  -assort idx value
  -minK idx value
  -maxK idx value
*/
void Neighbor_checkArgs(Neighborhood *neighbor, char **argv, int argc);

/* returns the number of neighborhood */
int Neighbor_nNeighborhood( Neighborhood *neighborhood );

/* roulette method */
int Neighbor_roulette( Neighborhood *neighborhood );

/* updates neighborhood intensities */
void Neighbor_getUpdatesIntensity(Neighborhood *neighborhood);

/* returns the vector intensity or the intensity of each neighborhood by idxNeighbor */
long double *Neighbor_getIntensity( Neighborhood *neighborhood );
long double Neighbor_getIdIntensity(Neighborhood *neighborhood, int idxNeighbor);

/* returns the vector assortment or the neighbor at idx position in vector assort  */
int *Neighbor_getAssortment(Neighborhood *neighborhood);
int Neighbor_getIdxAssortment(Neighborhood *neighborhood, int idx);

/* returns the min or max K of each neighborhood*/
int Neighbor_getMinK(Neighborhood *neighborhood, int idx);
int Neighbor_getMaxK(Neighborhood *neighborhood, int idx);

/* set the numbers of iterations on neighborhood*/
void Neighbor_setIt(Neighborhood *neig);
void Neighbor_setTime( Neighborhood *neighborhood);
void Neighbor_incrementI(Neighborhood *neighbor, int i);
void Neighbor_incrementEQ(Neighborhood *neighbor, int i);
void Neighbor_setF(Neighborhood *neighborhood);
void Neighbor_normF(Neighborhood *neighborhood);
void Neighbor_updatesIntensities(Neighborhood *neighborhood, int sw);
void Neighbor_clearImpTime(Neighborhood *neighborhood);
void Neighbor_setTE(Neighborhood *neighbor, int i, long double value);
void Neighbor_setTIV(Neighborhood *neighbor, int i, long double value);
void Neighbor_setTI(Neighborhood *neighbor, int i, long double value);

void Neighbor_setLastNeigh(Neighborhood *neigh, int value);


/*Verifies if the move infringes the rules of precedence.*/
int Neighbor_verifyPredInv( Solution *sol, int p1, int p2 );
int Neighbor_verifyPredShiftJob( Solution *sol, int p1, int p2 );
int Neighbor_verifyPredSwapJob( Solution *sol, int p1, int p2 );
int Neighbor_verifyPredSwapJobFILS( Solution *sol, int p1, int p2 );
int Neighbor_verifyPredInsertJobFILS( Solution *sol, int p1, int p2 );
void Neighbor_checkMoveSW(Neighborhood *neighborhood, Solution * sol, Solution * oldSol);
int Neighbor_getMoveSW(Neighborhood *neighborhood);

long double Neighbor_getI(Neighborhood *neighbor, int i);
int Neighbor_getLastNeigh(Neighborhood *neigh);
long double Neighbor_getEQ(Neighborhood *neighbor, int i);
long double Neighbor_getTE(Neighborhood *neighbor, int i);
long double Neighbor_getTI(Neighborhood *neighbor, int i);
long double Neighbor_getTIV(Neighborhood *neighbor, int i);

void Neighbor_setTI(Neighborhood *neighbor, int i, long double value);
int Neighbor_getLastJobModify(Neighborhood *neigh, int idx);
void Neighbor_setNullLastJobModify(Neighborhood *neigh);
int Neighbor_getContLastJ(Neighborhood *neigh);
int Neighbor_getLastJobModify(Neighborhood *neigh, int idx);
/* runs vnd whit neighborhood deterministic: first or best improvement */
int Neighbor_callDetLS( Neighborhood *neighborhood,
                        Solution *sol,
                        Solution *bestSol,
                        Solution *current,
                        //   int *nChanges,
                        //   int **nChangesJobsModes,
                        int kN,
                        int firstImp,
                        double timeRem,
                        Test *test);

/* runs the neighborhood method for Solution whit neighborhood stochastic*/
int Neighbor_callStocLS( Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *current,  LearningAutomata *la, int learning);
int Neighbor_callStocChosen( Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *current, int kN, Test *test);
//int Neighbor_Shake( VNS *vns, Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *current, int nMoves, Test *test);

int Neighbor_random_changeNModes(Neighborhood *neighborhood, Solution *sol,  Solution *current, int nModes, int **residency);

/*Change one mode on ModeSet */
int Neighbor_random_changeOneMode(Neighborhood *neighborhood,Solution *sol,  Solution *current);
/*Changes a pairs of modes on ModeSet */
int Neighbor_random_changeTwoMode(Neighborhood *neighborhood,Solution *sol, Solution *current);
/*Changes a three of modes on ModeSet */
int Neighbor_random_changeThreeMode(Neighborhood *neighborhood,Solution *sol, Solution *current);
/*Changes a four of modes on ModeSet */
int Neighbor_random_changeFourMode(Neighborhood *neighborhood,Solution *sol,  Solution *current);

/*After randomly selecting a position i of the sequence and a window size k, all values in positions
{i, . . . , p + k − 1} are removed and reinserted in inverse order on sequence if it is valid*/
int Neighbor_random_inv(Neighborhood *neighborhood,Solution *current);

/*After randomly selecting a job, one direction (left/right) and a number of positions k,
a job is moved k positions in the selected direction if it is valid.*/
int Neighbor_random_shiftJob(Neighborhood *neighborhood, Solution *current);

/* Two random positions of the sequence separated by at most k positions are randomly exchanged if it is valid*/
int Neighbor_random_swapJob(Neighborhood *neighborhood, Solution *current);

/*select a job j in a sequence received as a parameter.
A window W of random length l is set at random between the positions of
the nearest ancestor j (posd) and the closest successor of j (poss), so that start position
of W > posd and the end position of W < poss. Next, j is changed sequentially with
each job of W, until in an improvement in the objective function happens, but just in case it is valid.
*/
int Neighbor_random_swapJobFILS(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *current );

/* works with window predecessors and successors of a job j chosen randomly,
j is included at each position, shifting the others forward, until an
improved solution is obtained, but just in case it is valid.*/
int Neighbor_random_insertJobFILS(Neighborhood *neighborhood,Solution *sol, Solution *bestSol, Solution *current  );


/*Moves forward or backward all jobs of a project p on k positions */
int Neighbor_random_shiftProj(Neighborhood *neighborhood, Solution *current);

/*This movement receives as a parameter perc in (0, 1] determines percentage of jobs which
will be compressed. The project to be compressed is selected randomly. */
int Neighbor_random_compactProj(Neighborhood *neighborhood, Solution *current);

/* swaps two projects received as a parameter */
int Neighbor_random_swapTwoProj(Neighborhood *neighborhood,Solution *current);

/*This neighborhood  compresses a given project starting at some position i.
From this position on, all jobs will be adjacent on the sequence. */
int Neighbor_random_compactOnExtreme(Neighborhood *neighborhood, Solution *current);

/* after randomly selecting a position of the sequence, the project related to
the job in the selected position and k subsequent projects which appear in sequence
are considered. Jobs from all these projects are compacted in the selected position,
starting from the left or from the right. */
int Neighbor_random_moveProj(Neighborhood *neighborhood, Solution *current);

/*Neighborhood to call deterministic*/
int Neighbor_changeFourMode(Neighborhood *neighborhood,Solution *sol,int j1, int j2, int j3, int j4, int m1, int m2, int m3, int m4);
int Neighbor_changeThreeMode(Neighborhood *neighborhood,Solution *sol, int j1, int j2, int j3, int m1, int m2, int m3);
int Neighbor_changeTwoMode(Neighborhood *neighborhood,Solution *sol, int j1, int j2, int m1, int m2);
int Neighbor_changeOneMode(Neighborhood *neighborhood,Solution *sol, int j, int m);

int Neighbor_swapJob(Neighborhood *neighborhood,Solution *sol, int j1, int j2);
int Neighbor_shiftJob(Neighborhood *neighborhood,Solution *sol,  int i, int k, int dir);
int Neighbor_inv(Solution *sol, int i, int k);
int Neighbor_swapTwoProj(Neighborhood *neighborhood,Solution *sol, double p1, int p2);
int Neighbor_shiftProj(Neighborhood *neighborhood,Solution *sol, int k, int p, int dir);
int Neighbor_compactProj(Neighborhood *neighborhood,Solution *sol, double perc, int p);
int Neighbor_swapJobFILS(Neighborhood *neighborhood, Solution *current, Solution *sol, Solution *bestSol, int j1 );
int Neighbor_insertJobFILS(Neighborhood *neighborhood, Solution *current, Solution *sol, Solution *bestSol, int j1 );
int Neighbor_compactOnExtreme(Neighborhood *neighborhood,Solution *sol, int j1);
int Neighbor_moveProj(Neighborhood *neighborhood,Solution *sol, int nP, int dir);

int Neighbor_search_shiftProj(Neighborhood *neighborhood,Solution *sol, Solution *bestSol, Solution *current, int kN, int firstImp, double timeRem);
int Neighbor_search_swapJob(Neighborhood *neighborhood,Solution *sol, Solution *bestSol, Solution *current, int kN, int firstImp, double timeRem);
int Neighbor_search_shiftJob(Neighborhood *neighborhood,Solution *sol, Solution *bestSol, Solution *current, int kN, int firstImp, double timeRem);
int Neighbor_search_inv(Neighborhood *neighborhood,Solution *sol, Solution *bestSol, Solution *current, int kN, int firstImp, double timeRem);
int Neighbor_search_compactProj(Neighborhood *neighborhood,Solution *sol, Solution *bestSol, Solution *current, int kN, int firstImp, double timeRem);
int Neighbor_search_swapProj(Neighborhood *neighborhood,Solution *sol, Solution *bestSol, Solution *current, int kN, int firstImp, double timeRem);

int Neighbor_search_insertJobFILS(Neighborhood *neighborhood,Solution *sol, Solution *bestSol, Solution *current, int kN, int firstImp, double timeRem);
int Neighbor_search_swapJobFILS(Neighborhood *neighborhood,Solution *sol, Solution *bestSol, Solution *current, int kN, int firstImp, double timeRem);
int Neighbor_search_moveProj(Neighborhood *neighborhood,Solution *sol, Solution *bestSol, Solution *current, int kN, int firstImp, double timeRem);
int Neighbor_search_compactOneExtreme(Neighborhood *neighborhood,Solution *sol, Solution *bestSol, Solution *current, int kN, int firstImp, double timeRem);

/*int Neighbor_search_changeFourMode(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *current, int kN, int *nChanges,  int **nChangesJobsModes, int firstImp, double timeRem);
int Neighbor_search_changeThreeMode(Neighborhood *neighborhood,Solution *sol, Solution *bestSol, Solution *current, int kN, int *nChanges, int **nChangesJobsModes, int firstImp, double timeRem);
int Neighbor_search_changeTwoMode(Neighborhood *neighborhood,Solution *sol, Solution *bestSol, Solution *current, int kN, int *nChanges, int **nChangesJobsModes, int firstImp, double timeRem);
int Neighbor_search_changeOneMode(Neighborhood *neighborhood,Solution *sol, Solution *bestSol, Solution *current, int kN, int *nChanges, int **nChangesJobsModes, int firstImp, double timeRem);
*/

int Neighbor_search_changeFourMode(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *current, int kN, int firstImp, double timeRem);
int Neighbor_search_changeThreeMode(Neighborhood *neighborhood,Solution *sol, Solution *bestSol, Solution *current, int kN,  int firstImp, double timeRem);
int Neighbor_search_changeTwoMode(Neighborhood *neighborhood,Solution *sol, Solution *bestSol, Solution *current, int kN,  int firstImp, double timeRem);
int Neighbor_search_changeOneMode(Neighborhood *neighborhood,Solution *sol, Solution *bestSol, Solution *current, int kN,  int firstImp, double timeRem);


int Neighbor_callDetLS_Parallel( Neighborhood *neigh,
                                 Solution *sol,
                                 Solution *bestSol,
                                 Solution *currentT[],
                                 //  int *nChanges,
                                 int kN,
                                 int firstImp,
                                 double timeRem,
                                 int nThreads,
                                 Test *test);

int Neighbor_search_inv_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads);
int Neighbor_search_shiftJob_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *currentT[],int kN, int firstImp, double timeRem, int nTreads);
int Neighbor_search_swapJob_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads);
int Neighbor_search_shiftProj_parallel(Neighborhood *neighborhood,Solution *sol, Solution *bestSol, Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads);
int Neighbor_search_swapProj_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads);
int Neighbor_search_compactProj_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,  Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads);
int Neighbor_search_changeOneMode_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads);
int Neighbor_search_changeTwoMode_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads);
int Neighbor_search_changeThreeMode_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *currentT[],int kN, int firstImp, double timeRem, int nTreads);
int Neighbor_search_changeFourMode_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads);
int Neighbor_swapJobFILS_parallel(Neighborhood *neighborhood, Solution *current, Solution *sol, Solution *bestSol, int j1 );
int Neighbor_insertJobFILS_parallel(Neighborhood *neighborhood, Solution *current, Solution *sol, Solution *bestSol, int j1 );
int Neighbor_search_swapJobFILS_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads);
int Neighbor_search_insertJobFILS_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads);
int Neighbor_search_compactOnExtreme_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol,Solution *currentT[],  int kN, int firstImp, double timeRem, int nTreads);
int Neighbor_search_moveProj_parallel(Neighborhood *neighborhood, Solution *sol, Solution *bestSol, Solution *currentT[], int kN, int firstImp, double timeRem, int nTreads);

/* frees memory used by a Neighborhood */
void Neighbor_free( Neighborhood **_neighbor );

#endif // NEIGHBORHOOD_H













/* returns the numbers of iterations on neighborhood*/
//int Neighbor_getIt(Neighborhood *neig);
/* sets the vector assortment at idx position with a new neighborhood idxNeighbor*/
//void Neighbor_setIdxAssortment(Neighborhood *neighborhood, int idx, int idxNeighbor);
/* sets the vector intensity at idxNeighbor position with a new intensity priority*/
//void Neighbor_setIntensity(Neighborhood *neighborhood, int idxNeighbor, double intens);
/* returns an initialized object Neighborhood by default to MS*/
//Neighborhood *Neighbor_createMS( const Instance* inst);
/* runs the neighborhood call method for MS*/
//void Neighbor_callStocMS( Neighborhood *neighborhood, ModeSet *modeSet);
/* change one mode on ModeSet */
//void Neighbor_changeOneModeStocMS( Neighborhood *neighborhood,ModeSet *modeSet);
/* change two mode on ModeSet */
//void Neighbor_changeTwoModeStocMS( Neighborhood *neighborhood,ModeSet *modeSet );
