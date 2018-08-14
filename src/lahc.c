/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <time.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "lahc.h"
#include "vns.h"
#include <limits.h>

struct _LAHC {

    Solution *iniSol;
    Solution *bestSol;

    int lfa;
    Cost *f;
    int *nChangesModes;
    int *nChangesSequence;
    int **nTimesJobOnModes;
    int **nTimesJobOnSequence;
    int learning;
    int sideway;

    int nWOImprove;
    int nDiversification;
    int divResMode, divResJob, divTransMode, divTransJob;
    double perc, percRS;
    int nStayDiversification;
    int nCostList;
    float penaltySW;

    int nThread;
    int parallel;

    int online;
    int itUpdate;

    LearningAutomata *la;

    const struct _Instance *inst;

};

LearningAutomata *LAHC_getLA(LAHC *lahc)
{
    assert(lahc != NULL);
    return lahc->la;
}

Solution* LAHC_getBestSol(LAHC *lahc)
{

    assert(lahc != NULL);

    return lahc->bestSol;

}

void LAHC_checkArgs(LAHC *lahc, char **argv, int argc)
{

    assert(lahc != NULL);

    for(int n = 0; n< argc ; n++) {
        if (strcmp(argv[n],"-lfa") == 0) {
            n++;
            lahc->lfa =  atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-learning") == 0) {
            n++;
            lahc->learning =  atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-online") == 0) {
            n++;
            lahc->online =  atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-nwoimprove") == 0) {
            n++;
            lahc->nWOImprove =  atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-nThread") == 0) {
            n++;
            lahc->nThread = atoi( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-nDiversification") == 0) {
            n++;
            lahc->nDiversification = atoi( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-divResMode") == 0) {
            n++;
            lahc->divResMode = atoi( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-divResJob") == 0) {
            n++;
            lahc->divResJob = atoi( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-divTransMode") == 0) {
            n++;
            lahc->divTransMode = atoi( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-divTransJob") == 0) {
            n++;
            lahc->divTransJob = atoi( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-nStayDiversification") == 0) {
            n++;
            lahc->nStayDiversification = atoi( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-nCostList") == 0) {
            n++;
            lahc->nCostList = atoi( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-perc") == 0) {
            n++;
            lahc->perc = atof( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-percRS") == 0) {
            n++;
            lahc->percRS = atof( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-penaltySW") == 0) {
            n++;
            lahc->penaltySW = atof( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-parallel") == 0) {
            n++;
            lahc->parallel = atoi( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-sw") == 0) {
            n++;
            lahc->sideway = atoi( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-itUp") == 0) {
            n++;
            lahc->itUpdate = atoi( argv[n] );
            continue;
        }
    }
}

LAHC *LAHC_create( const Instance *inst, Solution* sol, Neighborhood* neighborhood, char **argv, int argc )
{
    assert(inst != NULL);
    assert(sol != NULL);
    assert(neighborhood != NULL);

    LAHC* lahc;

    ALLOCATE_INI( lahc, LAHC );

    lahc->iniSol = sol;
    lahc->bestSol = Sol_create(inst);
    Sol_cpy( lahc->bestSol, lahc->iniSol );

    /*parameter size list LAHC default*/
    lahc->lfa = 1700;

    /*parameter number of diversification, number that continuos diversifing default and percentage of diversification to residency or to transitivity memory */
    lahc->nDiversification = -1; //10000;
    lahc->nStayDiversification = -1; //20;
    lahc->perc =  0.25;
    lahc->percRS = 0.00001;


    /*parameter to activate the diversification of mode and/or modes to residency or to transitivity memory */
    lahc->divResMode = 0; //1;
    lahc->divResJob = 0; //1;
    lahc->divTransMode = 0;
    lahc->divTransJob = 0;


    /*parameter number of iterations without improvements to change intensities*/
    lahc->nWOImprove = 0;// 10000;

    /*parameter to activate learning with LA*/
    lahc->learning = 0;

    /*parameter to activate OnLine method 0 or 1*/
    lahc->online = 0;
    /*parameter that determines the number to modify intensities in OnLine method*/
    lahc->itUpdate = -1;

    /*parameter that determines the number to penalize cost list*/
    lahc->nCostList = -1;

    /*parameter to consider sideway 0 or 1 moves on OnLine method*/
    lahc->sideway = 0;
    /*parameter to determines the penalty of sideway*/
    lahc->penaltySW = 0.0;

    lahc->inst = inst;

    LAHC_checkArgs(lahc, argv, argc);

    ALLOCATE_VECTOR_INI( lahc->f, Cost, lahc->lfa );
    ALLOCATE_VECTOR_INI( lahc->nChangesModes, int, Inst_nJobs(inst) );
    ALLOCATE_VECTOR_INI( lahc->nChangesSequence, int, Inst_nJobs(inst) );

    ALLOCATE_VECTOR(lahc->nTimesJobOnModes, int*, Inst_nJobs(inst) );
    for(int i = 0; i < Inst_nJobs(inst) ; i++)
        ALLOCATE_VECTOR_INI( lahc->nTimesJobOnModes[i], int, 3);


    ALLOCATE_VECTOR(lahc->nTimesJobOnSequence, int*, Inst_nJobs(inst) );
    for(int i = 0; i < Inst_nJobs(inst) ; i++)
        ALLOCATE_VECTOR_INI( lahc->nTimesJobOnSequence[i], int, Inst_nJobs(inst));


    lahc->la = LA_create( Neighbor_nNeighborhood(neighborhood), Neighbor_getIntensity(neighborhood), argv, argc );

    return lahc;
}

void LAHC_increasingResidencyJobInMode(LAHC* lahc, Solution* current)
{

    assert(lahc != NULL);
    assert(current != NULL);

    const Job* job;
    int m;
    for(int i = 0 ; i < Inst_nJobs(lahc->inst) ; i++) {
        job = Inst_job(lahc->inst, i);
        m = Sol_getMode(current,Job_index(job));
        lahc->nTimesJobOnModes[Job_index(job)][m]++;
    }
}

void LAHC_increasingResidencyJobInSequence(LAHC* lahc, Solution* current)
{

    assert(lahc != NULL);
    assert(current != NULL);

    const Job* job;
    int s =0;
    for(int i = 0 ; i < Inst_nJobs(lahc->inst) ; i++) {
        job = Inst_job(lahc->inst, i);
        s = Sol_getPosJob(current,Job_index(job));
        lahc->nTimesJobOnSequence[Job_index(job)][s]++;
    }

}

Cost LAHC_penaltyModes(LAHC* lahc, Neighborhood* neighborhood, Cost foCurrent)
{

    assert(lahc != NULL);
    assert(neighborhood != NULL);

    Cost fator = 0, p = 0;
    fator = lahc->perc*foCurrent;

    int move = Neighbor_getLastNeigh(neighborhood), k = -1;

    if( move == changeOneMode) k = 1;
    else if( move == changeTwoMode) k = 2;
    else if( move == changeThreeMode) k = 3;
    else if( move == changeFourMode) k = 4;


    if(k>0) {

        int min=INT_MAX, max=0, sum = 0, it=0;
        double avg=0;

        for(int j = 0; j < Inst_nJobs(lahc->inst) ; j++) {
            for(int m = 0 ; m < 3 ; m++) {
                max = MAX(lahc->nTimesJobOnModes[j][m], max);
                min = MIN(lahc->nTimesJobOnModes[j][m], min);
                sum += lahc->nTimesJobOnModes[j][m];
                it++;
            }
        }

        avg = (double) sum/(double)it;

        for(int i = 0; i < k; i++) {

            int idxJob = Neighbor_getLastJob(neighborhood,i);
            int idxMode =Neighbor_getNewModes(neighborhood,i);

            p += ((lahc->nTimesJobOnModes[idxJob][idxMode] - avg) / avg)* fator;
#ifdef DEBUG
            //  printf("Job %d Mode %d Res %d \n Avg %f Fator %ld\n", idxJob, idxMode, lahc->nTimesJobOnModes[idxJob][idxMode], avg, fator);
#endif
        }
    }

#ifdef DEBUG
    //   if(p!=0)
    //     printf("Penalty Residency Job on Modes %ld\n", p);
#endif
    return p;

}

Cost LAHC_penaltyJobs(LAHC* lahc, Neighborhood* neighborhood, Cost foCurrent)
{

    assert(lahc != NULL);
    assert(neighborhood != NULL);

    Cost fator = 0, p = 0;
    fator = lahc->percRS*foCurrent;

    int min=INT_MAX, max=0, sum = 0, it=0;
    double avg=0;

    for(int j = 0; j < Inst_nJobs(lahc->inst) ; j++) {
        for(int s = 0 ; s < Inst_nJobs(lahc->inst) ; s++) {
            max = MAX(lahc->nTimesJobOnSequence[j][s], max);
            min = MIN(lahc->nTimesJobOnSequence[j][s], min);
            sum += lahc->nTimesJobOnSequence[j][s];
            it++;
        }
    }

    avg = (double) sum/(double)it;

    for(int j = 0; j < Neighbor_getContLastJ(neighborhood) ; j++) {
        int idxJob  = Neighbor_getLastJobModify(neighborhood,j);
        int idxSequence = Neighbor_getPosLastJobModify(neighborhood,j);
#ifdef DEBUG
        // printf("Job %d PosSequence %d Res %d Avg %f Factor %ld \n", idxJob, idxSequence, vns->nTimesJobOnSequence[idxJob][idxSequence], avg, fator);
#endif
        p += ((lahc->nTimesJobOnSequence[idxJob][idxSequence] - avg) / avg)* fator;
    }

#ifdef DEBUG
    //      if(p!=0)
    //      printf("Penalty Residency Job on Sequence %ld\n", p);
#endif

    return p;

}

void LAHC_increasingTransitivityOfModes(LAHC* lahc, Neighborhood* neighborhood)
{

    assert(lahc != NULL);
    assert(neighborhood != NULL);

    if( Neighbor_getLastNeigh(neighborhood) == changeOneMode)
        lahc->nChangesModes[Neighbor_getLastJob(neighborhood, 0)]++;
    else if ( Neighbor_getLastNeigh(neighborhood) == changeTwoMode) {
        lahc->nChangesModes[Neighbor_getLastJob(neighborhood, 0)]++;
        lahc->nChangesModes[Neighbor_getLastJob(neighborhood, 1)]++;
    } else if ( Neighbor_getLastNeigh(neighborhood) == changeThreeMode) {
        lahc->nChangesModes[Neighbor_getLastJob(neighborhood, 0)]++;
        lahc->nChangesModes[Neighbor_getLastJob(neighborhood, 1)]++;
        lahc->nChangesModes[Neighbor_getLastJob(neighborhood, 2)]++;
    } else if ( Neighbor_getLastNeigh(neighborhood) == changeFourMode) {
        lahc->nChangesModes[Neighbor_getLastJob(neighborhood, 0)]++;
        lahc->nChangesModes[Neighbor_getLastJob(neighborhood, 1)]++;
        lahc->nChangesModes[Neighbor_getLastJob(neighborhood, 2)]++;
        lahc->nChangesModes[Neighbor_getLastJob(neighborhood, 3)]++;
    }
}

Cost LAHC_penaltyTransModes(LAHC* lahc, Neighborhood* neighborhood, Cost foCurrent)
{

    assert(lahc != NULL);
    assert(neighborhood != NULL);

    Cost fator = 0, p = 0;
    fator = lahc->perc*foCurrent;

    int move = Neighbor_getLastNeigh(neighborhood), k = -1;

    if( move == changeOneMode) k = 1;
    else if( move == changeTwoMode) k = 2;
    else if( move == changeThreeMode) k = 3;
    else if( move == changeFourMode) k = 4;


    if(k>0) {

        int min=INT_MAX, max=0, sum = 0, it=0;
        double avg=0;

        for(int j = 0; j < Inst_nJobs(lahc->inst) ; j++) {
            max = MAX(lahc->nChangesModes[j], max);
            min = MIN(lahc->nChangesModes[j], min);
            sum += lahc->nChangesModes[j];
            it++;
        }

        avg = (double) sum/(double)it;

        for(int i = 0; i < k; i++) {
            int idxJob = Neighbor_getLastJob(neighborhood, i);
#ifdef DEBUG
            //  printf("Job %d Trans %d \n Avg %f Factor %ld \n", idxJob, vns->nChangesModes[idxJob], avg, fator);
#endif
            p += ((lahc->nChangesModes[idxJob]- avg) / avg)*fator;
        }
    }

#ifdef DEBUG
    //   if(p!=0)
    //     printf("Penalty transitivity modes %ld\n", p);
#endif
    return p;

}

void LAHC_increasingTransitivityOfSequence(LAHC* lahc, Neighborhood* neighborhood)
{

    assert(lahc != NULL);
    assert(neighborhood != NULL);

    for(int i = 0 ; i < Neighbor_getContLastJ(neighborhood) ; i++)
        lahc->nChangesSequence[Neighbor_getLastJobModify(neighborhood,i)]++;

}

Cost LAHC_penaltyTransSequence(LAHC* lahc, Neighborhood* neighborhood, Cost foCurrent)
{

    assert(lahc != NULL);
    assert(neighborhood != NULL);

    Cost fator = 0, p = 0;
    fator = lahc->perc*foCurrent;

    int min=INT_MAX, max=0, sum = 0, it=0;
    double avg=0;

    for(int j = 0; j < Inst_nJobs(lahc->inst) ; j++) {
        max = MAX(lahc->nChangesSequence[j], max);
        min = MIN(lahc->nChangesSequence[j], min);
        sum += lahc->nChangesSequence[j];
        it++;
    }

    avg = (double) sum/(double)it;

    for(int j = 0; j < Neighbor_getContLastJ(neighborhood) ; j++) {
        int idxJob = Neighbor_getLastJobModify(neighborhood,j);
#ifdef DEBUG
        // printf("Job %d Trans %d \n Avg %f Factor %ld \n", idxJob, lahc->nChangesSequence[idxJob], avg, fator);
#endif
        p += ((lahc->nChangesSequence[idxJob]- avg) / avg)* fator;
    }

#ifdef DEBUG
    //   if(p!=0)
    //   printf("Penalty Transitivity Sequence %ld\n", p);
#endif
    return p;

}

void LAHC_run(LAHC *lahc, Neighborhood* neighborhood, double timeRem, Test *test, char* nameInst)
{

    assert(lahc != NULL);
    assert(neighborhood != NULL);

    long double _time = 0, _timeP =0, oldT = omp_get_wtime();
    long double timeRemP = 10;
    long double _currentTime = 0;

    clock_t tStart = clock();

    int I = 0, v=0, valid=0;

    int nItLearning=0, nIt = 0, nItDiversification = 0, nItStay = -1, nItCostList = 0, nItUp = 0;

    Cost fo = Sol_getCost(lahc->iniSol);
    for(int l = 0; l < LAHC_getLfa(lahc); l++)
        LAHC_setF(lahc,l,fo);

    Solution *current = Sol_create(lahc->inst);
    Sol_cpy(current,lahc->iniSol);


    if( lahc->learning)
        LA_printProbabilities(lahc->la, lahc->penaltySW);

    _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
    if( lahc->online) {
        //Neighbor_printProbabilitiesTime( neighborhood, lahc->sideway, lahc->itUpdate, lahc->penaltySW, nameInst, _time);
        Neighbor_printProbabilities( neighborhood, lahc->sideway, lahc->itUpdate, lahc->penaltySW, nameInst);
    }


    while(_time < timeRem ) {

        _time = ( (long double)( clock()-tStart ));//CLOCKS_PER_SEC );

        /*increasing the residency of a job in a mode and a position in asequence*/
        if(lahc->divResMode)
            LAHC_increasingResidencyJobInMode(lahc, current);
        if(lahc->divResJob)
            LAHC_increasingResidencyJobInSequence(lahc, current);

        valid = Neighbor_callStocLS( neighborhood, lahc->iniSol, lahc->bestSol, current, lahc->la, lahc->learning);
        if(valid) {
            if(Neighbor_getLastNeigh(neighborhood) != seqSwapJobFILS || Neighbor_getLastNeigh(neighborhood) != seqInsertJobFILS)
                Sol_rebuild_opt(current, lahc->iniSol);
            //  Neighbor_checkMoveSW(neighborhood, current, lahc->iniSol);
        }


        if(lahc->learning) {
            _timeP = omp_get_wtime()-oldT;
            if ( _timeP >= timeRemP ) {
                oldT =  omp_get_wtime();
                LA_setTime(lahc->la);
                LA_printProbabilities(lahc->la,lahc->penaltySW);
            }
        }


        if(lahc->online) {
            _timeP = omp_get_wtime()-oldT;
            if ( _timeP >= timeRemP ) {
                oldT =  omp_get_wtime();
                Neighbor_setTime(neighborhood);
                Neighbor_printProbabilities( neighborhood, lahc->sideway, lahc->itUpdate,lahc->penaltySW, nameInst );
            }
        }

        v = I % LAHC_getLfa(lahc);

        Cost fo = Sol_getCost(current);

        /*penalizing residency*/
        if( nItDiversification == lahc->nDiversification || (nItStay < lahc->nStayDiversification && nItStay >= 0)) {

            if(lahc->divResMode) {
#ifdef DEBUG
                //             printf("\nPenalizing residency of modes...\n");
                fflush(stdout);
                fflush(stdin);
#endif
                fo += LAHC_penaltyModes(lahc, neighborhood, Sol_getCost(current) );
            }
            if(lahc->divResJob) {
#ifdef DEBUG
                //           printf("\nPenalizing residency of jobs...\n");
                fflush(stdout);
                fflush(stdin);
#endif
                fo += LAHC_penaltyJobs(lahc, neighborhood, Sol_getCost(current) );
            }
            if(lahc->divTransMode) {
#ifdef DEBUG
                //         printf("\nPenalizing transitivity of modes...\n");
                fflush(stdout);
                fflush(stdin);
#endif
                fo += LAHC_penaltyTransModes(lahc, neighborhood, Sol_getCost(current) );
            }
            if(lahc->divTransJob) {
                fo += LAHC_penaltyTransSequence(lahc, neighborhood, Sol_getCost(current) );
#ifdef DEBUG
                //       printf("\nPenalizing transitivity of jobs...\n");
                fflush(stdout);
                fflush(stdin);
#endif

            }

            nItStay++;

            if(nItStay == lahc->nStayDiversification-1) {
                nItDiversification = 0;
                nItStay=-1;
            }
        }

        //   printf("Value %ld", Sol_getCost(current));
        //Sol_getCost(current)

        if(  fo <= LAHC_getF(lahc, v) || fo < Sol_getCost(lahc->iniSol) ||  fo <  Sol_getCost(lahc->bestSol) ) {
            /*increasing the transitivity of modes*/
            if(lahc->divTransMode)
                LAHC_increasingTransitivityOfModes(lahc, neighborhood);
            if(lahc->divTransJob)
                LAHC_increasingTransitivityOfSequence(lahc, neighborhood);

            /*incrementing improvement and time to calculate the intensities of jobs by online method*/
            _currentTime = ( (long double)( clock()-tStart ));//CLOCKS_PER_SEC );
            if(lahc->online) {
                if(lahc->sideway) {
                    if(Sol_getCost(current) < Sol_getCost(lahc->iniSol)) {
                        int n = Neighbor_getLastNeigh(neighborhood);
                        Neighbor_incrementI( neighborhood, n-1);
                        Neighbor_setTI(neighborhood, n-1, (long double) (_currentTime - _time));
                    } else if( Sol_getCost(current) == Sol_getCost(lahc->iniSol)) {
                        if(Sol_getSolutionHash(current)!=Sol_getSolutionHash(lahc->iniSol)) {
                            //if(Neighbor_getMoveSW(neighborhood)) {
                            //  printf("Truly sideway, incrementing EQ: Hash Current %u  Hash Init %u\n\n", Sol_getSolutionHash(current), Sol_getSolutionHash(lahc->iniSol));
                            int n =  Neighbor_getLastNeigh(neighborhood);
                            Neighbor_incrementEQ(neighborhood, n-1);
                            Neighbor_setTE(neighborhood, n-1, (long double) (_currentTime - _time));
                        } else {
                            // printf("Fake sideway, incrementing TIV\n\n");
                            int n =  Neighbor_getLastNeigh(neighborhood);
                            Neighbor_setTIV(neighborhood, n-1, (long double) (_currentTime - _time));
                        }
                    }
                } else {
                    if(Sol_getCost(current) < Sol_getCost(lahc->iniSol)) {
                        int n =  Neighbor_getLastNeigh(neighborhood);
                        Neighbor_incrementI( neighborhood, n-1);
                        Neighbor_setTI(neighborhood, n-1, (long double) (_currentTime - _time));
                    } else {
                        int n =  Neighbor_getLastNeigh(neighborhood);
                        Neighbor_setTIV(neighborhood, n-1, (long double) (_currentTime - _time));
                    }
                }

            }


            /*update intensities of jobs by learning method*/
            if(lahc->learning) {
                if(lahc->sideway) {
                    if(Sol_getCost(current) < Sol_getCost(lahc->iniSol))
                        LA_update( lahc->la, 1 );
                    else if( Sol_getCost(current) == Sol_getCost(lahc->iniSol))
                        LA_update( lahc->la, lahc->penaltySW );
                } else {
                    if(Sol_getCost(current) < Sol_getCost(lahc->iniSol))
                        LA_update( lahc->la, 1 );
                }
            }

            Sol_cpy( lahc->iniSol,current );

            if( Sol_getCost(lahc->iniSol) < Sol_getCost(lahc->bestSol)) {
#ifdef DEBUG
                _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
                printf("\nLAHC %ld %.2Lf ", Sol_getCost(lahc->bestSol), _time);
                fflush(stdout);
                fflush(stdin);
#endif
                Sol_cpy( lahc->bestSol, lahc->iniSol );
                nIt = 0;
                nItLearning = 0;
                nItCostList = 0;
                nItDiversification = 0;

            }

        } else {

            nIt++;
            if(nIt == lahc->nWOImprove) {
#ifdef DEBUG
                printf("\nUpdating intensities offline...\n ");
                fflush(stdout);
                fflush(stdin);
#endif
                Neighbor_getUpdatesIntensity(neighborhood);
                nIt = 0;
            }

            _currentTime = ( (double)( clock()-tStart ));//CLOCKS_PER_SEC );
            if(lahc->online) {
                int n =  Neighbor_getLastNeigh(neighborhood);
                Neighbor_setTIV(neighborhood, n-1, (long double) (_currentTime - _time));
            }

            /*update Intensities Online*/
            if(nItUp == lahc->itUpdate) {
#ifdef DEBUG
                printf("\nUpdate Intensities online...");
#endif

                Neighbor_setF(neighborhood);
                Neighbor_normF(neighborhood);


                Neighbor_updatesIntensities(neighborhood, lahc->sideway);

                printf("\nIntens: ");
                for(int i=0; i<Neighbor_nNeighborhood(neighborhood); ++i)
                    printf("%.3Lf ", Neighbor_getIdIntensity(neighborhood, i));
                fflush(stdout);

                Neighbor_clearImpTime(neighborhood);
                _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
                //if( lahc->online)
                //Neighbor_printProbabilitiesTime( neighborhood, lahc->sideway, lahc->itUpdate, lahc->penaltySW, nameInst, _time);

                nItUp = 0;
            }

            /*update Intensities Learning*/
            if( lahc->learning) {
                nItLearning++;

                if(lahc->sideway) {
                    if( Sol_getCost(current) == Sol_getCost(lahc->iniSol))
                        LA_update( lahc->la, lahc->penaltySW );
                    else LA_update( lahc->la, 0 );
                } else
                    LA_update( lahc->la, 0 );

                /*resets the learning automata when it reaches the iteration limit*/
                if(nItLearning >= LA_getResetInterval(lahc->la)) {
#ifdef DEBUG
                    printf("\nReset Interval %d, %d\n", nItLearning, LA_getResetInterval(lahc->la));
                    fflush(stdout);
                    fflush(stdin);
#endif
                    LA_reset( lahc->la );
                    nItLearning=0;
                }
            }


            /*update cost of size list LFA of LAHC*/
            nItCostList++;
            if( nItCostList == lahc->nCostList) {
#ifdef DEBUG
                printf("\nUpdating Cost List...\n");
                fflush(stdout);
                fflush(stdin);
#endif
                LAHC_updateF(lahc);
                nItCostList= 0;
            }

            nItDiversification++;

            Sol_cpy( current, lahc->iniSol );
        }


        Neighbor_setNullLastJobModify(neighborhood);

        LAHC_setF(lahc,v, Sol_getCost(lahc->iniSol));
        nItUp++;
        I++;

        _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC ); //omp_get_wtime()-startT;
    }

    /*if(lahc->online) {
        Neighbor_printProbabilities( neighborhood, lahc->sideway, lahc->itUpdate,lahc->penaltySW, nameInst);
    }*/

    if(lahc->online) {
        //_time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
        //Neighbor_printProbabilitiesTime( neighborhood, lahc->sideway, lahc->itUpdate, lahc->penaltySW, nameInst, _time);
        if(Neighbor_getTimePrint(neighborhood) < 300) {
            Neighbor_setTime(neighborhood);
            Neighbor_printProbabilities( neighborhood, lahc->sideway, lahc->itUpdate,lahc->penaltySW, nameInst );
        }
    }

    Sol_cpy( lahc->iniSol, lahc->bestSol );

    Sol_free(&current);
}

void LAHC_run_nIt(LAHC *lahc, Neighborhood* neighborhood, double timeRem, int nIterations, Test *test, char* nameInst)
{

    assert(lahc != NULL);
    assert(neighborhood != NULL);

    long double _time = 0;

    clock_t tStart = clock();

    int I = 0, v=0, valid=0;

    Cost fo = Sol_getCost(lahc->iniSol);
    for(int l = 0; l < LAHC_getLfa(lahc); l++)
        LAHC_setF(lahc,l,fo);

    Solution *current = Sol_create(lahc->inst);
    Sol_cpy(current,lahc->iniSol);

    _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
    int nIte =0;

    while(nIte <= nIterations && _time < timeRem) {

        _time = ( (long double)( clock()-tStart ));//CLOCKS_PER_SEC );

        valid = Neighbor_callStocLS( neighborhood, lahc->iniSol, lahc->bestSol, current, lahc->la, lahc->learning);
        if(valid) {
            Sol_rebuild_opt(current, lahc->iniSol);
            //   Neighbor_checkMoveSW(neighborhood, current, lahc->iniSol);
        }

        v = I % LAHC_getLfa(lahc);

        Cost fo = Sol_getCost(current);

        if(  fo <= LAHC_getF(lahc, v) || fo < Sol_getCost(lahc->iniSol) ||  fo <  Sol_getCost(lahc->bestSol) ) {
            /*increasing the transitivity of modes*/
            Sol_cpy( lahc->iniSol,current );

            if( Sol_getCost(lahc->iniSol) < Sol_getCost(lahc->bestSol)) {
#ifdef DEBUG
                _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
                printf("\nLAHC %ld %.2Lf ", Sol_getCost(lahc->bestSol), _time);
                fflush(stdout);
                fflush(stdin);
#endif
                Sol_cpy( lahc->bestSol, lahc->iniSol );
                nIte = 0;
            }

        } else {

            nIte++;

            Sol_cpy( current, lahc->iniSol );
        }


        Neighbor_setNullLastJobModify(neighborhood);

        LAHC_setF(lahc,v, Sol_getCost(lahc->iniSol));
        I++;

        _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC ); //omp_get_wtime()-startT;
    }

    Sol_cpy( lahc->iniSol, lahc->bestSol );

    Sol_free(&current);
}

void LAHC_run_vns(LAHC *lahc, VNS* vns, Neighborhood* neighborhood, double timeRem,  int nIterations, char* nameInst)
{

    assert(lahc != NULL);
    assert(neighborhood != NULL);

    double _time = 0;

    clock_t tStart = clock();

    int I = 0, v=0, valid=0;
    int nIt = 0, nItUp = 0;

    Cost fo = Sol_getCost(lahc->iniSol);

    for(int l = 0; l < LAHC_getLfa(lahc); l++)
        LAHC_setF(lahc,l,fo);

    Solution *current = Sol_create(lahc->inst);
    Sol_cpy(current,lahc->iniSol);

    _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
    int nIte =0;

    while(nIte <= nIterations && _time < timeRem) {

        /*increasing the residency of a job in a mode and a position in asequence*/
        if(VNS_getDivRM(vns))
            VNS_increasingResidencyJobInMode(vns, current);
        if(VNS_getDivRJ(vns))
            VNS_increasingResidencyJobInSequence(vns, current);

        valid = Neighbor_callStocLS( neighborhood, lahc->iniSol, lahc->bestSol, current, lahc->la, lahc->learning);
        if(valid) {

            if(Neighbor_getLastNeigh(neighborhood) != seqSwapJobFILS || Neighbor_getLastNeigh(neighborhood) != seqInsertJobFILS)
                Sol_rebuild_opt(current, lahc->iniSol);


            v = I % LAHC_getLfa(lahc);

            Cost fo = Sol_getCost(current);


            if(  fo <= LAHC_getF(lahc, v) || fo < Sol_getCost(lahc->iniSol) ||  fo <  Sol_getCost(lahc->bestSol) ) {
                /*increasing the transitivity of modes*/
                if(VNS_getDivTM(vns))
                    VNS_increasingTransitivityOfModes(vns, neighborhood);
                if(VNS_getDivTJ(vns))
                    VNS_increasingTransitivityOfSequence(vns, neighborhood);


                Sol_cpy( lahc->iniSol,current );

                if( Sol_getCost(lahc->iniSol) < Sol_getCost(lahc->bestSol)) {
#ifdef DEBUG
                    _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
                    printf("\nLAHC %ld %.2f ", Sol_getCost(lahc->bestSol), _time);
                    fflush(stdout);
                    fflush(stdin);
#endif
                    Sol_cpy( lahc->bestSol, lahc->iniSol );
                    nIt = 0;
                    nIte = 0;
                }

            } else {

                nIt++;
                nIte++;
                if(nIt == lahc->nWOImprove) {
#ifdef DEBUG
                    printf("\nUpdating intensities offline...\n ");
                    fflush(stdout);
                    fflush(stdin);
#endif
                    Neighbor_getUpdatesIntensity(neighborhood);
                    nIt = 0;
                }

                Sol_cpy( current, lahc->iniSol );
            }
        } else
            Sol_cpy( current, lahc->iniSol );


        Neighbor_setNullLastJobModify(neighborhood);

        LAHC_setF(lahc,v, Sol_getCost(lahc->iniSol));
        nItUp++;
        I++;

        _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC ); //omp_get_wtime()-startT;
    }

    Sol_cpy( lahc->iniSol, lahc->bestSol );

    Sol_free(&current);
}

void LAHC_run_parallel(LAHC *lahc, Neighborhood* neighborhood, double timeRem, Solution *globalBest, Test *test)
{

    assert(lahc != NULL);
    assert(neighborhood != NULL);

    double _time = 0;
    double startT = omp_get_wtime();
    int I = 0, v=0, valid=0;

    int nItLearning=0, nIt = 0, nItDiversification =0, nItCostList=0;

    Cost fo = Sol_getCost(lahc->iniSol);
    for(int l = 0; l < LAHC_getLfa(lahc); l++)
        LAHC_setF(lahc,l,fo);

    Solution *current = Sol_create(lahc->inst);
    Sol_cpy(current,lahc->iniSol);

    /* synchronization of best solution between threads */
    #pragma omp critical
    {
        if ( Sol_getCost(lahc->bestSol) < Sol_getCost(globalBest) )
            Sol_cpy(globalBest,lahc->bestSol);
        else if ( Sol_getCost(globalBest) < Sol_getCost(current) )
            Sol_cpy(current,globalBest);
    }


    while(_time < timeRem )  {


        Test_setCurrentFO(test, Sol_getCost(lahc->bestSol));

        /*increasing the residency of a job in a mode*/
        LAHC_increasingResidencyJobInMode(lahc, current);
        LAHC_increasingResidencyJobInSequence(lahc, current);

        valid = Neighbor_callStocLS( neighborhood, lahc->iniSol, lahc->bestSol, current, lahc->la, lahc->learning);
        if(valid)
            Sol_rebuild_opt(current, lahc->iniSol);

        v = I % LAHC_getLfa(lahc);

        Cost fo = Sol_getCost(current);

        /*penalizing residency*/
        if( nItDiversification == lahc->nDiversification) {
#ifdef DEBUG
            printf("\nPenalizing...\n");
#endif
            /*valid = Neighbor_random_changeNModes(neighborhood, lahc->iniSol,  current, 4, lahc->nTimesJobOnModes );
            if(valid)
                Sol_rebuild_opt(current, lahc->iniSol);
            */
            fo += LAHC_penaltyModes(lahc, neighborhood, Sol_getCost(current) );
            //printf("Current FO %ld, foPenalized %ld, Penalty %ld", Sol_getCost(current),  fo, LAHC_penalty(lahc, neighborhood, Sol_getCost(current)));
            nItDiversification= 0;

        }

        if( fo <= LAHC_getF(lahc, v) || fo <= Sol_getCost(lahc->iniSol)) {

            /*increasing the transitivity of modes*/
            LAHC_increasingTransitivityOfModes(lahc, neighborhood);

            if(lahc->learning) {
                if( Sol_getCost(current) < LAHC_getF(lahc, v) || Sol_getCost(current) < Sol_getCost(lahc->iniSol))
                    LA_update( lahc->la, 1 );
                else
                    LA_update( lahc->la, 0.1 );
            }

            Sol_cpy( lahc->iniSol,current );
            if( Sol_getCost(lahc->iniSol) < Sol_getCost(lahc->bestSol)) {

#ifdef DEBUG
                _time = omp_get_wtime()-startT;
                printf("\n%ld %.2f %d ", Sol_getCost(lahc->bestSol), _time, omp_get_thread_num());
#endif
                Sol_cpy( lahc->bestSol, lahc->iniSol );
                /* checking best improvement */
                #pragma omp critical
                {
                    if (Sol_getCost(lahc->bestSol)  < Sol_getCost(globalBest) ) {
                        Sol_cpy(globalBest,lahc->bestSol);
#ifdef DEBUG
                        printf("\nThread %d: improved globalBest to: %ld. %.3f seconds. %d iterations. \n", omp_get_thread_num(), Sol_getCost(globalBest), omp_get_wtime()-startT, I);
                        fflush(stdout);
#endif
                    }
                }

                nIt = 0;
                nItLearning = 0;
                nItCostList = 0;
                nItDiversification = 0;

            }
        } else {

            nIt++;
            if(nIt == lahc->nWOImprove) {

#ifdef DEBUG
                printf("\nUpdating intensities...\n ");
#endif
                Neighbor_getUpdatesIntensity(neighborhood);
                /* updating LAHC bestSol after iterations without improvements */
                #pragma omp critical
                {
                    int th_id = omp_get_thread_num();
                    if (Sol_getCost(globalBest) < Sol_getCost(lahc->bestSol)  ) {
#ifdef DEBUG
                        printf("\nThread: %d updating LAHC bestSol %ld to globalBest %ld \n", th_id,  Sol_getCost(lahc->bestSol),Sol_getCost(globalBest) );
#endif
                        Sol_cpy(lahc->bestSol,globalBest);
                        Sol_cpy(lahc->iniSol,lahc->bestSol);
                    }
                }
                nIt = 0;
            }

            if( lahc->learning) {
                nItLearning++;
                LA_update( lahc->la, 0 );
                //resets the learning automata when it reaches the iteration limit
                if(nItLearning >= LA_getResetInterval(lahc->la)) {
#ifdef DEBUG
                    printf("\nReset learning interval %d, %d\n", nItLearning, LA_getResetInterval(lahc->la));
#endif
                    LA_reset( lahc->la );
                    nItLearning=0;
                }
            }

            nItCostList++;
            if( nItCostList == lahc->nCostList) {
#ifdef DEBUG
                printf("\nUpdating Cost List...\n");
#endif
                LAHC_updateF(lahc);
                nItCostList = 0;
            }

            nItDiversification++;

            Sol_cpy( current, lahc->iniSol );
        }

        Neighbor_setNullLastJobModify(neighborhood);

        LAHC_setF(lahc,v, Sol_getCost(lahc->iniSol));
        I++;

        _time = omp_get_wtime()-startT;
        Test_callTest(test, Test_getCurrentNeigh(test), Test_getCurrentTime(test), Sol_getCost(lahc->bestSol));

    }

    Sol_cpy( lahc->iniSol, lahc->bestSol );

    Sol_free(&current);
}

void LAHC_updateF(LAHC *lahc)
{

    assert(lahc!=NULL);
    for(int i = 0 ; i< lahc->lfa ; i++ )
        lahc->f[i] += (lahc->f[i] * lahc->perc);

}

void LAHC_setF(LAHC *lahc, int idxF, Cost value)
{

    assert(lahc!=NULL);
    assert(idxF < lahc->lfa);

    lahc->f[idxF] = value;

}

int LAHC_getLfa(LAHC *lahc)
{
    assert(lahc!=NULL);

    return lahc->lfa;
}

int LAHC_getItUp(LAHC *lahc)
{
    assert(lahc!=NULL);

    return lahc->itUpdate;
}

float LAHC_getPSW(LAHC *lahc)
{
    assert(lahc!=NULL);

    return lahc->penaltySW;
}

int LAHC_getRJ(LAHC *lahc)
{
    assert(lahc!=NULL);

    return lahc->divResJob;
}

int LAHC_getRM(LAHC *lahc)
{
    assert(lahc!=NULL);

    return lahc->divResMode;
}

int LAHC_getTJ(LAHC *lahc)
{
    assert(lahc!=NULL);

    return lahc->divTransJob;
}

int LAHC_getTM(LAHC *lahc)
{
    assert(lahc!=NULL);

    return lahc->divTransMode;
}

int LAHC_getF(LAHC *lahc, int idxF)
{

    assert(lahc!=NULL);
    assert(idxF >=0 && idxF < lahc->lfa);

    return lahc->f[idxF];

}

int LAHC_getNChangesModes(LAHC *lahc, int idx)
{

    assert(lahc!=NULL);

    return lahc->nChangesModes[idx];

}

int LAHC_getNDiversification(LAHC *lahc)
{
    assert(lahc!=NULL);
    return lahc->nDiversification;
}

int LAHC_getNStayDiversification(LAHC *lahc)
{
    assert(lahc!=NULL);
    return lahc->nStayDiversification;
}

int LAHC_getNWOImprove(LAHC *lahc)
{
    assert(lahc!=NULL);
    return lahc->nWOImprove;
}

int LAHC_getSW(LAHC *lahc)
{
    assert(lahc!=NULL);
    return lahc->sideway;
}

int LAHC_getNCostList(LAHC *lahc)
{
    assert(lahc!=NULL);
    return lahc->nCostList;
}

int LAHC_getNThread(LAHC *lahc)
{
    assert(lahc!=NULL);
    return lahc->nThread;
}

double LAHC_getPerc(LAHC *lahc)
{
    assert(lahc!=NULL);
    return lahc->perc;
}

void LAHC_free( LAHC **_lahc )
{

    LAHC *lahc = *_lahc;

    for(int i = 0 ; i < Inst_nJobs(Sol_inst(lahc->iniSol)); i++)
        free( lahc->nTimesJobOnModes[i] );
    free( lahc->nTimesJobOnModes );

    for(int i = 0; i < Inst_nJobs(Sol_inst(lahc->iniSol)) ; i++)
        free(lahc->nTimesJobOnSequence[i]);
    free( lahc->nTimesJobOnSequence );


    free( lahc->f );

    free( lahc->nChangesModes );
    free( lahc->nChangesSequence );


    Sol_free( &lahc->iniSol );
    Sol_free( &lahc->bestSol );
    LA_free( &lahc->la );

    free( lahc );

    *_lahc = NULL;

}


