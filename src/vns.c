/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <time.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <omp.h>
#include "test.h"
#include "vns.h"
#include "vnd.h"
#include "lahc.h"

struct _VNS {

    Solution *iniSol;
    Solution *bestSol;

    int nMoves;
    int sizeSamplingShake;
    LearningAutomata * la;
    int itRNA;
    int itLAHC;
    int lfa;


    int divResMode, divResJob, divTransMode, divTransJob;
    int *nChangesModes;
    int *nChangesSequence;
    int **nTimesJobOnModes;
    int **nTimesJobOnSequence;

    double perc;
    double percRS;

    const struct _Instance *inst;

};

int VNS_getDivRM(VNS *vns)
{

    assert(vns!=NULL);

    return vns->divResMode;
}


double VNS_getPerc(VNS *vns)
{

    assert(vns!=NULL);

    return vns->perc;
}

double VNS_getPercRS(VNS *vns)
{

    assert(vns!=NULL);

    return vns->percRS;
}

int VNS_getLfa(VNS *vns)
{

    assert(vns!=NULL);

    return vns->lfa;
}

int VNS_getDivRJ(VNS *vns)
{

    assert(vns!=NULL);

    return vns->divResJob;
}

int VNS_getDivTJ(VNS *vns)
{

    assert(vns!=NULL);

    return vns->divTransJob;
}

int VNS_getDivTM(VNS *vns)
{

    assert(vns!=NULL);

    return vns->divTransMode;
}

int VNS_getItLAHC(VNS *vns)
{

    assert(vns!=NULL);

    return vns->itLAHC;
}

int VNS_getItRNA(VNS *vns)
{

    assert(vns!=NULL);

    return vns->itRNA;
}

void VNS_checkArgs(VNS *vns, char **argv, int argc)
{

    assert( vns != NULL);

    int value;

    for(int n = 0 ; n < argc ; n++) {
        if (strcmp(argv[n],"-nMoves") == 0) {
            n++;
            value = atoi( argv[n] );
            vns->nMoves = value;
            continue;
        }
        if (strcmp(argv[n],"-nSize") == 0) {
            n++;
            value = atoi( argv[n] );
            vns->sizeSamplingShake = value;
            continue;
        }
        if (strcmp(argv[n],"-itRNA") == 0) {
            n++;
            value = atoi( argv[n] );
            vns->itRNA = value;
            continue;
        }
        if (strcmp(argv[n],"-itLAHC") == 0) {
            n++;
            value = atoi( argv[n] );
            vns->itLAHC = value;
            continue;
        }

        if (strcmp(argv[n],"-lfa") == 0) {
            n++;
            value = atoi( argv[n] );
            vns->lfa = value;
            continue;
        }
        if (strcmp(argv[n],"-divResMode") == 0) {
            n++;
            vns->divResMode = atoi( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-divResJob") == 0) {
            n++;
            vns->divResJob = atoi( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-divTransMode") == 0) {
            n++;
            vns->divTransMode = atoi( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-divTransJob") == 0) {
            n++;
            vns->divTransJob = atoi( argv[n] );
            continue;
        }

        if (strcmp(argv[n],"-perc") == 0) {
            n++;
            vns->perc = atof( argv[n] );
            continue;
        }
        if (strcmp(argv[n],"-percRS") == 0) {
            n++;
            vns->percRS = atof( argv[n] );
            continue;
        }

        //   if (strcmp(argv[n],"-allNeighbor") == 0) {
        //     n++;
        //   value = atoi( argv[n] );
        // vns->allNeighbor = value;
        //continue;
        //}

    }

}

VNS *VNS_create( const Instance *inst, Solution *sol, Neighborhood *neighborhood, char **argv, int argc )
{
    assert( inst != NULL);
    assert( neighborhood != NULL);
    assert( sol != NULL);

    VNS *vns;

    ALLOCATE_INI( vns, VNS );
    //  ALLOCATE_VECTOR_INI( vns->nChangesModes, int, Inst_nJobs( inst ) );
    /*ALLOCATE_VECTOR( vns->nTimesJobOnModes, int*, Inst_nJobs( inst ) );
    for(int i =0 ; i < Inst_nJobs( inst ); i++)
        ALLOCATE_VECTOR_INI( vns->nTimesJobOnModes[i], int, 3);
    */
    vns->bestSol = sol;
    vns->iniSol = Sol_create(inst);
    Sol_cpy( vns->iniSol, vns->bestSol );
    vns->la = LA_create( Neighbor_nNeighborhood(neighborhood), Neighbor_getIntensity(neighborhood), argv, argc );

    vns->inst = inst;
    vns->nMoves = 14;
    vns->itRNA = 0;
    vns->itLAHC = 0;


    /*parameter to activate the diversification of mode and/or modes to residency or to transitivity memory */
    vns->divResMode = 0; //1;
    vns->divResJob = 0; //1;
    vns->divTransMode = 0;
    vns->divTransJob = 0;
    vns->perc = 0.25;
    vns->percRS = 0.00001;

    VNS_checkArgs(vns,argv, argc);

    ALLOCATE_VECTOR_INI( vns->nChangesModes, int, Inst_nJobs(inst) );
    ALLOCATE_VECTOR_INI( vns->nChangesSequence, int, Inst_nJobs(inst) );

    ALLOCATE_VECTOR(vns->nTimesJobOnModes, int*, Inst_nJobs(inst) );
    for(int i = 0; i < Inst_nJobs(inst) ; i++)
        ALLOCATE_VECTOR_INI( vns->nTimesJobOnModes[i], int, 3);


    ALLOCATE_VECTOR(vns->nTimesJobOnSequence, int*, Inst_nJobs(inst) );
    for(int i = 0; i < Inst_nJobs(inst) ; i++)
        ALLOCATE_VECTOR_INI( vns->nTimesJobOnSequence[i], int, Inst_nJobs(inst));

    return vns;
}

//vns with  Sequential Variable Neighborhood Descent k = 14, // if we want a Reduced Variable Neighborhood Search just remove vnd VNS_run_reduced
void VNS_run_general_vnd(VNS *vns, Neighborhood* neighborhood, double timeRem, int nNeighbor, int firstImprovement, Test *test, char **argv, int argc )
{
    assert( vns != NULL);
    assert( neighborhood != NULL);

    double _time = 0;
    clock_t tStart = clock();

    int nMoves = VNS_getNMoves(vns);
    int it = 1, valid = 0;

    _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

    Solution *current = Sol_create(Sol_inst( vns->iniSol));
    Sol_cpy(current,vns->iniSol);
    VND *vnd = VND_create(vns->inst, current, argv, argc);

    while(_time < timeRem) {

        while( it <= nMoves && _time < timeRem) {

            valid = Neighbor_callStocChosen(neighborhood,vns->iniSol,vns->bestSol, current, it, test);
            if(!valid) {
                it++;
                continue;
            }

            Sol_rebuild_opt(current, vns->iniSol);

            _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            Neighbor_setNullLastJobModify(neighborhood);

            //VND_runDet(vnd, neighborhood, timeRem - _time, 6, firstImprovement, vns->nChangesModes, vns->nTimesJobOnModes,test);
            VND_runDet(vnd, neighborhood, timeRem - _time, nMoves, firstImprovement, test);

            if(Sol_getCost(current) < Sol_getCost(vns->iniSol) ) {
                Sol_cpy(vns->iniSol,current);
                if(Sol_getCost(current) < Sol_getCost(vns->bestSol)) {
#ifdef DEBUG
                    _time = ( (double)clock()/CLOCKS_PER_SEC );
                    printf("\nVNS %ld %f\n ", Sol_getCost(current), _time);
                    fflush(stdout);
                    fflush(stdin);
#endif
                    Sol_cpy( vns->bestSol,current);
                    it = 1;
                } else
                    it++;
            } else
                it++;

            _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
        }

        it = 1;
        _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

    }
    VND_free(&vnd);

}

/*no testing*/
void VNS_run_reduced(VNS *vns, Neighborhood* neighborhood, double timeRem, int nNeighbor, int firstImprovement, Test *test, char **argv, int argc )
{
    assert( vns != NULL);
    assert( neighborhood != NULL);

    double _time = 0;
    clock_t tStart = clock();

    int nMoves = VNS_getNMoves(vns);
    int it = 1, valid = 0;

    _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

    Solution *current = Sol_create(Sol_inst( vns->iniSol));
    Sol_cpy(current,vns->iniSol);

    while(_time < timeRem) {

        while( it <= nMoves && _time < timeRem) {

            valid = Neighbor_callStocChosen(neighborhood,vns->iniSol,vns->bestSol, current, it, test);
            if(!valid) {
                it++;
                continue;
            }

            Sol_rebuild_opt(current, vns->iniSol);

            _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
            Neighbor_setNullLastJobModify(neighborhood);

            if(Sol_getCost(current) < Sol_getCost(vns->iniSol) ) {
                Sol_cpy(vns->iniSol,current);
                if(Sol_getCost(current) < Sol_getCost(vns->bestSol)) {
#ifdef DEBUG
                    _time = ( (double)clock()/CLOCKS_PER_SEC );
                    printf("\nVNS %ld %f\n ", Sol_getCost(current), _time);
                    fflush(stdout);
                    fflush(stdin);
#endif
                    Sol_cpy( vns->bestSol,current);
                    it = 1;
                } else
                    it++;
            } else
                it++;

            _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
        }

        it = 1;
        _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

    }

}

void VNS_rna(VNS *vns,  Solution *current, Neighborhood* neighborhood, double timeRem, Test *test)
{
    int cont = 0, valid =0;
    double _time = 0;
    clock_t tStart = clock();

    Solution *newcurrent = Sol_create(Sol_inst( current));
    Sol_cpy(newcurrent,current);

    _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

    while(cont<vns->itRNA && _time < timeRem) {
        valid = Neighbor_callStocLS( neighborhood, current, current, newcurrent, vns->la, 0);
        if(!valid) {
            cont++;
            continue;
        }

        Sol_rebuild_opt(newcurrent,current);
        if(Sol_getCost(newcurrent) < Sol_getCost(current) ) {
            Sol_cpy(current,newcurrent);
#ifdef DEBUG
            _time = ( (double)clock()/CLOCKS_PER_SEC );
            printf("\nRNA %ld %f ", Sol_getCost(current), _time);
            fflush(stdout);
            fflush(stdin);
#endif
            cont=0;

        } else {
            Sol_cpy(newcurrent,current);
            cont++;
        }
        Neighbor_setNullLastJobModify(neighborhood);
        _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
    }

    Sol_free(&newcurrent);

}

void VNS_run_RNA(VNS *vns, Neighborhood* neighborhood, double timeRem, int nNeighbor,  Test *test, char **argv, int argc )
{

    assert( vns != NULL);
    assert( neighborhood != NULL);

    double _time = 0;
    clock_t tStart = clock();

    int nMoves = VNS_getNMoves(vns);
    int it = 1, valid = 0;
    _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

    while(_time < timeRem) {

        Solution *current = Sol_create(Sol_inst( vns->iniSol));
        Sol_cpy(current, vns->iniSol);

        while( it <= nMoves && _time < timeRem) {

            valid = Neighbor_callStocChosen(neighborhood,vns->iniSol,vns->bestSol, current, it, test);
            if(!valid) {
                it++;
                continue;
            }
            Sol_rebuild_opt(current, vns->iniSol);

            _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

            Neighbor_setNullLastJobModify(neighborhood);

            VNS_rna(vns, current, neighborhood, timeRem-_time, test);

            if(Sol_getCost(current) < Sol_getCost(vns->iniSol) ) {
                Sol_cpy(vns->iniSol,current);
                if(Sol_getCost(current) < Sol_getCost(vns->bestSol)) {
#ifdef DEBUG
                    _time = ( (double)clock()/CLOCKS_PER_SEC );
                    printf("\n---VNS %ld %f\n ", Sol_getCost(current), _time);
                    fflush(stdout);
                    fflush(stdin);
#endif
                    Sol_cpy( vns->bestSol,current);
                }
                it = 1;
            } else {
#ifdef DEBUG
                printf("\n---VNS N %ld %f\n ", Sol_getCost(current), _time);
                fflush(stdout);
                fflush(stdin);
#endif
                Sol_cpy(current,vns->iniSol);
                it++;
            }
            _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
        }
        it = 1;

        _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

        Sol_free(&current);
    }

}

void VNS_run_RNA_shake(VNS *vns, Neighborhood* neighborhood, double timeRem, int nNeighbor,  Test *test, char **argv, int argc )
{

    assert( vns != NULL);
    assert( neighborhood != NULL);

    double _time = 0;
    clock_t tStart = clock();

    int nMoves = VNS_getNMoves(vns);
    _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

    while(_time < timeRem) {

        Solution *current = Sol_create(Sol_inst( vns->iniSol));
        Sol_cpy(current, vns->iniSol);


        Neighbor_Shake(vns, neighborhood, vns->iniSol, vns->bestSol, current, nMoves, test);

#ifdef DEBUG
        //    printf("\n---SHAKE FINAL %ld \n ", Sol_getCost(current));
#endif

        _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );


        VNS_rna(vns, current, neighborhood, timeRem-_time, test);
#ifdef DEBUG
        //  printf("\n---RNA FINAL %ld \n ", Sol_getCost(current));
#endif

        if(Sol_getCost(current) < Sol_getCost(vns->iniSol) ) {
            Sol_cpy(vns->iniSol,current);
            if(Sol_getCost(current) < Sol_getCost(vns->bestSol)) {
#ifdef DEBUG
                _time = ( (double)clock()/CLOCKS_PER_SEC );
                printf("\n---VNS %ld %f\n ", Sol_getCost(current), _time);
                fflush(stdout);
                fflush(stdin);
#endif
                Sol_cpy( vns->bestSol,current);
            }
        } else {
            //#ifdef DEBUG
            // printf("\n---VNS N %ld %f\n ", Sol_getCost(current), _time);
            // fflush(stdout);
            // fflush(stdin);
            //#endif
            Sol_cpy(current,vns->iniSol);
        }


        _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

        Sol_free(&current);
    }

}

void VNS_run_LAHC_smartshake(VNS *vns, Neighborhood* neighborhood, double timeRem, int nNeighbor, char **argv, int argc )
{

    assert( vns != NULL);
    assert( neighborhood != NULL);

    double _time = 0;
    clock_t tStart = clock();

    int nMovesShake = VNS_getNMoves(vns);
    int sizeSamplingShake = VNS_getNSizeSamplingShake(vns);
    _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

    Solution *current = Sol_create(Sol_inst( vns->iniSol));
    Sol_cpy(current, vns->iniSol);

    while(_time < timeRem) {

        _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
        LAHC *lahc = LAHC_create(Sol_inst( vns->iniSol), current, neighborhood, argv, argc);
        LAHC_run_vns(lahc, vns, neighborhood, timeRem-_time, vns->itLAHC, argv[2]);

#ifdef DEBUG
        printf("\n---LAHC FINAL %ld \n ", Sol_getCost(current));
#endif

        if(Sol_getCost(current) < Sol_getCost(vns->iniSol) ) {
            Sol_cpy(vns->iniSol,current);
            if(Sol_getCost(current) < Sol_getCost(vns->bestSol)) {
#ifdef DEBUG
                _time = ( (double)clock()/CLOCKS_PER_SEC );
                printf("\n---VNS %ld %f\n ", Sol_getCost(current), _time);
                fflush(stdout);
                fflush(stdin);
#endif
                Sol_cpy( vns->bestSol, current);
            }
        } else
            Sol_cpy(current, vns->iniSol);

        _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

        Neighbor_SmartShake(vns, neighborhood, current, vns->bestSol, nMovesShake, sizeSamplingShake, timeRem-_time);

#ifdef DEBUG
        printf("\n---SHAKE FINAL %ld \n ", Sol_getCost(current));
#endif

        _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
    }

    Sol_free(&current);
}

void VNS_run_LAHC_shake2(VNS *vns, Neighborhood* neighborhood, double timeRem, int nNeighbor,  Test *test, char **argv, int argc )
{

    assert( vns != NULL);
    assert( neighborhood != NULL);

    double _time = 0;
    clock_t tStart = clock();

    int nMoves = VNS_getNMoves(vns);
    int it = 1;
    _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

    while(_time < timeRem) {

        Solution *current = Sol_create(Sol_inst( vns->iniSol));
        Sol_cpy(current, vns->iniSol);

        while( it <= 14 && _time < timeRem) {

            Neighbor_Shake2(vns, neighborhood, vns->iniSol, vns->bestSol, current, it, nMoves, test);

#ifdef DEBUG
            //    printf("\n---SHAKE FINAL %ld \n ", Sol_getCost(current));
#endif


            LAHC *lahc = LAHC_create(Sol_inst( vns->iniSol), current, neighborhood, argv, argc);

            _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

            LAHC_run_nIt(lahc, neighborhood, timeRem-_time, vns->itLAHC, test, argv[2]);

#ifdef DEBUG
            //  printf("\n---RNA FINAL %ld \n ", Sol_getCost(current));
#endif

            if(Sol_getCost(current) < Sol_getCost(vns->iniSol) ) {
                Sol_cpy(vns->iniSol,current);
                if(Sol_getCost(current) < Sol_getCost(vns->bestSol)) {
#ifdef DEBUG
                    _time = ( (double)clock()/CLOCKS_PER_SEC );
                    printf("\n---VNS %d %ld %f\n ", Neighbor_getLastNeigh(neighborhood), Sol_getCost(current), _time);
                    fflush(stdout);
                    fflush(stdin);
#endif
                    Sol_cpy( vns->bestSol,current);
                }
                it = 1;
            } else {
                //#ifdef DEBUG
                // printf("\n---VNS N %ld %f\n ", Sol_getCost(current), _time);
                // fflush(stdout);
                // fflush(stdin);
                //#endif
                Sol_cpy(current,vns->iniSol);
                it++;
            }
            _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
        }
        it = 1;

        _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

        Sol_free(&current);
    }

}

/*no testing*/
void VNS_run_RNA_shake2(VNS *vns, Neighborhood* neighborhood, double timeRem, int nNeighbor,  Test *test, char **argv, int argc )
{

    assert( vns != NULL);
    assert( neighborhood != NULL);

    double _time = 0;
    clock_t tStart = clock();

    int nMoves = VNS_getNMoves(vns);
    int it = 1;
    _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

    while(_time < timeRem) {

        Solution *current = Sol_create(Sol_inst( vns->iniSol));
        Sol_cpy(current, vns->iniSol);

        while( it <= 14 && _time < timeRem) {

            Neighbor_Shake2(vns, neighborhood, vns->iniSol, vns->bestSol, current, it, nMoves, test);

#ifdef DEBUG
            //    printf("\n---SHAKE FINAL %ld \n ", Sol_getCost(current));
#endif

            _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );


            VNS_rna(vns, current, neighborhood, timeRem-_time, test);
#ifdef DEBUG
            //  printf("\n---RNA FINAL %ld \n ", Sol_getCost(current));
#endif

            if(Sol_getCost(current) < Sol_getCost(vns->iniSol) ) {
                Sol_cpy(vns->iniSol,current);
                if(Sol_getCost(current) < Sol_getCost(vns->bestSol)) {
#ifdef DEBUG
                    _time = ( (double)clock()/CLOCKS_PER_SEC );
                    printf("\n---VNS %ld %f\n ", Sol_getCost(current), _time);
                    fflush(stdout);
                    fflush(stdin);
#endif
                    Sol_cpy( vns->bestSol,current);
                }
                it = 1;
            } else {
                //#ifdef DEBUG
                // printf("\n---VNS N %ld %f\n ", Sol_getCost(current), _time);
                // fflush(stdout);
                // fflush(stdin);
                //#endif
                Sol_cpy(current,vns->iniSol);
                it++;
            }
            _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
        }
        it = 1;

        _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

        Sol_free(&current);
    }

}

Cost VNS_penaltyModes(VNS* vns, Neighborhood* neighborhood, Cost foCurrent)
{

    assert(vns != NULL);
    assert(neighborhood != NULL);

    Cost fator = 0, p = 0;
    fator = vns->perc*foCurrent;

    int move = Neighbor_getLastNeigh(neighborhood), k = -1;

    if( move == changeOneMode) k = 1;
    else if( move == changeTwoMode) k = 2;
    else if( move == changeThreeMode) k = 3;
    else if( move == changeFourMode) k = 4;


    if(k>0) {

        int min=INT_MAX, max=0, sum = 0, it=0;
        double avg=0;

        for(int j = 0; j < Inst_nJobs(vns->inst) ; j++) {
            for(int m = 0 ; m < 3 ; m++) {
                max = MAX(vns->nTimesJobOnModes[j][m], max);
                min = MIN(vns->nTimesJobOnModes[j][m], min);
                sum += vns->nTimesJobOnModes[j][m];
                it++;
            }
        }

        avg = (double) sum/(double)it;

        for(int i = 0; i < k; i++) {

            int idxJob = Neighbor_getLastJob(neighborhood,i);
            int idxMode =Neighbor_getNewModes(neighborhood,i);

            p += ((vns->nTimesJobOnModes[idxJob][idxMode] - avg) / avg)* fator;
#ifdef DEBUG
            //  printf("Job %d Mode %d Res %d \n Avg %f Fator %ld\n", idxJob, idxMode, vns->nTimesJobOnModes[idxJob][idxMode], avg, fator);
#endif
        }
    }

#ifdef DEBUG
    if(p!=0)
        printf("Penalty Residency Job on Modes %ld\n", p);
#endif
    return p;

}

Cost VNS_penaltyJobs(VNS* vns, Neighborhood* neighborhood, Cost foCurrent)
{

    assert(vns != NULL);
    assert(neighborhood != NULL);

    Cost fator = 0, p = 0;
    fator = vns->percRS*foCurrent;

    int min=INT_MAX, max=0, sum = 0, it=0;
    double avg=0;

    for(int j = 0; j < Inst_nJobs(vns->inst) ; j++) {
        for(int s = 0 ; s < Inst_nJobs(vns->inst) ; s++) {
            max = MAX(vns->nTimesJobOnSequence[j][s], max);
            min = MIN(vns->nTimesJobOnSequence[j][s], min);
            sum += vns->nTimesJobOnSequence[j][s];
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
        p += ((vns->nTimesJobOnSequence[idxJob][idxSequence] - avg) / avg)* fator;
    }

#ifdef DEBUG
    if(p!=0)
        printf("Penalty Residency Job on Sequence %ld\n", p);
#endif

    return p;

}

void VNS_increasingTransitivityOfModes(VNS* vns, Neighborhood* neighborhood)
{

    assert(vns != NULL);
    assert(neighborhood != NULL);

    if( Neighbor_getLastNeigh(neighborhood) == changeOneMode)
        vns->nChangesModes[Neighbor_getLastJob(neighborhood, 0)]++;
    else if ( Neighbor_getLastNeigh(neighborhood) == changeTwoMode) {
        vns->nChangesModes[Neighbor_getLastJob(neighborhood, 0)]++;
        vns->nChangesModes[Neighbor_getLastJob(neighborhood, 1)]++;
    } else if ( Neighbor_getLastNeigh(neighborhood) == changeThreeMode) {
        vns->nChangesModes[Neighbor_getLastJob(neighborhood, 0)]++;
        vns->nChangesModes[Neighbor_getLastJob(neighborhood, 1)]++;
        vns->nChangesModes[Neighbor_getLastJob(neighborhood, 2)]++;
    } else if ( Neighbor_getLastNeigh(neighborhood) == changeFourMode) {
        vns->nChangesModes[Neighbor_getLastJob(neighborhood, 0)]++;
        vns->nChangesModes[Neighbor_getLastJob(neighborhood, 1)]++;
        vns->nChangesModes[Neighbor_getLastJob(neighborhood, 2)]++;
        vns->nChangesModes[Neighbor_getLastJob(neighborhood, 3)]++;
    }
}

Cost VNS_penaltyTransModes(VNS* vns, Neighborhood* neighborhood, Cost foCurrent)
{

    assert(vns != NULL);
    assert(neighborhood != NULL);

    Cost fator = 0, p = 0;

    fator = vns->perc*foCurrent;

    int move = Neighbor_getLastNeigh(neighborhood), k = -1;

    if( move == changeOneMode) k = 1;
    else if( move == changeTwoMode) k = 2;
    else if( move == changeThreeMode) k = 3;
    else if( move == changeFourMode) k = 4;


    if(k>0) {

        int min=INT_MAX, max=0, sum = 0, it=0;
        double avg=0;

        for(int j = 0; j < Inst_nJobs(vns->inst) ; j++) {
            max = MAX(vns->nChangesModes[j], max);
            min = MIN(vns->nChangesModes[j], min);
            sum += vns->nChangesModes[j];
            it++;
        }

        avg = (double) sum/(double)it;

        for(int i = 0; i < k; i++) {
            int idxJob = Neighbor_getLastJob(neighborhood, i);
#ifdef DEBUG
            //  printf("Job %d Trans %d \n Avg %f Factor %ld \n", idxJob, vns->nChangesModes[idxJob], avg, fator);
#endif
            p += ((vns->nChangesModes[idxJob]- avg) / avg)*fator;
        }
    }

#ifdef DEBUG
    if(p!=0)
        printf("Penalty %ld\n", p);
#endif

    return p;

}

void VNS_increasingTransitivityOfSequence(VNS* vns, Neighborhood* neighborhood)
{

    assert(vns != NULL);
    assert(neighborhood != NULL);

    for(int i = 0 ; i < Neighbor_getContLastJ(neighborhood) ; i++)
        vns->nChangesSequence[Neighbor_getLastJobModify(neighborhood,i)]++;

}

//transSequence
Cost VNS_penaltyTransJobs(VNS* vns, Neighborhood* neighborhood, Cost foCurrent)
{

    assert(vns != NULL);
    assert(neighborhood != NULL);

    Cost fator = 0, p = 0;
    fator = vns->perc*foCurrent;

    int min=INT_MAX, max=0, sum = 0, it=0;
    double avg=0;

    for(int j = 0; j < Inst_nJobs(vns->inst) ; j++) {
        max = MAX(vns->nChangesSequence[j], max);
        min = MIN(vns->nChangesSequence[j], min);
        sum += vns->nChangesSequence[j];
        it++;
    }

    avg = (double) sum/(double)it;

    for(int j = 0; j < Neighbor_getContLastJ(neighborhood) ; j++) {
        int idxJob = Neighbor_getLastJobModify(neighborhood,j);
#ifdef DEBUG
        // printf("Job %d Trans %d \n Avg %f Factor %ld \n", idxJob, vns->nChangesSequence[idxJob], avg, fator);
#endif
        p += ((vns->nChangesSequence[idxJob]- avg) / avg)* fator;
    }

#ifdef DEBUG
    if(p!=0)
        printf("Penalty Transitivity Sequence %ld\n", p);
#endif
    return p;

}

void VNS_increasingResidencyJobInMode(VNS* vns, Solution* current)
{

    assert(vns != NULL);
    assert(current != NULL);

    const Job* job;
    int m;
    for(int i = 0 ; i < Inst_nJobs(vns->inst) ; i++) {
        job = Inst_job(vns->inst, i);
        m = Sol_getMode(current,Job_index(job));
        vns->nTimesJobOnModes[Job_index(job)][m]++;
    }
}

void VNS_increasingResidencyJobInSequence(VNS* vns, Solution* current)
{

    assert(vns != NULL);
    assert(current != NULL);

    const Job* job;
    int s =0;
    for(int i = 0 ; i < Inst_nJobs(vns->inst) ; i++) {
        job = Inst_job(vns->inst, i);
        s = Sol_getPosJob(current,Job_index(job));
        vns->nTimesJobOnSequence[Job_index(job)][s]++;
    }

}

int VNS_getNMoves(VNS *vns)
{
    assert( vns != NULL );

    return vns->nMoves;

}

int VNS_getNSizeSamplingShake(VNS *vns)
{
    assert( vns != NULL );

    return vns->sizeSamplingShake;

}

void VNS_free( VNS **_vns )
{

    VNS *vns = *_vns;

    for(int i =0 ; i < Inst_nJobs( Sol_inst(vns->iniSol) ); i++)
        free( vns->nTimesJobOnModes[i] );
    free( vns->nTimesJobOnModes);


    for(int i = 0; i < Inst_nJobs(Sol_inst(vns->iniSol)) ; i++)
        free(vns->nTimesJobOnSequence[i]);
    free( vns->nTimesJobOnSequence );


    Sol_free( &vns->iniSol );
    Sol_free( &vns->bestSol );
    LA_free( &vns->la );

    free( vns->nChangesModes );
    free( vns->nChangesSequence );

    free( vns );

    *_vns = NULL;

}












/*
void VNS_run_allNeigh(VNS *vns, double timeRem, int nNeighbor, int firstImprovement, Test *test, char **argv, int argc )
{

    assert( vns != NULL);

    double _time = 0;
    clock_t tStart = clock();

    int kN = 0;

    Solution *current = Sol_create(Sol_inst( vns->iniSol));
    _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

    while(_time < timeRem) {
        kN=0;
        while(kN < nNeighbor) {
            Sol_cpy(current, vns->iniSol);

            Test_setCurrentFO(test, Sol_getCost(current));

            Neighbor_callStocChosen(vns->neighborhood, vns->iniSol,vns->bestSol,current,kN, vns->nChangesModes, test);
            Sol_rebuild_opt(current, vns->iniSol);

            Test_callTest(test, Test_getCurrentNeigh(test), Test_getCurrentTime(test), Sol_getCost(current));

            _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

            VND *vnd = VND_create(vns->inst, current, argv, argc);
            VND_runDet(vnd, vns->neighborhood, timeRem - _time, 6, firstImprovement, vns->nChangesModes, test);

            if(Sol_getCost(current) < Sol_getCost(vns->iniSol) ) {
                Sol_cpy(vns->iniSol,current);
                kN=0;
                if(Sol_getCost(current) < Sol_getCost(vns->bestSol)) {
                    //_time = ( (double)clock()/CLOCKS_PER_SEC );
                    //printf("\n%ld %f\n ", Sol_getCost(current) , _time);
                    Sol_cpy( vns->bestSol,current);
                }
            } else {
                Sol_cpy(current,vns->iniSol);
                kN++;
            }

        }
        _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

    }

    Sol_free(&current);

}
*/

/*
int VNS_getNChanges(VNS *vns, int j)
{

    assert( vns != NULL );

    return vns->nChangesModes[j];

}

void VNS_setNChanges(VNS *vns, int j)
{

    assert( vns != NULL );

    vns->nChangesModes[j]++;

}

int VNS_getAllNeighbor(VNS *vns)
{

    assert( vns != NULL );

    return vns->allNeighbor;

}
*/

