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
#include "sa.h"

struct _SA {

    Solution *iniSol;
    Solution *bestSol;
    int SAmax;
    double alpha;
    int T;
    int *nChangesModes;
    int **nTimesJobOnModes;
    int learning;

    LearningAutomata *la;

    const struct _Instance *inst;

};


void SA_checkArgs(SA *sa, char **argv, int argc)
{

    for(int n=0 ; n < argc; n++) {
        if (strcmp(argv[n],"-T") == 0) {
            n++;
            sa->T =  atoi(argv[n]);
            continue;
        }

        if (strcmp(argv[n],"-SAmax") == 0) {
            n++;
            sa->SAmax =  atoi(argv[n]);
            continue;
        }

        if (strcmp(argv[n],"-alpha") == 0) {

            n++;
            sa->alpha = (double) atof(argv[n]);
            continue;
        }

        if (strcmp(argv[n],"-learning") == 0) {
            n++;
            sa->learning = atoi(argv[n]);
            continue;
        }
    }
}


SA *SA_create( const Instance *inst, Neighborhood* neighborhood, Solution* sol, char **argv, int argc )
{
    SA* sa;

    ALLOCATE_INI( sa, SA );

    sa->iniSol = sol;
    sa->bestSol = Sol_create(inst);
    Sol_cpy( sa->bestSol, sa->iniSol );

    sa->T = 10000;
    sa->SAmax = 100;
    sa->alpha = 0.6;
    sa->learning = 0;

    ALLOCATE_VECTOR_INI( sa->nChangesModes, int, Inst_nJobs( inst ) );

    ALLOCATE_VECTOR( sa->nTimesJobOnModes, int*, Inst_nJobs( inst ) );
    for(int i = 0 ; i < Inst_nJobs( inst ) ; i++)
        ALLOCATE_VECTOR_INI( sa->nTimesJobOnModes[i], int, 3);

    sa->inst = inst;

    sa->la = LA_create( Neighbor_nNeighborhood(neighborhood), Neighbor_getIntensity(neighborhood), argv, argc);

    SA_checkArgs(sa, argv, argc);

    return sa;
}

void SA_increasingResidence(SA* sa, Solution* current)
{

    assert(sa != NULL);
    assert(current != NULL);

    const Job* job;
    int m;
    for(int i = 0 ; i < Inst_nJobs(sa->inst) ; i++) {
        job = Inst_job(sa->inst, i);
        m = Sol_getMode(current,Job_index(job));
        sa->nTimesJobOnModes[Job_index(job)][m]++;
    }
}

void SA_increasingTransitivity(SA* sa, Neighborhood* neighborhood)
{

    assert(sa != NULL);
    assert(neighborhood != NULL);

    if( Neighbor_getLastNeigh(neighborhood) == changeOneMode)
        sa->nChangesModes[Neighbor_getLastJob(neighborhood, 0)]++;
    else if ( Neighbor_getLastNeigh(neighborhood) == changeTwoMode) {
        sa->nChangesModes[Neighbor_getLastJob(neighborhood, 0)]++;
        sa->nChangesModes[Neighbor_getLastJob(neighborhood, 1)]++;
    } else if ( Neighbor_getLastNeigh(neighborhood) == changeThreeMode) {
        sa->nChangesModes[Neighbor_getLastJob(neighborhood, 0)]++;
        sa->nChangesModes[Neighbor_getLastJob(neighborhood, 1)]++;
        sa->nChangesModes[Neighbor_getLastJob(neighborhood, 2)]++;
    } else if ( Neighbor_getLastNeigh(neighborhood) == changeFourMode) {
        sa->nChangesModes[Neighbor_getLastJob(neighborhood, 0)]++;
        sa->nChangesModes[Neighbor_getLastJob(neighborhood, 1)]++;
        sa->nChangesModes[Neighbor_getLastJob(neighborhood, 2)]++;
        sa->nChangesModes[Neighbor_getLastJob(neighborhood, 3)]++;
    }
}


void SA_run(SA *sa, Neighborhood* neighborhood, double timeRem, Test *test)
{

    double _time = 0;
    clock_t tStart = clock();

    Solution *current = Sol_create(sa->inst);


    int delta;
    double xrand;
    int T, iterT=0, valid = 0;
    int nUnimproved=0;
    T = sa->T;

    while(_time < timeRem ) {

        while( iterT < sa->SAmax &&  _time < timeRem ) {
            iterT++;

            Sol_cpy(current, sa->iniSol);

            Test_setCurrentFO(test, Sol_getCost( sa->bestSol ));

            /*increasing the permanence of a job in a mode*/
            SA_increasingResidence(sa, current);

            valid = Neighbor_callStocLS( neighborhood, sa->iniSol, sa->bestSol, current, sa->la, sa->learning, test);
            if(valid) Sol_rebuild_opt(current, sa->iniSol);

            delta = Sol_getCost(current) - Sol_getCost(sa->iniSol);
            if( delta < 0) {
                if(sa->learning)
                    LA_update( sa->la, 1 );

                /*increasing the transitivity of modes*/
                SA_increasingTransitivity(sa, neighborhood);

                Sol_cpy(sa->iniSol, current );
                if( Sol_getCost(current) < Sol_getCost(sa->bestSol)) {
#ifdef DEBUG
                    _time = ( (double)clock()/CLOCKS_PER_SEC );
                    printf("\n%ld %f", Sol_getCost(sa->bestSol), _time);
#endif // DEBUG
                    Sol_cpy( sa->bestSol, sa->iniSol );
                }
            } else {
                xrand = ((double)rand()) / (((double)RAND_MAX)+1);
                if( (xrand < pow(E,-(delta/T)))) {
                    if(sa->learning)
                        LA_update( sa->la, 1 );
                    Sol_cpy( sa->iniSol, current );
                }
                nUnimproved++;
                if(sa->learning) {
                    LA_update( sa->la, 0 );
                    if(nUnimproved >= LA_getResetInterval(sa->la)) {
                        //  printf("\nReset Interval %d, %d\n", nUnimproved, LA_getResetInterval(sa->la));
                        LA_reset( sa->la );
                        nUnimproved=0;
                    }
                }
                Sol_cpy( current, sa->iniSol );
            }
            _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );

            Neighbor_setNullLastJobModify(neighborhood);
            Test_callTest(test, Test_getCurrentNeigh(test), Test_getCurrentTime(test), Sol_getCost(sa->bestSol));
        }

        T *= sa->alpha;
        iterT = 0;
        if(T <= 0) T = sa->T;
        _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );



    }

    Sol_cpy( sa->iniSol, sa->bestSol );
    Sol_free(&current);

}

double SA_getAlpha(SA *_sa)
{
    assert( _sa != NULL );

    return _sa->alpha;
}

int SA_getT(SA *_sa)
{
    assert( _sa != NULL );

    return _sa->T;
}

int SA_getSAmax(SA *_sa)
{
    assert( _sa != NULL );

    return _sa->SAmax;
}


void SA_free( SA **_sa )
{

    SA *sa = *_sa;

    for(int i = 0 ; i < Inst_nJobs( Sol_inst(sa->bestSol) ) ; i++)
        free( sa->nTimesJobOnModes[i] );
    free( sa->nTimesJobOnModes);

    free( sa->nChangesModes );

    Sol_free( &sa->iniSol );
    Sol_free( &sa->bestSol );

    LA_free( &sa->la);

    free( sa );

    *_sa = NULL;

}
