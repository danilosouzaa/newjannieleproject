/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                              of Santos, H.G., Toffolo, T.A.M. and Baltar, D.
 */

#include "learning.h"
#include "macros.h"
#include <string.h>
#include "neighborhood.h"
#include <assert.h>
#include "list_int.h"

//#define PRINT_LEARNING

const double EPSILON = 0; //0
const double LEARNING_RATE = 1e-6; //4,3,2...
const int N_REWARDED_ACTIONS = 1;

struct _LearningAutomata {
    /* Probability distribution (one for each neighborhood) */
    double *initialProbabilities, *probabilities;
    /* Number of possible actions (neighborhoods) */
    int nProbabilities;
    /* Index of last selected action */
    int lastSelection;
    /* Total number of iterations */
    int iters;
    /* Linked list with the last N_REWARDED_ACTIONS or rewardedActions selected actions */
    ListInt *lastNSelections;
    /* Parameters to specify the learning rate(s) */
    double learningRate, learningRate2;

    int timePrint;
    /*To reset learning*/
    int  resetInterval;
    int rewardedActions;
    int sideway;


};


void LA_checkArgs(LearningAutomata *la, char **argv, int argc)
{

    assert(la!=NULL);

    for(int n=0 ; n < argc; n++) {
        if (strcmp(argv[n],"-learningRate") == 0) {
            n++;
            la->learningRate =  atof(argv[n]);
            la->learningRate2 = la->learningRate * EPSILON;
            continue;
        }

        if (strcmp(argv[n],"-resetInterval") == 0) {
            n++;
            la->resetInterval =  atoi(argv[n]);
            continue;
        }

        if (strcmp(argv[n],"-rewardedActions") == 0) {
            n++;
            la->rewardedActions =  atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-sw") == 0) {
            n++;
            la->sideway = atoi( argv[n] );
            continue;
        }
    }
}

/* creates and returns an empty LearningAutomata */
LearningAutomata *LA_create( const int nIntensities, const long double intensities[], char **argv, int argc )
{
    LearningAutomata *la;
    ALLOCATE( la, LearningAutomata );

    ALLOCATE_VECTOR( la->initialProbabilities, double, nIntensities );
    ALLOCATE_VECTOR( la->probabilities, double, nIntensities );

    la->nProbabilities = nIntensities;
    la->lastSelection = 0;
    la->iters = 0;
    la->lastNSelections = LInt_create();
    la->learningRate = LEARNING_RATE; //1e-3, 1e-4, 1e-5
    la->learningRate2 = LEARNING_RATE * EPSILON;

    double sumWeights = 0;
    for ( int i = 0; i < la->nProbabilities; ++i )
        sumWeights += intensities[i];

    for ( int i = 0; i < la->nProbabilities; ++i ) {
        la->probabilities[i] = intensities[i] / sumWeights;
        la->initialProbabilities[i] = intensities[i] / sumWeights;
    }

    la->resetInterval = INT_MAX_M; //1000, 10000
    la->rewardedActions = N_REWARDED_ACTIONS; //1,2
    la->timePrint = 0;
    LA_checkArgs(la,argv, argc);

    //procurar sobre reset no texto do tony

    return la;
}

/* returns the next neighborhood index based on the what has been learned */
int LA_next( LearningAutomata *la )
{
    assert(la!=NULL);
    // return action with highest probability
    double randomProb = DOUBLE_RANDOM( 0.0, 1.0 );
    double cumProb = 0.0; // cumulative probability

    for ( int i = 0; i < la->nProbabilities; ++i ) {
        if ( cumProb <= randomProb && randomProb <= cumProb + la->probabilities[i] ) {
            la->lastSelection = i;
            break;
        }
        cumProb = cumProb + la->probabilities[i];
    }

    la->iters++;

    LInt_push( la->lastNSelections, la->lastSelection );
    if ( LInt_size( la->lastNSelections ) > la->rewardedActions ) // la->rewardedAction or N_REWARDED_ACTIONS
        LInt_poll( la->lastNSelections );

    return la->lastSelection;
}

/* resets the learning automata to its original state */
void LA_reset( LearningAutomata *la )
{
    assert(la!=NULL);
    la->iters = 0;
    for ( int i = 0; i < la->nProbabilities; ++i )
        la->probabilities[i] = la->initialProbabilities[i];

    //LInt_clear( la->lastNSelections );
    LInt_free(la->lastNSelections);
    la->lastNSelections = LInt_create();

    la->lastSelection = 0;

}

/* updates a single action (neighborhood) */
/* this should be a 'private' method, i.e. invisible to others */
void LA_updateProbability( LearningAutomata *la, int usedAction, double obtainedReinforcement )
{
    assert(la!=NULL);

    for ( int i = 0; i < la->nProbabilities; ++i ) {
        if ( usedAction == i ) {
            la->probabilities[i] = la->probabilities[i] + la->learningRate * obtainedReinforcement * ( 1 - la->probabilities[i] )
                                   - la->learningRate2 * ( 1 - obtainedReinforcement ) * la->probabilities[i];
        } else {
            la->probabilities[i] = la->probabilities[i] - la->learningRate * obtainedReinforcement * la->probabilities[i]
                                   + la->learningRate2 * ( 1 - obtainedReinforcement )
                                   * ( ( ( double ) 1 / ( ( double ) la->nProbabilities - ( double ) 1 ) ) - la->probabilities[i] );
        }
    }
}

/* updates (gives feedback to) the learning automata */
/* the obtainedReinforcement x should be a value such that 0.0 <= x <= 1.0 */
void LA_update( LearningAutomata *la, const double obtainedReinforcement )
{
    assert(la!=NULL);

    int i = LInt_size( la->lastNSelections ) - 1;
    ListIntIter *iterator = LInt_iterator( la->lastNSelections );
    while ( LIntIter_hasNext( iterator ) ) {
        int usedAction = LIntIter_next( &iterator );
        LA_updateProbability( la, usedAction, obtainedReinforcement / pow( 2, i ) ); // TODO: check discount function
        i--;
    }
}

/* returns the number of iterations of the learning automata */
int LA_getIters( const LearningAutomata *la )
{
    assert(la!=NULL);

    return la->iters;
}

/* returns the number of reset interval of the learning automata */
int LA_getResetInterval( const LearningAutomata *la )
{
    assert(la!=NULL);
    return la->resetInterval;
}

/* returns the last selection (choice) of the learning automata */
int LA_getLastSelection( const LearningAutomata *la )
{
    assert(la!=NULL);
    return la->lastSelection;
}

/* returns the current learning rate */
double LA_getLearningRate( const LearningAutomata *la )
{
    assert(la!=NULL);
    return la->learningRate;
}

/* sets the probabilities of the learning using the intensities value */
void LA_setIntensities( LearningAutomata *la, double intensities[] )
{
    assert(la!=NULL);
    double sumWeights = 0;
    for ( int i = 0; i < la->nProbabilities; ++i )
        sumWeights += intensities[i];

    for ( int i = 0; i < la->nProbabilities; ++i ) {
        la->probabilities[i] = intensities[i] / sumWeights;
        la->initialProbabilities[i] = intensities[i] / sumWeights;
    }
}

/* sets the learning rate */
void LA_setLearningRate( LearningAutomata *la, double learningRate )
{
    assert(la!=NULL);
    la->learningRate = learningRate;
}

/* sets the probabilities of the learning automata */
void LA_setProbabilities( LearningAutomata *la, double probabilities[] )
{
    assert(la!=NULL);
    for ( int i = 0; i < la->nProbabilities; ++i )
        la->probabilities[i] = probabilities[i];
}

/* sets the probabilities of the learning automata */
void LA_setTime( LearningAutomata *la)
{
    assert(la!=NULL);
    la->timePrint += 10;
}


void LA_printProbabilities( LearningAutomata *la, float pesw)
{
    char fil[30] = "probabilities_";
    char bufferII[60];
    char side[60];
    char psw[60];
    sprintf(bufferII, "%f", la->learningRate);
    strcat (fil,bufferII);
    strcat (fil,"_");
    sprintf(side, "%d", la->sideway);
    strcat (fil,side);
    strcat (fil,"_");
    sprintf(psw, "%f", pesw);
    strcat (fil,psw);
    strcat (fil,".txt");


    FILE* ap = fopen(fil, "a");

    for ( int i = 0; i < la->nProbabilities; ++i ) {
        char* sigla= "";
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

        //printf("\n%d %s %.8f", la->timePrint, sigla, la->probabilities[i]  );
        fprintf(ap,"\n%d %s %.8f", la->timePrint, sigla, la->probabilities[i] );
    }
    fclose(ap);
}

void LA_free( LearningAutomata **_la )
{

    LearningAutomata *la = *_la;

    LInt_free( la->lastNSelections );
    free( la->initialProbabilities);
    free( la->probabilities);

    free( la );

    *_la = NULL;

}

