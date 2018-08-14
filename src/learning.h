/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                              of Santos, H.G., Toffolo, T.A.M. and Baltar, D.
 */

#ifndef LEARNING_HEADER
#define LEARNING_HEADER

/* opaque data structures */
typedef struct _LearningAutomata LearningAutomata;

/* creates and returns an empty LearningAutomata */
LearningAutomata *LA_create( const int nIntensities, const long double intensities[], char **argv, int argc );

void LA_checkArgs(LearningAutomata *la, char **argv, int argc);

/* returns the next neighborhood index based on the what has been learned */
int LA_next( LearningAutomata *la );

/* resets the learning automata to its original state */
void LA_reset( LearningAutomata *la );

/* updates (gives feedback to) the learning automata */
/* the obtainedReinforcement x should be a value such that 0.0 <= x <= 1.0 */
void LA_update( LearningAutomata *la, const double obtainedReinforcement );
double LA_getLearningRate( const LearningAutomata *la );

/* returns the number of iterations of the learning automata */
int LA_getIters( const LearningAutomata *la );

/* returns the last selection (choice) of the learning automata */
int LA_getLastSelection( const LearningAutomata *la );

/* returns the current learning rate */
double LA_getLearningRate( const LearningAutomata *la );

int LA_getResetInterval( const LearningAutomata *la );
/* sets the probabilities of the learning using the intensities value */
void LA_setIntensities( LearningAutomata *la, double intensities[] );

/* sets the learning rate */
void LA_setLearningRate( LearningAutomata *la, double learningRate );

/* sets the probabilities of the learning automata */
void LA_setProbabilities( LearningAutomata *la, double probabilities[] );

/* runs a simple test to check whether the learning automata is working well */
void LA_runSimpleTest();

void LA_printProbabilities( LearningAutomata *la, float pesw);
void LA_setTime( LearningAutomata *la);
void LA_free( LearningAutomata **_la );

#endif /* LEARNING_HEADER */
