#ifndef SOLUTION_H
#define SOLUTION_H

#include "instance.h"
#include "mode_set.h"
#include "rrusage.h"

typedef struct _Solution Solution;

/* returns an initialized object solution */
Solution *Sol_create( const Instance *inst );

/* Print the solution */
void Sol_print( const Solution *sol );

/* Write the solution to a file */
void Sol_write( const Solution *sol, char *file );

/* read a initial solution from a file */
//void Sol_read( Solution *sol, char *file );
int Sol_read( Solution *sol, char *file );

/* frees memory used by a solution */
void Sol_free( Solution **_solution );

/* returns the job, on sequence position (idx) */
int Sol_getSequence( const Solution *sol, int idx );

/* Returns the ModeSet of the solution */
ModeSet* Sol_getModeSet( const Solution *sol );

/* sets the allocation start time of a job */
void Sol_setStartJob(Solution *sol, int idxJob, int time);

/* returns the allocation start time of a job */
int Sol_getStartTime( const Solution *sol, int idxJob );

RRUsage* Sol_getRRUsage( Solution *sol);

unsigned int Sol_getSolutionHash( const Solution *sol);

const Instance *Sol_inst( const Solution *sol );

int *Sol_sequence( const Solution *sol );

int *Sol_startJobs( const Solution *sol );

void Sol_rebuild( Solution *solution);

void Sol_rebuild_opt( Solution *current, const Solution *solOLD);

void Sol_topSort( Solution *solution, int sequence[] );

int Sol_getPosJob(Solution *sol, int idxJob);

/* Clean the variables needed to rebuild a solution from the sequence */
void Sol_clearBuildSequence(Solution *sol);

/* Returns the cost of the solution */
Cost Sol_getCost( const Solution *sol );

/* Change the cost of the solution */
void Sol_setCost( Solution *sol, Cost cost );

Cost Sol_getTPD( const Solution *sol );
Cost Sol_getTMS( const Solution *sol );
/* Calculates the cost of the solution */
void Sol_calcCost( Solution *sol);

void Sol_fillSequence( Solution *sol );

int Sol_getMode( const Solution *sol, int job );

void Sol_cpy( Solution *target, const Solution *sol );

void Sol_setPosJob(Solution *sol, int idxJob, int value);

void Sol_setMinTimeSucc(Solution *sol, const Job* job, int time);

int Sol_getMinTime(Solution *sol, const Job* job, const Mode* mode, int minTime);

int Sol_getMinTime2(Solution *sol, int idxJob);


/*Evaluate a solution to GA*/
Cost Sol_evaluate( Solution *sol );

/* Function of mutation to GA*/
void Sol_firstMutation( Solution *sol );

/* Function of crossover to GA*/
void Sol_firstCrossover( Solution *sol );

#endif // SOLUTION_H





















/* change the entire set of modes solution */
//void Sol_setModesByProj( Solution *sol, const ModeSet *modeSet );

/* change the mode of one job in a solution */
//void Sol_setMode( Solution *sol, int job, int mode );

/* returns the mode of job in a solution */
//int Sol_getMode( const Solution *sol, int job );

/* sets the jobs on sequence*/
//void Sol_setSequence( Solution *sol, int idx, int idJob );
//void Sol_setPosJob(Solution *sol, int idxJob, int value);
//void Sol_cpySequence(const Solution *sol, int sequenceOLD[]);
//void Sol_cpyStartJobs(const Solution *sol, int startJobsOLD[]);

//void Sol_reconstruct( Solution *solution, int sequence[], int startJobs[], int modes[] );
/* Changes a job position in the sequence vector */
//void Sol_changeSequence( Solution *sol, int pos, int job, int sequence[] );

/* Creates a vector of random sequence */
//int *Sol_randomSequence( const Solution *sol);

/* Returns the job in a particular position of the sequence vector */
//int Sol_getJobOnSequence( const Solution *sol, int idx);

//RRUsage *Sol_rru( const Solution *sol );

/* Construct an initial solution based on the ordering of the heap by priority predecessors */
//void Sol_build( Solution *solution );

/* Construct an initial solution based on vector of sequence,
 * topological sort is applied to this vector if validTopSort
 * is false ( != 1)*/
//void Sol_buildBySequence( Solution *solution, int sequence[], ...);
//void Sol_clear( Solution *sol );
//int Sol_verifyPred( Solution *sol, int p1, int p2 ) ;
