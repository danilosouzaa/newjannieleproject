/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 *
 * rrusage: controls the usage of renewable resources in a solution
 *
 */

#ifndef RRUSAGE_H
#define RRUSAGE_H

#include "instance.h"

typedef struct _RRUsage RRUsage;


//int RRU_getTsUsed(RRUsage *rru);
//void RRU_setTsCap(RRUsage *rru, int value);

/* creates a RRUsage object to control the usage of renewable resources */
RRUsage *RRU_create( const Instance *inst );

/* finds the first time instant where a job with the following
 * mode can be allocated from startingTime on */
int RRU_find( const RRUsage *rru, const Mode *mode, int startingTime );
//int RRU_getTsCap(RRUsage *rru);

/* allocate resources for mode at time start */
void RRU_allocate( RRUsage *rru, const Mode *mode, int start );
void  RRU_deallocate( RRUsage *rru, const Mode *mode, int start);
void RRU_clear_opt( RRUsage *rru, int j, const int seq[], const int starts[], const int modes[], int nJobs);

/* returns usage of resource idxRR at time */
int RRU_usage( const RRUsage *rru, int time, int idxRR);

void RRU_cpy( RRUsage *target, const RRUsage *rru);
void RRU_copy_part( RRUsage *rru, int j, const int seq[], const int starts[], const int modes[]);
void RRU_clear( RRUsage *rru);

/* frees data structure */
void RRU_free( RRUsage **_rru );

#endif /* !RRUSAGE_H */

