/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef SA_H
#define SA_H

#include "instance.h"
#include "neighborhood.h"
#include "solution.h"
#include "sa.h"


typedef struct _SA SA;

/*creates a solver to make the allocation of jobs*/
SA *SA_create( const Instance *inst, Neighborhood* neighborhood,  Solution* sol, char **argv, int argc );

/* -iterT -T -SAmax, -alpha*/
void SA_checkArgs(SA *sa, char **argv, int argc);

double SA_getAlpha(SA *_sa);

int SA_getT(SA *_sa);

int SA_getSAmax(SA *_sa);

void SA_increasingTransitivity(SA* sa, Neighborhood* neighborhood);
void SA_increasingResidence(SA* sa, Solution* current);

/*run the solver*/
void SA_run(SA *sa, Neighborhood* neighborhood, double timeRem, Test *test);

/* frees memory used by SA */
void SA_free( SA **_sa );

#endif

