/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef CONSTRUCTIVE_H
#define CONSTRUCTIVE_H

#include "instance.h"
#include "solution.h"

typedef struct _Constructive Constructive;

/*creates a solver to make the initial allocation of jobs*/
Constructive *Cons_create( const Instance *inst, Solution* sol, char** argv );

/*run the solver, creating the initial solution by random sequence whit the set of mode min*/
void Cons_run(Constructive *cons);
void Cons_runByProj(Constructive *cons);

void Cons_checkArgs(Constructive *cons, char **argv);

/*returns the LFA to lahc */
int Cons_getLfa(Constructive *cons);

/*returns the number of iterations without improvements to lahc */
int Cons_getIt(Constructive *cons);

/* frees memory used by Constructive */
void Cons_free( Constructive **_cons );

Cost Cons_getSolInitial(Constructive *cons);

#endif
