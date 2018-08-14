/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "instance.h"
#include "ms_solver_lahc.h"
#include "mode_set.h"
#include "constructive.h"
#include "solution.h"

struct _Constructive {
    Solution *solInitial;

    int lfa; //will be removed after find the best parameters.
    int it;

    const struct _Instance *inst;
};

void Cons_checkArgs(Constructive *cons, char **argv)
{

    int n=0;
    while(argv[++n] != NULL) {
        if (strcmp(argv[n],"-lfa") == 0) {
            n++;
            cons->lfa = atoi(argv[n]);
            continue;
        }
        if (strcmp(argv[n],"-it") == 0) {
            n++;
            cons->it = atoi(argv[n]);
            continue;
        }
    }
}



Constructive *Cons_create( const Instance *inst, Solution* sol, char** argv )
{
    Constructive* cons;

    ALLOCATE_INI( cons, Constructive );

    cons->inst = inst;

    cons->solInitial = sol;

    cons->lfa = 10;

    cons->it = 1000;

    Cons_checkArgs(cons, argv);

    return cons;
}

void Cons_free( Constructive **_cons )
{

    Constructive *cons = *_cons;

    free( cons );

    *_cons = NULL;
}

int Cons_getLfa(Constructive *cons)
{

    assert(cons != NULL);

    return cons->lfa;
}

int Cons_getIt(Constructive *cons)
{

    assert(cons != NULL);

    return cons->it;
}

const Instance *Cons_getInst(Constructive *cons)
{

    assert(cons != NULL);

    return cons->inst;
}

void Cons_runByProj(Constructive *cons)
{

    /*creating the set of mode min*/
    MS_runByProj(cons->inst, cons->solInitial, Cons_getLfa(cons), Cons_getIt(cons));

    /* creating the initial solution by random sequence */
    int *sequence = Sol_randomSequence(cons->solInitial);
    Sol_buildBySequence(cons->solInitial, sequence);

    free(sequence);
}

void Cons_run(Constructive *cons)
{

    /*creating the set of mode min*/
    MSSolverLAHC *msLahc = MS_create(Cons_getInst(cons), Cons_getLfa(cons));
    MS_run(msLahc, Cons_getIt(cons));

    Modes_cpy(Sol_getModeSet(cons->solInitial), MS_bestModes(msLahc));

    /* creating the initial solution by random sequence */
    int *sequence = Sol_randomSequence(cons->solInitial);
    Sol_buildBySequence(cons->solInitial, sequence);

    MS_free( &msLahc );
    free(sequence);
}

Cost Cons_getSolInitial(Constructive *cons)
{
    return Sol_getCost(cons->solInitial);
}
