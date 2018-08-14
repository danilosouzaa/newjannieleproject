/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <time.h>
#include <assert.h>
#include "ils.h"

struct _ILS {

    Solution *solIni;
    Neighborhood *neighborhood;

    const struct _Instance *inst;

};


ILS *ILS_create( const Instance *inst, Solution* sol )
{
    ILS* ils;

    ALLOCATE_INI( ils, ILS );

    ils->solIni = sol;

    ils->neighborhood = Neighbor_createMS(inst);

    return ils;
}

void ILS_run(ILS *ils)
{

}

void ILS_free( ILS **_ils )
{

    ILS *ils = *_ils;

    Sol_free( &ils->solIni );

    Neighbor_free( &ils->neighborhood );

    *_ils = NULL;

}
