/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef MODE_SET_H
#define MODE_SET_H

#include "instance.h"
#include "macros.h"

typedef struct _ModeSet ModeSet;

ModeSet *Modes_create( const Instance* inst );

/* returns the mode of a job on modeSet, the idxJob corresponds to the instance*/
int Modes_job( const ModeSet *modeSet, int idxJob );

int *Modes_getModes(const ModeSet *modeSet);

/* changes the mode of a job on modeSet, the idxJob corresponds to the instance,
   compute de cost (min total duration), it isn't the cost of the solution.
   don't count the infeasibility */
void Modes_modify( ModeSet *modeSet, int idxJob, int idxNewMode);
int Modes_verify( ModeSet *modeSet, int idxJob, int idxNewMode );
void Modes_isntEmpty( ModeSet *modeSet );
/* changes the mode of a job on modeSet, the idxJob corresponds to the instance,
   count the infeasibility,
   compute the cost (min total duration), it isn't the cost of the solution.
   penalizes (1000) the cost by infeasibility,
   return if is feasible or infeasible */
int Modes_modifyAndVerify( ModeSet *modeSet, int idxJob, int idxNewMode );


/* Makes a copy of modeSet  to target */
void Modes_cpy( ModeSet *target, const ModeSet *modeSet );


/*returns the cost allocations of modes of jobs on project,
considering duration and infeasibility of non-renewable resources*/
Cost Modes_cost( const ModeSet *modeSet );

/*returns the infeasibility of non-renewable resources*/
int Modes_inf( const ModeSet *modeSet );

/* frees memory used by a ModeSet */
void Modes_free( ModeSet **_modeSet );

#endif


















/* creates and initializes a set of index modes for jobs.
This structure allows working with intervals corresponding to the jobs on project,
but can be used in solution with intervals of total jobs on instance,
already calculates the infeasibility.*/
//ModeSet *Modes_createForMS( const Instance* inst);
//ModeSet *Modes_createForMSP( const Instance* inst, int firstJob, int lastJob, int nResN );

/* just creates a set of index modes for jobs, but doesn't initializes the set.
This structure allows working with intervals corresponding to the jobs on project,
but can be used in solution with intervals of total jobs on instance
*/
//ModeSet *Modes_create( const Instance* inst, int firstJob, int lastJob, int nResN );

//void Modes_isntEmpty( ModeSet *modeSet );

/* For MODE SET
   changes the mode of a job on modeSet, the idxJob corresponds to the instance,
   count the infeasibility,
   compute the cost (min total duration), it isn't the cost of the solution.
   penalizes (1000) the cost by infeasibility*/
//void Modes_modifyCount( ModeSet *modeSet, int idxJob, int idxNewMode );
//void Modes_cpyByProj( ModeSet *target, const ModeSet *modeSet );
/* Makes a copy of modeSet  to modeSet on Solution in the correct position. */
//void Modes_cpyToSolByProj( ModeSet *modeSet, const ModeSet *bestModes );

/*clean the element*/
//void Modes_clearByProjec( ModeSet *modeSet );
//void Modes_clear( ModeSet *modeSet );

/*returns the number of jobs on modeSet */
//int Modes_nJobs( const ModeSet *modeSet );

/*returns the first jobs on modeSet */
//int Modes_firstJob( const ModeSet *modeSet );

/*returns the last jobs on modeSet */
//int Modes_lastJob( const ModeSet *modeSet );

/*returns the instance of modeSet */
//const Instance *Modes_inst( const ModeSet *modeSet );

/* extracts the order of modes for an extra vector. */
//void Modes_fillModes( const ModeSet *modeSet, int modes[] );

/* reconstructs a ModeSet from a set of modes, compared to current ModeSet. */
//void Modes_reconstructsModes( ModeSet *modeSet, const int modes[] );

