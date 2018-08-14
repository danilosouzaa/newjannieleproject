
/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <assert.h>
#include <string.h>
#include "mode_set.h"
#include "macros.h"

struct _ModeSet {
    /* Modes for each job on project*/
    int *modes;
    /* uses of each Non Renewable Resource on project*/
    int *usageNonRenewable;
    /* index of first job of project on instance*/
    int firstJob;
    int lastJob;
    int isEmpty;
    /* index of start non renewable resource of project on instance*/
    int startResN;
    /* number of non renewable resource of project*/
    int nResN;

    int nJobs;

    Cost cost;
    int infeasibility;

    const struct _Instance *inst;
};

ModeSet *Modes_create( const Instance* inst )
{
    assert( inst!=NULL );

    int nJobs = Inst_nJobs(inst);

    ModeSet *modeSet;
    ALLOCATE_INI(modeSet, ModeSet);
    ALLOCATE_VECTOR_INI( modeSet->modes, int, nJobs );
    ALLOCATE_VECTOR_INI( modeSet->usageNonRenewable, int, Inst_nResN(inst) );

    memset(modeSet->modes, -1, sizeof(int)*nJobs);

    modeSet->nResN = Inst_nResN(inst);
    modeSet->startResN = Inst_idxResNProj( inst, 0 );
    modeSet->inst = inst;
    modeSet->firstJob = 0;
    modeSet->lastJob = nJobs;
    modeSet->nJobs = nJobs;
    modeSet->cost = 0;
    modeSet->infeasibility = 0;
    modeSet->isEmpty = 1;


    /* for(int j = 0 ; j < nJobs ; j++)
         Modes_modify( modeSet, j, 0 );
     modeSet->isEmpty = 0;
    */
    return modeSet;
}

int Modes_job( const ModeSet *modeSet, int idxJob )
{
    assert( modeSet!=NULL );
    assert( idxJob >= 0 );
    assert( idxJob < Inst_nJobs( modeSet->inst ) );

    return modeSet->modes[idxJob-modeSet->firstJob];
}

int Modes_modifyAndVerify( ModeSet *modeSet, int idxJob, int idxNewMode )
{
    assert( modeSet!=NULL );
    assert( idxJob >= 0 );
    assert( idxJob < Inst_nJobs( modeSet->inst ) );

    const Job *job = Inst_job( modeSet->inst, idxJob );

    assert( idxNewMode >= 0 );
    assert( idxNewMode < Job_nModes(job) );

    int idxJobOnModeSet = idxJob - modeSet->firstJob;
    int mode = Modes_job(modeSet, idxJobOnModeSet);
    const Mode *modeOld = Job_mode( job, mode );
    const Mode *modeNew = Job_mode( job, idxNewMode );

    modeSet->cost -= Mode_duration( modeOld );
    modeSet->cost += Mode_duration( modeNew );

    int idxResource, idxResOnModeSet, capN=0;


    for( int rn=0; ( rn<Mode_nResN( modeOld ) ); ++rn ) {
        int overOld = 0, overNew = 0;

        idxResource = Mode_idxResN( modeOld, rn );
        idxResOnModeSet = idxResource-modeSet->startResN;
        capN = Inst_capResN( modeSet->inst, idxResource );

        overOld = ( ( modeSet->usageNonRenewable[idxResOnModeSet] - capN ) > 0 ) ?
                  ( modeSet->usageNonRenewable[idxResOnModeSet] - capN ) : 0;

        modeSet->usageNonRenewable[idxResOnModeSet] -= Mode_useResN( modeOld, rn );

        overNew = ( ( modeSet->usageNonRenewable[idxResOnModeSet] - capN ) > 0 ) ?
                  ( modeSet->usageNonRenewable[idxResOnModeSet] - capN ) : 0;

        modeSet->infeasibility -= ( overOld-overNew );
        modeSet->cost -= ( 1000* ( overOld-overNew ) );
    }

    modeSet->modes[idxJobOnModeSet] = Mode_index(modeNew);

    for( int rn=0; ( rn<Mode_nResN( modeNew ) ); ++rn ) {
        int overNew = 0, overOld =0;

        idxResource = Mode_idxResN( modeNew, rn );
        idxResOnModeSet = idxResource-modeSet->startResN;
        capN = Inst_capResN( modeSet->inst, idxResource );

        overOld = ( ( modeSet->usageNonRenewable[idxResOnModeSet] - capN ) > 0 ) ?
                  ( modeSet->usageNonRenewable[idxResOnModeSet] - capN ) : 0;

        modeSet->usageNonRenewable[idxResOnModeSet] += Mode_useResN( modeNew, rn );

        overNew = ( ( modeSet->usageNonRenewable[idxResOnModeSet] - capN ) > 0 ) ?
                  ( modeSet->usageNonRenewable[idxResOnModeSet] - capN ) : 0;

        modeSet->infeasibility += ( overNew-overOld );
        modeSet->cost += ( 1000* ( overNew-overOld ) );
    }

    if(modeSet->infeasibility) return 0;

    return 1;

}


int Modes_verify( ModeSet *modeSet, int idxJob, int idxNewMode )
{
    assert( modeSet!=NULL );
    assert( idxJob >= 0 );
    assert( idxJob < Inst_nJobs( modeSet->inst ) );

    const Job *job = Inst_job( modeSet->inst, idxJob );

    assert( idxNewMode >= 0 );
    assert( idxNewMode < Job_nModes(job) );
    int* resources;

    ALLOCATE_VECTOR_INI( resources, int, Inst_nResN(modeSet->inst) );
    memcpy(resources,modeSet->usageNonRenewable, sizeof(int)* Inst_nResN(modeSet->inst));

    int infeasibility = 0;

    int idxJobOnModeSet = idxJob - modeSet->firstJob;
    int mode = Modes_job(modeSet, idxJobOnModeSet);

    if(mode!=-1) {
        const Mode *modeOld = Job_mode( job, mode );
        const Mode *modeNew = Job_mode( job, idxNewMode );
        int idxResource, idxResOnModeSet, capN=0;

        for( int rn=0; ( rn<Mode_nResN( modeOld ) ); ++rn ) {
            int overOld = 0, overNew = 0;

            idxResource = Mode_idxResN( modeOld, rn );
            idxResOnModeSet = idxResource-modeSet->startResN;
            capN = Inst_capResN( modeSet->inst, idxResource );

            overOld = ( ( resources[idxResOnModeSet] - capN ) > 0 ) ?
                      ( resources[idxResOnModeSet] - capN ) : 0;

            resources[idxResOnModeSet] -= Mode_useResN( modeOld, rn );

            overNew = ( ( resources[idxResOnModeSet] - capN ) > 0 ) ?
                      ( resources[idxResOnModeSet] - capN ) : 0;

            // printf("OverOld %d, OverNew %d", overOld, overNew );
            infeasibility -= ( overOld-overNew );
        }

        for( int rn=0; ( rn<Mode_nResN( modeNew ) ); ++rn ) {
            int overNew = 0, overOld =0;

            idxResource = Mode_idxResN( modeNew, rn );
            idxResOnModeSet = idxResource-modeSet->startResN;
            capN = Inst_capResN( modeSet->inst, idxResource );

            overOld = ( ( resources[idxResOnModeSet] - capN ) > 0 ) ?
                      ( resources[idxResOnModeSet] - capN ) : 0;

            resources[idxResOnModeSet] += Mode_useResN( modeNew, rn );

            overNew = ( ( resources[idxResOnModeSet] - capN ) > 0 ) ?
                      ( resources[idxResOnModeSet] - capN ) : 0;

            // printf("OverOld %d, OverNew %d", overOld, overNew );
            infeasibility += ( overNew-overOld );
        }

    } else {

        const Mode *modeNew = Job_mode( job, idxNewMode );


        for( int rn=0; ( rn<Mode_nResN( modeNew ) ); ++rn ) {
            int overNew = 0, overOld =0;
            int idxResource, idxResOnModeSet, capN=0;

            idxResource = Mode_idxResN( modeNew, rn );
            idxResOnModeSet = idxResource-modeSet->startResN;
            capN = Inst_capResN( modeSet->inst, idxResource );

            overOld = ( ( resources[idxResOnModeSet] - capN ) > 0 ) ?
                      ( resources[idxResOnModeSet] - capN ) : 0;

            //printf("\nOld resources[idxResOnModeSet] %d - capN%d : %d\n", resources[idxResOnModeSet],capN, resources[idxResOnModeSet] - capN );
            resources[idxResOnModeSet] += Mode_useResN( modeNew, rn );

            overNew = ( ( resources[idxResOnModeSet] - capN ) > 0 ) ?
                      ( resources[idxResOnModeSet] - capN ) : 0;

            //  printf("\nNew resources[idxResOnModeSet] %d - capN%d : %d\n", resources[idxResOnModeSet],capN, resources[idxResOnModeSet] - capN );

            //  printf("OverOld %d, OverNew %d", overOld, overNew );
            infeasibility += ( overNew-overOld );
        }

    }

    free(resources);

    //  printf("\ninfeasibility of Non-Renewable Resources %d \n", infeasibility);

    if(infeasibility) return 0;

    return 1;

}

void Modes_modify( ModeSet *modeSet, int idxJob, int idxNewMode)
{
    assert( modeSet!=NULL );
    assert( idxJob >= 0 );
    assert( idxJob < Inst_nJobs( modeSet->inst ) );

    const Job *job = Inst_job( modeSet->inst, idxJob );

    assert( idxNewMode >= 0 );
    assert( idxNewMode < Job_nModes(job) );

    int idxJobOnModeSet = idxJob - modeSet->firstJob;

    int mode = Modes_job(modeSet, idxJobOnModeSet);

    if(mode!=-1) {

        const Mode *modeOld = Job_mode( job, mode );
        const Mode *modeNew = Job_mode( job, idxNewMode );

        modeSet->cost -= Mode_duration( modeOld );
        modeSet->cost += Mode_duration( modeNew );

        int idxResource, idxResOnModeSet;


        for( int rn=0; ( rn<Mode_nResN( modeOld ) ); ++rn ) {

            idxResource = Mode_idxResN( modeOld, rn );
            idxResOnModeSet = idxResource-modeSet->startResN;
            modeSet->usageNonRenewable[idxResOnModeSet] -= Mode_useResN( modeOld, rn );

        }

        modeSet->modes[idxJobOnModeSet] = Mode_index(modeNew);

        for( int rn=0; ( rn<Mode_nResN( modeNew ) ); ++rn ) {

            idxResource = Mode_idxResN( modeNew, rn );
            idxResOnModeSet = idxResource-modeSet->startResN;
            modeSet->usageNonRenewable[idxResOnModeSet] += Mode_useResN( modeNew, rn );

            int capN = Inst_capResN( modeSet->inst, idxResource );
            if(modeSet->usageNonRenewable[idxResOnModeSet] > capN)
                printf("Infeasible to non renewable resource value %d capacity %d", modeSet->usageNonRenewable[idxResOnModeSet], capN );
        }


    } else {

        const Mode *modeNew = Job_mode( job, idxNewMode );
        modeSet->cost += Mode_duration( modeNew );

        modeSet->modes[idxJobOnModeSet] = Mode_index(modeNew);

        for( int rn=0; ( rn<Mode_nResN( modeNew ) ); ++rn ) {

            int idxResource, idxResOnModeSet;

            idxResource = Mode_idxResN( modeNew, rn );
            idxResOnModeSet = idxResource-modeSet->startResN;
            modeSet->usageNonRenewable[idxResOnModeSet] += Mode_useResN( modeNew, rn );

            int capN = Inst_capResN( modeSet->inst, idxResource );
            if(modeSet->usageNonRenewable[idxResOnModeSet] > capN)
                printf("Infeasible to non renewable resource value %d capacity %d", modeSet->usageNonRenewable[idxResOnModeSet], capN );
        }
    }
}

int *Modes_getModes(const ModeSet *modeSet)
{
    assert( modeSet != NULL );

    return modeSet->modes;

}

void Modes_cpy( ModeSet *target, const ModeSet *modeSet )
{
    assert( target!=NULL );
    assert( modeSet!=NULL );

    target->nResN = modeSet->nResN;
    target->nJobs = modeSet->nJobs;
    target->startResN = modeSet->startResN;
    target->cost = modeSet->cost;
    target->infeasibility = modeSet->infeasibility;
    //target->inst = modeSet->inst;
    target->firstJob = modeSet->firstJob;
    target->isEmpty = modeSet->isEmpty;

    COPY_VECTOR( target->modes, modeSet->modes, int, modeSet->nJobs );
    COPY_VECTOR( target->usageNonRenewable, modeSet->usageNonRenewable, int, modeSet->nResN);

}

int Modes_inf(const ModeSet *modeSet)
{
    assert( modeSet!=NULL );

    return modeSet->infeasibility;
}

Cost Modes_cost(const ModeSet *modeSet)
{
    assert( modeSet!=NULL );

    return modeSet->cost;
}

void Modes_free( ModeSet **_modeSet )
{
    ModeSet *modeSet = *_modeSet;

    free( modeSet->modes );
    free( modeSet->usageNonRenewable );

    free( modeSet );

    *_modeSet = NULL;
}













/*
const Instance *Modes_inst( const ModeSet *modeSet )
{
    assert( modeSet!=NULL );

    return modeSet->inst;
}

void Modes_reconstructsModes( ModeSet *modeSet, const int modes[] )
{
    assert( modeSet!=NULL );
    assert( sizeof( modeSet->modes ) == sizeof( modes ) );

    int idxJob;
    for(int j = 0; ( j<Modes_nJobs( modeSet ) ); ++j ) {
        idxJob = j+Modes_firstJob( modeSet );
        if( Modes_job( modeSet, idxJob ) != modes[j] )
            Modes_modifyCount( modeSet, idxJob, modes[j]);
    }
}
*/





/*
ModeSet *Modes_createForMSP( const Instance* inst, int firstJob, int lastJob, int nResN )
{
    assert( inst != NULL );
    assert( firstJob >= 0 );
    assert( lastJob > firstJob );
    assert( lastJob < Inst_nJobs( inst ) );
    assert( firstJob < lastJob );

    int nJobs = lastJob-firstJob+1;

    ModeSet *modeSet;
    ALLOCATE_INI( modeSet, ModeSet );
    ALLOCATE_VECTOR( modeSet->modes, int, nJobs );
    ALLOCATE_VECTOR_INI( modeSet->usageNonRenewable, int, nResN );

    modeSet->nResN = nResN;

    const Job *job = Inst_job( inst, firstJob );
    int p = Job_project( job );

    modeSet->startResN = Inst_idxResNProj( inst, p );
    modeSet->inst = inst;
    modeSet->firstJob = firstJob;
    modeSet->lastJob = lastJob;
    modeSet->nJobs = nJobs;

    for ( int j=firstJob; ( j<=lastJob ); ++j ) {
        const Job *job = Inst_job( inst, j );

        for(int m = 0; m < Job_nModes(job); m++) {
            const Mode *mode = Job_mode(job,m);
            if(Mode_isFeasible(inst, mode)) {

                const Mode *mode = Job_mode( job, m );
                modeSet->modes[j-firstJob] = Mode_index(mode);
                modeSet->cost += Mode_duration( mode );
                for( int rn=0; ( rn<Mode_nResN( mode ) ); ++rn) {
                    int idxResource = Mode_idxResN( mode, rn );
                    modeSet->usageNonRenewable[idxResource-modeSet->startResN] += Mode_useResN( mode, rn );
                }
                break;
            }

        }

    }

    int capN = 0;
    modeSet->infeasibility = 0;
    int idxResource = modeSet->startResN;
    for( int rn = 0; ( rn < nResN ); ++rn, ++idxResource ) {
        capN = Inst_capResN( inst, idxResource );
        modeSet->infeasibility += ( ( modeSet->usageNonRenewable[rn] ) - capN > 0 ) ?
                                  ( modeSet->usageNonRenewable[rn] - capN ) : 0;
    }

    modeSet->cost += ( 1000*modeSet->infeasibility );
    modeSet->isEmpty = 0;

    return modeSet;
}

ModeSet *Modes_createForMS( const Instance* inst)
{
    assert( inst != NULL );

    int nJobs = Inst_nJobs(inst);

    ModeSet *modeSet;
    ALLOCATE_INI( modeSet, ModeSet );
    ALLOCATE_VECTOR_INI( modeSet->modes, int, nJobs );
    ALLOCATE_VECTOR_INI( modeSet->usageNonRenewable, int, Inst_nResN(inst) );

    modeSet->nResN = Inst_nResN(inst);

    modeSet->startResN = Inst_idxResNProj( inst, 0 );
    modeSet->inst = inst;
    modeSet->firstJob = 0;
    modeSet->lastJob = nJobs-1;
    modeSet->nJobs = nJobs;

    for ( int j=0; ( j<= modeSet->lastJob ); ++j ) {
        const Job *job = Inst_job( inst, j );

        for(int m = 0; m < Job_nModes(job); m++) {
            const Mode *mode = Job_mode(job,m);
            if(Mode_isFeasible(inst, mode)) {

                const Mode *mode = Job_mode( job, m );
                modeSet->modes[j] = Mode_index(mode);
                modeSet->cost += Mode_duration( mode );
                for( int rn=0; ( rn<Mode_nResN( mode ) ); ++rn) {
                    int idxResource = Mode_idxResN( mode, rn );
                    modeSet->usageNonRenewable[idxResource] += Mode_useResN( mode, rn );
                }
                break;
            }

        }

    }

    int capN = 0;
    modeSet->infeasibility = 0;
    for( int rn = 0; ( rn < Inst_nResN(inst) ); ++rn ) {
        capN = Inst_capResN( inst, rn );
        modeSet->infeasibility += ( ( modeSet->usageNonRenewable[rn] ) - capN > 0 ) ?
                                  ( modeSet->usageNonRenewable[rn] - capN ) : 0;
    }

    modeSet->cost += ( 1000*modeSet->infeasibility );

    modeSet->isEmpty = 0;

    return modeSet;
}
*/

/*
int Modes_firstJob( const ModeSet *modeSet )
{
    assert( modeSet!=NULL );

    return modeSet->firstJob;
}

int Modes_lastJob( const ModeSet *modeSet )
{
    assert( modeSet!=NULL );

    return modeSet->lastJob;
}

int Modes_nJobs( const ModeSet *modeSet )
{
    assert( modeSet!=NULL );

    return modeSet->nJobs;
}


int Modes_isEmpty( const ModeSet *modeSet )
{
    assert( modeSet!=NULL );

    return modeSet->isEmpty;
}


*/
void Modes_isntEmpty( ModeSet *modeSet )
{
    assert( modeSet!=NULL );

    modeSet->isEmpty = 0;
}

/*
void Modes_modifyCount( ModeSet *modeSet, int idxJob, int idxNewMode )
{
    assert( modeSet!=NULL );
    assert( idxJob >= 0 );

    const Job *job = Inst_job( modeSet->inst, idxJob );

    assert( idxNewMode >= 0 );
    assert( idxNewMode < Job_nModes(job) );

    int idxJobOnModeSet = idxJob - modeSet->firstJob;
    int mode = Modes_job(modeSet, idxJob);
    const Mode *modeOld = Job_mode( job, mode );
    const Mode *modeNew = Job_mode( job, idxNewMode );

    modeSet->cost -= Mode_duration( modeOld );
    modeSet->cost += Mode_duration( modeNew );

    int idxResource, idxResOnModeSet, capN=0;


    for( int rn=0; ( rn<Mode_nResN( modeOld ) ); ++rn ) {
        int overOld = 0, overNew = 0;

        idxResource = Mode_idxResN( modeOld, rn );
        idxResOnModeSet = idxResource-modeSet->startResN;
        capN = Inst_capResN( modeSet->inst, idxResource );

        overOld = ( ( modeSet->usageNonRenewable[idxResOnModeSet] - capN ) > 0 ) ?
                  ( modeSet->usageNonRenewable[idxResOnModeSet] - capN ) : 0;

        modeSet->usageNonRenewable[idxResOnModeSet] -= Mode_useResN( modeOld, rn );

        overNew = ( ( modeSet->usageNonRenewable[idxResOnModeSet] - capN ) > 0 ) ?
                  ( modeSet->usageNonRenewable[idxResOnModeSet] - capN ) : 0;

        modeSet->infeasibility -= ( overOld-overNew );
        modeSet->cost -= ( 1000* ( overOld-overNew ) );
    }

    modeSet->modes[idxJobOnModeSet] = Mode_index(modeNew);

    for( int rn=0; ( rn<Mode_nResN( modeNew ) ); ++rn ) {
        int overNew = 0, overOld =0;

        idxResource = Mode_idxResN( modeNew, rn );
        idxResOnModeSet = idxResource-modeSet->startResN;
        capN = Inst_capResN( modeSet->inst, idxResource );

        overOld = ( ( modeSet->usageNonRenewable[idxResOnModeSet] - capN ) > 0 ) ?
                  ( modeSet->usageNonRenewable[idxResOnModeSet] - capN ) : 0;

        modeSet->usageNonRenewable[idxResOnModeSet] += Mode_useResN( modeNew, rn );

        overNew = ( ( modeSet->usageNonRenewable[idxResOnModeSet] - capN ) > 0 ) ?
                  ( modeSet->usageNonRenewable[idxResOnModeSet] - capN ) : 0;

        modeSet->infeasibility += ( overNew-overOld );
        modeSet->cost += ( 1000* ( overNew-overOld ) );
    }


}
*/


/*
void Modes_fillModes(const ModeSet *modeSet, int modes[])
{

    for( int j = 0; j<Modes_nJobs( modeSet ); j++ ) {
        int idxJob = Modes_firstJob( modeSet )+j;
        modes[j] = Modes_job( modeSet, idxJob );
    }
}
*/
/*
void Modes_cpyByProj( ModeSet *target, const ModeSet *modeSet )
{
    assert( modeSet!=NULL );

    const Job *job = Inst_job( modeSet->inst, modeSet->firstJob );
    const Project *p = Inst_project( modeSet->inst, Job_project( job ) );

    target->nResN = modeSet->nResN;
    target->nJobs = modeSet->nJobs;
    target->startResN = modeSet->startResN;
    target->cost = modeSet->cost;
    target->infeasibility = modeSet->infeasibility;
    target->inst = modeSet->inst;
    target->firstJob = modeSet->firstJob;
    target->isEmpty = modeSet->isEmpty;

    for( int j=0; j<Project_nJobs( p ); ++j )
        target->modes[j] = modeSet->modes[j];

    for( int rn=0; rn<modeSet->nResN; ++rn )
        target->usageNonRenewable[rn] = modeSet->usageNonRenewable[rn];
}


void Modes_cpyToSolByProj( ModeSet *modeSet, const ModeSet *bestModes )
{
    assert( modeSet!=NULL );

    ModeSet *target = modeSet;

    int idxOnInstance = Modes_firstJob( bestModes );
    int idxLastOnInstance = Modes_lastJob( bestModes );
    int idxJobOnMode = 0;

    target->cost += bestModes->cost;
    target->infeasibility += bestModes->infeasibility;
    target->isEmpty = bestModes->isEmpty;

    for( int j=idxOnInstance; j<idxLastOnInstance; ++j, ++idxJobOnMode )
        target->modes[j] = bestModes->modes[idxJobOnMode];
    for( int rn=0; rn<bestModes->nResN; ++rn )
        target->usageNonRenewable[bestModes->startResN+rn] = bestModes->usageNonRenewable[rn];
}

void Modes_clearByProjec( ModeSet *modeSet )
{
    assert( modeSet!=NULL );

    memset( modeSet->modes, 0, (int) sizeof(int) * (int) Modes_nJobs( modeSet ) );
    memset( modeSet->usageNonRenewable, 0, (int) sizeof(int) * ( Inst_nResN( modeSet->inst )
            / Inst_nProjects( modeSet->inst ) ) );

    modeSet->cost = 0;
    modeSet->infeasibility = 0;
}

*/

/*
void Modes_clear( ModeSet *modeSet )
{
    assert( modeSet!=NULL );

    memset( modeSet->modes, 0, (int) sizeof(int) * (int) Modes_nJobs( modeSet ) );
    memset( modeSet->usageNonRenewable, 0, (int) sizeof(int) *Inst_nResN( modeSet->inst ) );

    modeSet->cost = 0;
    modeSet->isEmpty = 1;
    modeSet->infeasibility = 0;
}
*/

