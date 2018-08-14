/*
 * proj_decomp.c
 * Copyright (C) 2015 H.G. Santos <haroldo.santos@gmail.com>
 */

#include <stdio.h>
#include <stdlib.h>
#include "proj_decomp.h"
#include "macros.h"
#include "lp.h"
#include "vec_str.h"
#include "mip_project.h"

struct _ProjDecomp {
    const Instance *inst;
    LinearProgram *master;
    int **idxRowCap;
    int tpdSum;

    /* pricing subproblems */
    MIPProject **mipP;

    int maxTime;
};

static void prjdc_create_dummy_variables( ProjDecomp *prjdc );
static void prjdc_create_constraint_convexification( ProjDecomp *prjdc );
static void prjdc_create_constraint_resources( ProjDecomp *prjdc );
static void prjdc_create_master( ProjDecomp *prjdc );
static double prjc_dummy_cost( const ProjDecomp *prjdc );
static void create_pricing_problems( ProjDecomp *prjdc );

ProjDecomp *prjdc_create( const Instance *instance, int tpdSum )
{
    ProjDecomp *prjdc;
    ALLOCATE_INI( prjdc, ProjDecomp );
    prjdc->inst = instance;
    prjdc->tpdSum = tpdSum;

    prjdc->maxTime = 0;
    for ( int i=0 ; (i<Inst_nProjects(instance)) ; ++i ) {
        const Project *proj = Inst_project( instance, i );
        prjdc->maxTime = MAX( prjdc->maxTime, Project_releaseDate(proj)+Project_criticalPath(proj)+tpdSum );
    }

    prjdc_create_master( prjdc );
    create_pricing_problems( prjdc );


    lp_write_lp( prjdc->master, "master.lp" );

    return prjdc;
}

void prjdc_create_master( ProjDecomp *prjdc )
{
    prjdc->master = lp_create();
    prjdc_create_dummy_variables( prjdc );
    prjdc_create_constraint_convexification( prjdc );
    prjdc_create_constraint_resources( prjdc );
}


void prjdc_create_dummy_variables( ProjDecomp *prjdc )
{
    VecStr *names = VStr_create( STR_SIZE );

    for ( int i=0 ; i<Inst_nProjects(prjdc->inst) ; ++i ) {
        char vname[STR_SIZE];
        sprintf( vname, "lbda(%d)", i );
        VStr_pushBack( names, vname );
    }

    double obj[ Inst_nProjects(prjdc->inst) ];
    FILL( obj, 0, Inst_nProjects(prjdc->inst), prjc_dummy_cost(prjdc) );

    lp_add_bin_cols( prjdc->master, Inst_nProjects(prjdc->inst), obj, VStr_ptr(names) );

    VStr_free( &names );
}

void prjdc_create_constraint_convexification( ProjDecomp *prjdc )
{
    for ( int j=0 ; (j<Inst_nProjects(prjdc->inst)) ; ++j ) {
        char rname[STR_SIZE];
        sprintf( rname, "project(%d)", j );
        int idx = j;
        double coef = 1.0;
        lp_add_row( prjdc->master, 1, &idx, &coef, rname, 'E', 1.0 );
    }
}

void prjdc_free( ProjDecomp **_prjdc )
{
    ProjDecomp *prjdc = *_prjdc;
    lp_free( &prjdc->master );

    const int nProjects = Inst_nProjects( prjdc->inst );
    for ( int i=0 ; (i<nProjects) ; ++i )
        MipP_free( &(prjdc->mipP[i]) );

    free( prjdc->idxRowCap[0] );
    free( prjdc->idxRowCap );

    free( prjdc->mipP );

    free( *_prjdc );
    *_prjdc = NULL;
}

double prjc_dummy_cost( const ProjDecomp *prjdc )
{
    return 10000;
}

void create_pricing_problems( ProjDecomp *prjdc )
{
    const Instance *inst = prjdc->inst;
    const int nProjects = Inst_nProjects(inst);
    ALLOCATE_VECTOR( prjdc->mipP, MIPProject*, nProjects );

    for ( int i=0 ; (i<nProjects) ; ++i )
        prjdc->mipP[i] = MipP_create( inst, i, prjdc->tpdSum );
}

void prjdc_create_constraint_resources( ProjDecomp *prjdc )
{
    const Instance *inst = prjdc->inst;
    LinearProgram *master = prjdc->master;

    ALLOCATE_VECTOR( prjdc->idxRowCap, int*, Inst_nResRGlobal(inst) );
    ALLOCATE_VECTOR_INI( prjdc->idxRowCap[0], int, Inst_nResRGlobal(inst)*prjdc->maxTime );
    for ( int i=1 ; (i<Inst_nResRGlobal(inst)) ; ++i )
        prjdc->idxRowCap[i] = prjdc->idxRowCap[i-1] + prjdc->maxTime;

    FILL( prjdc->idxRowCap[0], 0, Inst_nResRGlobal(inst)*prjdc->maxTime, -1 );

    for ( int i=0 ; (i<Inst_nResRGlobal(inst)) ; ++i ) {
        for ( int j=0 ; (j<prjdc->maxTime) ; ++j ) {
            char rname[STR_SIZE];
            sprintf( rname, "capacity(%d,%d)", i,j );
            int idx;
            double coef;
            prjdc->idxRowCap[i][j] = lp_rows( master );
            lp_add_row( master, 0, &idx, &coef, rname, 'L', Inst_capResR(inst,i) );
        }
    }
}
