/*
 * mip.c
 * Copyright (C) 2016 H.G. Santos <haroldo.santos@gmail.com>
 */

#include <stdio.h>
#include <stdlib.h>
#include "instance.h"
#include "solution.h"
#include "mip_compact.h"

int main( int argc, const char **argv )
{
    if (argc<6) {
        fprintf( stderr, "usage: mip instDir inst solDir sol tpdSum\n");
        exit( EXIT_FAILURE );
    }

    Instance *inst = Inst_create( argv[1], argv[2] );
    int tpdSum = atoi(argv[5]);

    Solution *sol = Sol_create( inst );
    {
        char solFName[FILE_NAME_SIZE];
        sprintf( solFName, "%s/%s", argv[3], argv[4] );
        Sol_read( sol, solFName );
        printf("loaded solution with tpd %ld\n", Sol_getTPD(sol) );
    }

    MIPCompact *mipC = MipC_create( inst, -1, tpdSum );

    MipC_setInitialSolution( mipC, sol );

    MipC_solve( mipC );

    MipC_writeLP( mipC, "pspmip" );

    Sol_free( &sol );

    MipC_free( &mipC );

    Inst_free( &inst );
}


