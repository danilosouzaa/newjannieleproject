/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include "instance.h"
#include "solution.h"
#include "macros.h"
#include "time.h"
#include "lp.h"
#include "vec_str.h"


void Read_fileoutput( const char *dir, const char *file )
{
    char fname[FILE_NAME_SIZE];
    sprintf( fname, "%s/%s", dir, file );
    FILE *f = fopen( fname, "r" );

    VecStr *linesInst = VStr_create( STR_SIZE );
    VStr_readFrom( linesInst, fname, False );

    int linesFile = VStr_size(linesInst);

    char line[LINE_SIZE];
    char *s;
    int cont =0;
    while ( (s=fgets( line, LINE_SIZE, f)) ){
        if(cont == linesFile-8 )
            printf("%s %s\n", file, s);
        cont++;
    }

    VStr_free( &linesInst );

}


void lp_insert_init_sol(LinearProgram *lp, Solution *initsol)
{

    int nJobs = Inst_nJobs(Sol_inst(initsol));
    char **namesvar;
    ALLOCATE_VECTOR_INI(namesvar,char*,nJobs);
    double *value;
    ALLOCATE_VECTOR_INI(value,double*,nJobs);
    for(int s = 0 ; s < nJobs; s++) {
        ALLOCATE_VECTOR_INI(namesvar[s],char,256);
        int j = Sol_getSequence(initsol,s);
        int m = Sol_getMode(initsol,j);
        int t = Sol_getStartTime(initsol,j);
        sprintf(namesvar[s],"x(%d,%d,%d)",j,m,t);
        // printf("%s\n", namesvar[s]);
        value[s] = 1.0;
    }

    lp_load_mip_start(lp,nJobs,namesvar,value);
}


void lp_fix_init_sol(LinearProgram *lp, Solution *initsol)
{

    int nJobs = Inst_nJobs(Sol_inst(initsol));

    for(int s = 0 ; s < nJobs; s++) {
        char *namesvar;
        ALLOCATE_VECTOR_INI(namesvar,char,256);
        int j = Sol_getSequence(initsol,s);
        int m = Sol_getMode(initsol,j);
        int t = Sol_getStartTime(initsol,j);
        sprintf(namesvar,"x(%d,%d,%d)",j,m,t);
        // printf("%s\n", namesvar);
        int c = lp_col_index(lp,namesvar);
        lp_fix_col(lp,c,1.0);
        free(namesvar);
    }
}


int parseName( const char *name, char *prefix, int *idx )
{
#define MAX_COMMAS 64
    int nCommas=0;
    int commaPos[MAX_COMMAS];
    /* open and close  () */
    int pO=-1, pC=-1;
    int l = strlen(name);
    int i;
    for ( i=0 ; (i<l) ; ++i )
        switch (name[i]) {
            case '(' :
                pO = i;
                break;
            case ')' :
                pC = i;
                break;
            case ',' :
                commaPos[nCommas++] = i;
                break;
        }

    /* not in the propper format */
    if ( (pO==-1) || (pC==-1) )
        return 0;


    strcpy( prefix, name );

    //printf("%s\n", name );

    prefix[pO] = '\0';

    char str[STR_SIZE];

    for ( i=0 ; (i<nCommas+1) ; ++i ) {
        int pStart = pO;
        if (i>=1)
            pStart = commaPos[i-1];

        int pEnd   = pC;
        if (i<nCommas)
            pEnd = commaPos[i];


        //printf("X%d %d %d\n", pStart, pEnd, i);

        strcpy( str, name+pStart+1 );
        str[pEnd-pStart-1] = '\0';
        idx[i] = atoi(str);
    }

    return nCommas+1;
#undef MAX_COMMAS
}

int main( int argc, char **argv )
{

    double startT = omp_get_wtime();

    if(argc < 4) {
        Read_fileoutput(argv[1], argv[2]);
        exit(0);
    }

    if ( argc < 6 ) {
        fprintf( stderr, "ERROR! Enter path, instance, file lp, tpd, time, initial solution.");
        fflush( stderr );
        exit( EXIT_FAILURE );
    }


    printf("instance: %s\n",argv[2]);
    printf("lp: %s\n",argv[3]);
    printf("initial sol: %s\n",argv[6]);
    fflush(stdout);
    fflush(stderr);

    int seed = 100000;
    int sumtpd = atoi( argv[5] );
    int tmp = atoi( argv[4] );
    double tmpd = atof(argv[4]);
    tmp = (double) (tmpd-tmp) >= 0.5 ? ceil(tmpd) : tmp;
    srand(seed);
    Instance *inst = Inst_read( argv, argc );
    Solution *sol = Sol_create(inst);

    //Inst_print(inst);
    //getchar();

    LinearProgram *lp = lp_create();
    lp_read(lp,argv[3]);
    char prefix[256];
    int idx[3];
    int ncols = lp_cols(lp);
    int priorities[ncols];
    //lp_fix_init_sol(lp,initSol);

    if(atoi(argv[6])!=-1) {
        Solution *initSol = Sol_create(inst);
        Sol_read(initSol,argv[6]);
        Sol_print(initSol);
        lp_insert_init_sol (lp,initSol);
    }

    int status=-1;
    if(atoi(argv[7])==2) {
        for(int c = 0 ; c < ncols; c++) {
            char name[256];
            lp_col_name(lp, c, name );
            parseName( name, prefix, idx );
            if (prefix[0]=='x') {
                int j = idx[0];
                int t = idx[2];
                int m = idx[1];
                const Job*job = Inst_job(inst,j);
                priorities[c] = 100000.0 * ( 1.0/((double)Job_est(job) ) ) +  1000 * (1.0/((double)t));
                //printf("\n name %s, est %d priorities %d\n", name, Job_est(job), priorities[c]);
            } else
                priorities[c] = 1;
        }

        lp_set_branching_direction( lp, 1 );
        lp_set_branching_priorities( lp,priorities );
        //  lp_add_cutoff( lp, (double)sumtpd, 1 );

        lp_set_max_seconds(lp,tmp);
        lp_set_print_messages(lp,1);
        status = lp_optimize(lp);

    } else {
        if(atoi(argv[7])==0 || atoi(argv[7])==-1 ) {
            //solve a given lp as optimal on root with cuts and without parallism
            lp_set_max_seconds(lp,3600);
            lp_set_print_messages(lp,1);
            lp_set_max_nodes(lp,1);
            lp_set_cuts(lp,True);
            lp_set_parallel(lp,False);
            lp_set_integer(lp, lp_cols(lp), lp_original_colummns(lp) );
            status = lp_optimize(lp);
        }
         if(atoi(argv[7])==1 ) {
            //solve a given lp as optimal with cuts and without parallism
            lp_set_print_messages(lp,1);
            lp_set_cuts(lp,True);
            lp_set_parallel(lp,False);
            lp_set_integer(lp, lp_cols(lp), lp_original_colummns(lp) );
            status = lp_optimize(lp);
        }
         if(atoi(argv[7])==-3) {
            //solve a given lp as optimal on root without cuts and parallism
            lp_set_max_seconds(lp,3600);
            lp_set_print_messages(lp,1);
            lp_set_cuts(lp,False);
            lp_set_max_nodes(lp,1);
            lp_set_parallel(lp,False);
            lp_set_integer(lp, lp_cols(lp), lp_original_colummns(lp) );
            status = lp_optimize(lp);
         }
          if(atoi(argv[7])==3) {
            //solve a given lp as optimal with gurobi cuts without parallelism
            lp_set_max_seconds(lp,86395);
            lp_set_print_messages(lp,1);
            lp_set_parallel(lp,False);
            lp_set_cuts(lp,True);
            lp_set_integer(lp, lp_cols(lp), lp_original_colummns(lp) );
            status = lp_optimize(lp);
         }
         if(atoi(argv[7])==4) {
            //solve a given lp as optimal on root with cuts and without parallism
            lp_set_max_seconds(lp,3600);
            lp_set_print_messages(lp,1);
            lp_set_max_nodes(lp,1);
            lp_set_parallel(lp,False);
            lp_set_cuts(lp,True);
            lp_set_integer(lp, lp_cols(lp), lp_original_colummns(lp) );
            status = lp_optimize(lp);
         }

    }

    char *namelp;
    ALLOCATE_VECTOR_INI(namelp,char,256);
    strcat(namelp,argv[2]);


    if(atoi(argv[7])==0) {
        strcat(namelp,"_ip_lp.lp");
        printf("%s", namelp);
    }
    if(atoi(argv[7])==-1) {
        strcat(namelp,"_ip_lp_first_phase.lp");
        printf("%s", namelp);
    }
    if(atoi(argv[7])==1) {
        strcat(namelp,"_ip_lp_optimality.lp");
        printf("%s", namelp);
    }
    if(atoi(argv[7])==3) {
        strcat(namelp,"_cut_lp_optimality.lp");
        printf("%s", namelp);
    }
    if(atoi(argv[7])==-3) {
        strcat(namelp,"_cut_lp.lp");
        printf("%s", namelp);
    }
    if(atoi(argv[7])==4) {
        strcat(namelp,"_cut_and_cut_gurobi_lp.lp");
        printf("%s", namelp);

    }
    lp_write_lp(lp,namelp);
    double best_bound_root = lp_best_bound(lp);
    double time_best_bound_root =  (double) omp_get_wtime()-startT;
    double obj_value_root = lp_obj_value(lp);

    if (status == LP_OPTIMAL || status == LP_FEASIBLE) {

        char *namesol;
        ALLOCATE_VECTOR_INI(namesol,char,256);
        strcat(namesol,argv[2]);

        if(atoi(argv[7])==0) {
            strcat(namesol,"_ip_lp.sol");
            printf("%s", namesol);
        }
        if(atoi(argv[7])==1) {
            strcat(namesol,"_ip_lp_optimality.sol");
            printf("%s", namesol);
        }
        if(atoi(argv[7])==-1) {
            strcat(namesol,"_ip_lp_first_phase.sol");
            printf("%s", namesol);
        }
         if(atoi(argv[7])==3) {
            strcat(namesol,"_cut_lp_optimality.sol");
            printf("%s", namesol);
        }
        if(atoi(argv[7])==-3) {
            strcat(namesol,"_cut_lp.sol");
            printf("%s", namesol);
        }
        if(atoi(argv[7])==4) {
            strcat(namesol,"_cut_and_cut_gurobi_lp.sol");
            printf("%s", namesol);
        }
        lp_write_sol(lp,namesol);

        const double *x = lp_x(lp);
        for ( int i=0 ; (i<lp_cols(lp)) ; ++i ) {
            if ( x[i] <= 1.0-EPS )  // not an active binary var
                continue;
            char name[STR_SIZE];
            lp_col_name(lp,i,name);
            parseName( name, prefix, idx );
            if (prefix[0]=='x') {
                int j = idx[0];
                int t = idx[2];
                int m = idx[1];
                ModeSet *ms = Sol_getModeSet( sol );
                Modes_modify( ms, j, m);
                Sol_setStartJob(sol, j, t);
                Sol_calcCost(sol);
            } // x variables
        }//all variables
        printf("LP_STATUS %d (OPT%d,FEA%d). built sol tpd: %ld tms: %ld\n", status, LP_OPTIMAL, LP_FEASIBLE, Sol_getTPD(sol), Sol_getTMS(sol) ) ;


        char *namesolution;
        ALLOCATE_VECTOR_INI(namesolution,char,256);
        strcat(namesolution,argv[2]);
        if(status == LP_OPTIMAL){
            if(atoi(argv[7])==0)
                strcat(namesolution,"_opt_ip.sol");
            if(atoi(argv[7])==1)
                strcat(namesolution,"_opt_ip_optimality.sol");
            if(atoi(argv[7])==3)
                strcat(namesolution,"_opt_cut_optimality.sol");
            if(atoi(argv[7])==-1)
                strcat(namesolution,"_opt_ip_first_phase.sol");
            if(atoi(argv[7])==-3)
                strcat(namesolution,"_opt_cut.sol");
            if(atoi(argv[7])==4)
                strcat(namesolution,"_opt_cut_gurobi_cut.sol");
        }
        if(status == LP_FEASIBLE){
             if(atoi(argv[7])==0)
                strcat(namesolution,"_fea_ip.sol");
            if(atoi(argv[7])==1)
                strcat(namesolution,"_fea_ip_optimality.sol");
            if(atoi(argv[7])==3)
                strcat(namesolution,"_fea_cut_optimality.sol");
            if(atoi(argv[7])==-1)
                strcat(namesolution,"_fea_ip_first_phase.sol");
            if(atoi(argv[7])==-3)
                strcat(namesolution,"_fea_cut.sol");
            if(atoi(argv[7])==4)
                strcat(namesolution,"_fea_cut_gurobi_cut.sol");
        }

        printf("%s", namesolution);
        Sol_write(sol,namesolution);

        free(namesol);
        free(namesolution);
    } else
        printf("infeasible first phase. status %d", status);

        char fileName[300] = {"results.txt"};
    for(int n = 0; n < argc ; n++)
        if (strcmp(argv[n], "-filename") == 0) {
            strcpy(fileName, argv[n+1]);
            break;
        }

    FILE *fp = fopen( fileName, "a" );
    if ( fp == NULL ) {
        printf( "File was not opened. : path: %s\n", fileName);
        exit( 0 );
    }


    if(atoi(argv[7])==-1) {

        LinearProgram *lp_second_phase = lp_create();
        lp_read(lp_second_phase,namelp);
        int ncols = lp_cols(lp_second_phase);

        //solve a first phase lp as optimal on time with cuts and without parallism
        lp_set_max_seconds(lp_second_phase,tmp);
        lp_set_print_messages(lp_second_phase,1);
        lp_set_cuts(lp_second_phase,True);
        lp_set_max_nodes(lp_second_phase, INT_MAX_M);
        lp_set_parallel(lp_second_phase,False);
        lp_set_integer(lp_second_phase, lp_cols(lp_second_phase), lp_original_colummns(lp_second_phase) );
        status = lp_optimize(lp_second_phase);

        char *namelp;
        ALLOCATE_VECTOR_INI(namelp,char,256);
        strcat(namelp,argv[2]);
        strcat(namelp,"_ip_lp_second_phase.lp");
        printf("%s", namelp);
        lp_write_lp(lp_second_phase,namelp);
        free(namelp);

        if (status == LP_OPTIMAL || status == LP_FEASIBLE) {

            char *namesol;
            ALLOCATE_VECTOR_INI(namesol,char,256);
            strcat(namesol,argv[2]);
            strcat(namesol,"_ip_lp_second_phase.sol");
            printf("%s", namesol);
            lp_write_sol(lp_second_phase,namesol);

            const double *x = lp_x(lp_second_phase);
            for ( int i=0 ; (i<lp_cols(lp_second_phase)) ; ++i ) {
                if ( x[i] <= 1.0-EPS )  // not an active binary var
                    continue;
                char name[STR_SIZE];
                lp_col_name(lp,i,name);
                parseName( name, prefix, idx );
                if (prefix[0]=='x') {
                    int j = idx[0];
                    int t = idx[2];
                    int m = idx[1];
                    ModeSet *ms = Sol_getModeSet( sol );
                    Modes_modify( ms, j, m);
                    Sol_setStartJob(sol, j, t);
                    Sol_calcCost(sol);
                } // x variables
            }//all variables
            printf("LP_STATUS %d (OPT%d,FEA%d). built sol tpd: %ld tms: %ld\n", status, LP_OPTIMAL, LP_FEASIBLE, Sol_getTPD(sol), Sol_getTMS(sol) ) ;


            char *namesolution;
            ALLOCATE_VECTOR_INI(namesolution,char,256);
            strcat(namesolution,argv[2]);
            if(status == LP_OPTIMAL){
                strcat(namesolution,"_opt_ip_second_phase.sol");
            }
            if(status == LP_FEASIBLE){
                strcat(namesolution,"_fea_ip_second_phase.sol");
            }

            printf("%s", namesolution);
            Sol_write(sol,namesolution);

            free(namesol);
            free(namesolution);
        } else
            printf("infeasible second phase. status %d", status);

        fprintf(fp, " %s ; %f ; %f; %f; %f;  %f; %f ; %ld ; %ld  \n", argv[2], obj_value_root, lp_obj_value(lp_second_phase), best_bound_root, lp_best_bound(lp_second_phase), time_best_bound_root, (double) omp_get_wtime()-startT,  Sol_getTPD(sol), Sol_getTMS(sol));
    }else{
        fprintf(fp, " %s ; %f ; %f; %f; %f;  %f; %f ; %ld ; %ld  \n", argv[2], obj_value_root, lp_obj_value(lp), best_bound_root, lp_best_bound(lp), time_best_bound_root, (double) omp_get_wtime()-startT,  Sol_getTPD(sol), Sol_getTMS(sol));
    }

    free(namelp);

    fclose( fp );
    lp_close_env();


    Inst_free( &inst );


    return 0;
}
