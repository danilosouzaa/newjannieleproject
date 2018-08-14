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
#include "cut_pool.h"

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
    for(int s = 0 ; s < nJobs; s++)
            free(namesvar[s]);
    free(namesvar);
    free(value);


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


int main( int argc, char **argv )
{

    double startT = omp_get_wtime();


    if ( argc < 6 ) {
        fprintf( stderr, "ERROR! Enter path, instance, file lp,  time, tpd, tms, initial solution.");
        fflush( stderr );
        exit( EXIT_FAILURE );
    }


    printf("instance: %s\n",argv[2]);
    printf("lp: %s\n",argv[3]);
    printf("initial sol: %s\n",argv[7]);
    fflush(stdout);
    fflush(stderr);

    int seed = 100000;
    int sumtpd = atoi( argv[5] );
    int sumtms = atoi( argv[6] );
    int tmp = atoi( argv[4] );
    srand(seed);
    const Instance *inst = Inst_read( argv, argc );
    Solution *sol = Sol_create(inst);
    Parameters *par = Par_create(inst,argv,argc);

    //Inst_print(inst);
   // printf("nJobs out %d", Inst_nJobs(inst));// getchar();

    LinearProgram *lp = lp_create();

    lp_read(lp,argv[3]);


    int ncols = lp_cols(lp);
   // int priorities[ncols];
    //lp_fix_init_sol(lp,initSol);

    if(atoi(argv[7])!=-1) {
        Solution *initSol = Sol_create(inst);
        int read = Sol_read(initSol,argv[7]);
      //  Sol_print(initSol);
        if(read==1) {
            lp_insert_init_sol (lp,initSol);
        }
        Sol_free(&initSol);
    }

    int status=-1;

#ifdef GRB
    DataGRBCallback *usrdata = DataGRB_create( inst,tmp, Par_getLifting(par), Par_getMaxCut(par));
    lp_grb_cut_getcallback(lp, usrdata);
#endif // GRB

    double cutoffvalue = sumtpd * 100000 + sumtms;

    lp_add_cutoff( lp, (double)cutoffvalue, 0 );

    lp_set_cuts(lp,True);

    //lp_set_parallel(lp,False);
   // lp_set_nThreads(lp,1);
    lp_set_integer(lp, lp_cols(lp), lp_original_colummns(lp) );
    lp_set_max_seconds(lp,tmp);
    lp_set_print_messages(lp,1);

    double createCgraphTime = omp_get_wtime();
    _time = ( (double) tmp - (omp_get_wtime()-startT) );
    CGraph *cgraph;
    if(cc || co || cg ||cgcpu || cggpur2 || cggpu)
        cgraph = CutP_compute_conflicts_create_complete_graph(lp, origCols, lp, inst,_time, Par_getContinuous(par), Par_getMaxCut(par));
    _time = ( (double) tmp - (omp_get_wtime()-startT) );


   //                                                                                                                                                                                                                                                                                                                                                                                                                        lp_set_numeric_focus(lp,3);
 //   lp_set_mip_emphasis( lp, LP_ME_OPTIMALITY );//  LP_ME_FEASIBILITY
    //lp_config_grb_params(lp);
    status = lp_optimize(lp);


    char *namelp;
    ALLOCATE_VECTOR_INI(namelp,char,256);
    strcat(namelp,argv[2]);
    strcat(namelp,"_cut_lp.lp");
    printf("%s", namelp);

    if(atoi(argv[7])==4) {
        strcat(namelp,"_cut_and_cut_gurobi_lp.lp");
        printf("%s", namelp);

    }
    lp_write_lp(lp,namelp);

    double best_bound_root = lp_best_bound(lp);
    double time_best_bound_root =  (double) omp_get_wtime()-startT;
    double obj_value_root = lp_obj_value(lp);

    char prefix[256];
    int idx[3];

    if (status == LP_OPTIMAL || status == LP_FEASIBLE) {

        char *namesol;
        ALLOCATE_VECTOR_INI(namesol,char,256);
        strcat(namesol,argv[2]);
        strcat(namesol,"_cut_lp.sol");
        printf("%s", namesol);

        if(atoi(argv[7])==4) {
            strcat(namesol,"_cut_and_cut_gurobi_lp.sol");
            printf("%s", namesol);
        }
        lp_write_sol(lp,namesol);

        const double *x = lp_x(lp);
        for ( int i=0 ; (i<lp_cols(lp)) ; ++i ) {
            if ( fabs(x[i]) < EPS )  // not an active binary var
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
            } // x variables
        }//all variables
        Sol_calcCost(sol);
        printf("LP_STATUS %d (OPT%d,FEA%d). built sol tpd: %ld tms: %ld\n", status, LP_OPTIMAL, LP_FEASIBLE, Sol_getTPD(sol), Sol_getTMS(sol) ) ;


        char *namesolution;
        ALLOCATE_VECTOR_INI(namesolution,char,256);
        strcat(namesolution,argv[2]);
        if(status == LP_OPTIMAL){
            strcat(namesolution,"_opt_cut.sol");
        }
        if(status == LP_FEASIBLE){
            strcat(namesolution,"_fea_cut.sol");
        }

        printf("%s", namesolution);
        Sol_write(sol,namesolution);

        free(namesol);
        free(namesolution);
    }

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


    fprintf(fp, " %s ; %f ; %f ; %f ; %f ; %f ; %f ; %d ; %f ; %f ; %ld ; %ld  \n", argv[2], obj_value_root, lp_obj_value(lp), best_bound_root, lp_best_bound(lp), lp_get_gap(lp), lp_get_number_nodes_explored(lp), DataGRB_getNCuts(usrdata), time_best_bound_root, (double) omp_get_wtime()-startT,  Sol_getTPD(sol), Sol_getTMS(sol));


    free(namelp);

    fclose( fp );
    lp_free(&lp);
    lp_close_env();
    free(lp);
    Sol_free(&sol);
    Par_free(&par);
//    Res_free(&res);
    Inst_free( &inst );
    free(usrdata);



    return 0;
}
