
/* GOAL PSP Solver II
 * Resource-Constrained Project Scheduling Problems (RCPSP)
 *
 * Develop as part of the D.Sc. thesis of Araujo, Janniele A. S., with collaboration
 *                                   of Santos, H.G.
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <memory.h>
#include <omp.h>
#include <time.h>
#include "lp.h"
#include "cut_pool.h"
#ifdef GRB
#include "gurobi_c.h"
#endif

#include "parameters.h"
//#include "cut_gpu.h"
//#include "solutionGpu.h"
#include "prepareGpu.h"
#include "prepareCPU.h"

#include "cut_cg.h"
#include "cut_clique.h"
#include "cut_oddholes.h"
#include "cut_cover.h"
#include "cut_precedence.h"
//#include "cut_jobset.h"


#define VERBOSE 1


static const unsigned int hashval[] = { 11, 269, 3, 7, 31, 37, 131, 13, 17, 647, 653, 89, 97, 101, 39, 149, 151, 157, 821, 257, 263, 389, 397, 457, 461, 463, 331, 337, 347, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 9, 53, 59  };

static const unsigned int nHashvalues = sizeof(hashval)/sizeof(int);

int parseName( const char *name, char *prefix, int *idx )
{
#define MAX_COMMAS 64
    int nCommas=0;
    int commaPos[MAX_COMMAS];
    //open and close  ()
    int pO=-1, pC=-1;
    int l = strlen(name);
    int i;
    for ( i=0 ; (i<l) ; ++i )
        switch (name[i])
        {
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

    //not in the propper format
    if ( (pO==-1) || (pC==-1) )
        return 0;

    assert( pO<pC );

    strcpy( prefix, name );

    //printf("%s\n", name );

    prefix[pO] = '\0';

    char str[STR_SIZE];

    for ( i=0 ; (i<nCommas+1) ; ++i )
    {
        int pStart = pO;
        if (i>=1)
            pStart = commaPos[i-1];

        int pEnd   = pC;
        if (i<nCommas)
            pEnd = commaPos[i];


        //printf("X%d %d %d\n", pStart, pEnd, i);
        assert( pStart<pEnd );

        strcpy( str, name+pStart+1 );
        str[pEnd-pStart-1] = '\0';
        idx[i] = atoi(str);
    }

    return nCommas+1;
#undef MAX_COMMAS
}

int hash_code_vint( int n, const int *v, int hashSize )
{
    size_t code = 0;
    int i;
    for ( i=0 ; i<n ; ++i )
        code += ((size_t)v[i])*hashval[i%nHashvalues];

    int res = (code % ((size_t)INT_MAX));
    res = res % hashSize;

    assert( res>=0 && res<hashSize );

    return res;
}

struct _CutPool
{

    VecDbl  **hashCuts;
    VecStr  **namesCuts;
    int *nCuts;
    int nHash;
    double totalTimeSeparation;
    LinearProgram *mip;

    const Instance* inst;
};

/*create and initialize the structure of cutPool*/
CutPool *CutP_create( const Instance *inst, LinearProgram *mip)
{
    CutPool *cutP;
    ALLOCATE_INI( cutP, CutPool );
    cutP->nHash = lp_cols(mip);
    cutP->inst = inst;

    VecDbl  **hashCuts;
    ALLOCATE_VECTOR( hashCuts, VecDbl*, cutP->nHash  );
    ALLOCATE_VECTOR( cutP->namesCuts, VecStr*, cutP->nHash  );
    ALLOCATE_VECTOR_INI( cutP->nCuts, int, cutP->nHash  );

    for ( int i=0; (i<cutP->nHash) ; ++i )
    {
        hashCuts[i] = VDbl_create();
        cutP->namesCuts[i] = VStr_create(256);
    }
    cutP->hashCuts = hashCuts;
    cutP->totalTimeSeparation =0;
    cutP->mip = mip;
    return cutP;
}

/*create the conflict graph by searching conflicts between pair of activated variables also add a clique of jobs into this graph*/
CGraph *CutP_compute_conflicts_create_graph( LinearProgram *lp, const int *origCols, LinearProgram *origLP, const Instance *inst, double timeLeft)
{

    double tinit = omp_get_wtime();

    assert(lp);
    assert(inst);
    assert(origLP);

    // int nround = Res_getRound(res);

    const double *xf = lp_x(lp);
    double maxRC = lp_get_max_reduced_cost(origLP);
    const double *rdc = lp_reduced_cost(lp);
    char name[256];

    /* maximum number of indexes in a var name */
#define MAX_IDX 16
    char prefix[STR_SIZE];
    int idx[MAX_IDX];
    int maxT = 0;

    int nJobs = Inst_nJobs(inst);
    int nModes = Inst_nMaxModes(inst);

    // checking maximum time  to all variables //where there is some allocation ** if
    for ( int i=0 ; (i<lp_cols(lp)) ; ++i )
    {
        if (fabs(xf[i])>1e-5 || maxRC != -1.0)
        {
            if(maxRC != -1.0)
                if(rdc[i] > maxRC) continue;

            lp_col_name( origLP, origCols[i], name );
            if (tolower(name[0])=='x')
            {
                parseName( name, prefix, idx );
                if (prefix[0]=='x')
                {
                    int t = idx[2];
                    maxT = MAX( maxT, t );
                    //  printf( "x(%d,%d,%d) %g\n", j, m, t, xf[i] );
                }
            }
        }
    }
    int nTimes = maxT+1;
    // allocation of x
    double ***x;
    ALLOCATE_VECTOR( x, double **, nJobs );
    ALLOCATE_VECTOR( x[0], double *, nJobs*nModes );
    ALLOCATE_VECTOR( x[0][0], double, nJobs*nModes*nTimes );
    for ( int j=1 ; (j<nJobs) ; ++j )
        x[j] = x[j-1] + nModes;

    for ( int j=0; (j<nJobs) ; ++j )
        for ( int m=0 ; (m<nModes) ; ++m )
            x[j][m] = x[0][0] + j*nModes*nTimes + m*nTimes;

    // filling default values
    {
        double *it, *itEnd;
        it = x[0][0];
        itEnd = it + (nJobs*nModes*nTimes);
        for ( ; (it<itEnd) ; it++ )
            *it = 0.0;
    }

    // filling x from fractional solution
    for ( int i=0 ; (i<lp_cols(lp)) ; ++i )
    {
        if (fabs(xf[i])>1e-5 || maxRC != -1.0)
        {
            if(maxRC != -1.0)
                if(rdc[i] > maxRC) continue;

            lp_col_name( origLP, origCols[i], name );
            if (tolower(name[0]=='x'))
            {
                parseName( name, prefix, idx );
                if (prefix[0]=='x')
                    x[idx[0]][idx[1]][idx[2]] = xf[i];

            }
        }
    }

    /* x indexes pre-processed problem*/
    int ***xIdx;
    ALLOCATE_VECTOR( xIdx, int **, nJobs );
    ALLOCATE_VECTOR( xIdx[0], int *, nJobs*nModes );
    ALLOCATE_VECTOR( xIdx[0][0], int, nJobs*nModes*nTimes );
    for ( int j=1 ; (j<nJobs) ; ++j )
        xIdx[j] = xIdx[j-1] + nModes;

    for ( int j=0; (j<nJobs) ; ++j )
        for ( int m=0 ; (m<nModes) ; ++m )
            xIdx[j][m] = xIdx[0][0] + j*nModes*nTimes + m*nTimes;
    // filling default values
    {
        int *it, *itEnd;
        it = xIdx[0][0];
        itEnd = it + (nJobs*nModes*nTimes);
        for ( ; (it<itEnd) ; it++ )
            *it = -1;
    }

    /* filling indexes from variables in this pre-processed problem */
    for ( int i=0 ; (i<lp_cols(lp)) ; ++i )
    {
        if (fabs(xf[i])>1e-5 || maxRC != -1.0)
        {
            if(maxRC != -1.0)
                if(rdc[i] > maxRC) continue;

            lp_col_name( origLP, origCols[i], name );
            if (tolower(name[0])=='x')
            {
                parseName( name, prefix, idx );
                if (prefix[0]=='x')
                {
                    int jlp = idx[0];
                    int mlp =idx[1];
                    int tlp = idx[2];
                    xIdx[jlp][mlp][tlp] = i;
                }
            }
        }
    }

    int ncols = lp_cols(lp);

    int **conf;
    ALLOCATE_VECTOR(conf, int*, ncols);

    int *nConf;
    ALLOCATE_VECTOR_INI(nConf, int, ncols);

    int **varInConfs;
    ALLOCATE_VECTOR_INI(varInConfs, int*, ncols);


    int nProj = Inst_nProjects(inst);
    int nResources = Inst_nResR(inst);

    int **vectorJobs;
    ALLOCATE_VECTOR(vectorJobs,int*,nJobs);
    for(int i = 0 ; i < nJobs ; i++)
        ALLOCATE_VECTOR_INI(vectorJobs[i],int,lp_cols(lp));

    IntTriple **vectorProjects;
    ALLOCATE_VECTOR(vectorProjects,IntTriple*,nProj);
    for(int i = 0 ; i < nProj ; i++)
        ALLOCATE_VECTOR_INI(vectorProjects[i],IntTriple,ncols);

    IntTriple ***vectorTimes;
    ALLOCATE_VECTOR(vectorTimes,IntTriple**,nTimes);
    for(int t = 0 ; t < nTimes ; t++)
    {
        ALLOCATE_VECTOR(vectorTimes[t],IntTriple*,nResources);
        for(int r = 0 ; r < nResources ; r++)
            ALLOCATE_VECTOR(vectorTimes[t][r],IntTriple,ncols);
    }

    int *nVJobs;
    ALLOCATE_VECTOR_INI(nVJobs,int,nJobs);
    int *nVProj;
    ALLOCATE_VECTOR_INI(nVProj,int,nProj);
    int **nVTime;
    ALLOCATE_VECTOR(nVTime,int*,nTimes);
    for(int t = 0 ; t < nTimes; t++)
        ALLOCATE_VECTOR_INI(nVTime[t],int,nResources);

    VecInt** conflicts;
    ALLOCATE_VECTOR(conflicts,VecInt*,ncols);

    double _time;
    for(int l = 0; l < ncols ; l++)
    {
        ALLOCATE_VECTOR_INI(varInConfs[l], int, ncols);
        conflicts[l] = VInt_create();
        if (fabs(xf[l])>1e-5 || maxRC != -1.0)
        {
            if(maxRC != -1.0)
                if(rdc[l] > maxRC) continue;

            lp_col_name( origLP, origCols[l], name );
            if (tolower(name[0])=='x')
            {
                parseName( name, prefix, idx );
                if (prefix[0]=='x')
                {

                    int jlp = idx[0];
                    int mlp =idx[1];
                    int tlp = idx[2];

                    const Job *job = Inst_job(inst,jlp);
                    const Mode *mode = Job_mode(job,mlp);
                    if(Mode_duration(mode)==0)continue;

                    int plp = Job_project(job);

                    if(xIdx[jlp][mlp][tlp] == -1) continue;

                    vectorJobs[jlp][nVJobs[jlp]] = l;

                    for(int r = 0; r < nResources ; r++)
                    {

                        int idxROnInstD = Mode_idxResROnMode(inst,mode,r);
                        if( idxROnInstD == -1) continue;
                        int size = nVTime[tlp][r];
                        vectorTimes[tlp][r][size].idx = l;
                        vectorTimes[tlp][r][size].j = jlp;
                        vectorTimes[tlp][r][size].m = mlp;
                        vectorTimes[tlp][r][size].t = tlp;
                        vectorTimes[tlp][r][size].value = xf[l];

                        nVTime[tlp][r]++;

                    }

                    vectorProjects[plp][nVProj[plp]].idx = l;
                    vectorProjects[plp][nVProj[plp]].j = jlp;
                    vectorProjects[plp][nVProj[plp]].m = mlp;
                    vectorProjects[plp][nVProj[plp]].t = tlp;
                    vectorProjects[plp][nVProj[plp]].value = xf[l];

                    //   printf("(%d, %d, %d) %d | ", jlp,mlp,tlp,vectorJobs[jlp][nVJobs[jlp]]  );
                    nVJobs[jlp]++;
                    nVProj[plp]++;

                }
            }
        }
    }

    for(int t = 0 ; t < nTimes ; t++)
    {
        for(int r = 0 ; r < nResources ; r++)
        {
            if(nVTime[t][r] == 0) continue;
            for(int elem = 0 ; elem < nVTime[t][r] ; elem++)
            {
                int j = vectorTimes[t][r][elem].j, m = vectorTimes[t][r][elem].m;// t1 = vectorTimes[t][r][elem].t;
                int l = vectorTimes[t][r][elem].idx;
                const Job *job = Inst_job(inst,j);
                const Mode *mode = Job_mode(job,m);
                if(Mode_duration(mode)==0)continue;
                //   printf("\nE(%d,%d,%d)", j, m, t1);
                for(int elem2 = 0 ; elem2 < nVTime[t][r] ; elem2++)
                {
                    if(elem == elem2 ) continue;
                    int j2 = vectorTimes[t][r][elem2].j, m2 = vectorTimes[t][r][elem2].m,  t2 = vectorTimes[t][r][elem2].t;
                    if(j==j2 && m==m2 && t ==t2) continue;
                    const Job *job2 = Inst_job(inst,j2);
                    const Mode *mode2 = Job_mode(job2,m2);
                    if(Mode_duration(mode2)==0)continue;

                    int res = Mode_idxResROnMode(inst,mode,r);
                    int res2 = Mode_idxResROnMode(inst,mode2,r);
                    int sumUse = Mode_useResR(mode,res)+ Mode_useResR(mode2,res2);
                    if( sumUse > Inst_capResR(inst,r))
                    {
                        int l2 = vectorTimes[t][r][elem2].idx;
                        if(varInConfs[l][l2]==1) continue;
                        //         printf("R (%d,%d,%d) %d ", j2, m2, t2, l2);
                        VInt_pushBack(conflicts[l],l2);
                        varInConfs[l][l2] = 1;
                    }
                }
            }
        }
    }

    for(int p = 0 ; p < nProj ; p++)
    {
        for(int elem = 0 ; elem < nVProj[p] ; elem++)
        {
            int l = vectorProjects[p][elem].idx;
            int j = vectorProjects[p][elem].j, m = vectorProjects[p][elem].m, t =  vectorProjects[p][elem].t;
            const Job *job = Inst_job(inst,j);
            const Mode *mode = Job_mode(job,m);
            if(Mode_duration(mode)==0)continue;
            //    printf("\nE(%d,%d,%d)", j, m, t);
            for(int elem2 = 0 ; elem2 < nVProj[p] ; elem2++)
            {
                if(elem == elem2 ) continue;
                int j2 = vectorProjects[p][elem2].j, m2 = vectorProjects[p][elem2].m, t2 = vectorProjects[p][elem2].t;
                if(j==j2 && m==m2 && t ==t2) continue;

                const Job *job2 = Inst_job(inst,j2);
                const Mode *mode2 = Job_mode(job2,m2);
                if(Mode_duration(mode2)==0)continue;

                int durationE = Mode_duration(mode);
                int timeEndPred = t+durationE;

                int over = 0;
                if( Job_hasIndSucc(job,j2) )
                {
                    int value = Inst_getMaxDIJ(inst,j,j2)- Job_minDuration(job);
                    int winTime = value ;
                    if(timeEndPred>t2 || ((t2-timeEndPred) < winTime))
                    {
                        //     printf("timeEndPred = t %d + durationE %d \n", t, durationE );
                        //   printf("\ntimeEndPred %d > t2 %d || ((t2 %d - timeEndPred %d) < winTime %d)", timeEndPred,t2,t2,timeEndPred, winTime);
                        int l2 = vectorProjects[p][elem2].idx;
                        if(varInConfs[l][l2]==1) continue;
                        //    printf("timed D(%d,%d,%d) %d ", j2, m2, t2, l2);
                        VInt_pushBack(conflicts[l],l2);
                        varInConfs[l][l2] = 1;
                        continue;
                    }

                    int part1 = 0, cp1 = 0;
                    const Project *project = Inst_project(inst,p);
                    cp1 = Project_releaseDate(project)+Project_criticalPath(project) -  Inst_getMaxDIMJM(inst, job,Mode_index(mode),j2,Mode_index(mode2));
                    part1 = t - cp1 < 0 ? 0 :
                            t - cp1 ;
                    over =  part1;
                    if(over > Inst_getSumTPD(inst))
                    {
                        int l2 = vectorProjects[p][elem2].idx;
                        if(varInConfs[l][l2]==1) continue;
                        //   printf("tpd(%d,%d,%d) %d ", j2, m2, t2, l2);
                        VInt_pushBack(conflicts[l],l2);
                        varInConfs[l][l2] = 1;
                    }
                }
            }

            for(int p2 = 0 ; p2 < nProj ; p2++)
            {
                if(p==p2) continue;
                for(int elem2 = 0 ; elem2 < nVProj[p2] ; elem2++)
                {
                    int j2 = vectorProjects[p2][elem2].j, m2 = vectorProjects[p2][elem2].m, t2 = vectorProjects[p2][elem2].t;
                    if(j==j2 && m==m2 && t ==t2) continue;

                    const Job *job2 = Inst_job(inst,j2);
                    const Mode *mode2 = Job_mode(job2,m2);
                    if(Mode_duration(mode2)==0)continue;

                    const Project *project = Inst_project(inst,p);
                    int cp1 = Project_releaseDate(project)+Project_criticalPath(project) - Inst_getMaxDIM(inst,job,Mode_index(mode));
                    const Project *project2 = Inst_project(inst,Job_project(job2));
                    int cp2 = Project_releaseDate(project2)+Project_criticalPath(project2) - Inst_getMaxDIM(inst,job2,Mode_index(mode2));

                    int part1 = t - cp1 < 0 ? 0 :
                                t - cp1 ;
                    int part2 = t2 - cp2 < 0 ? 0 :
                                t2 - cp2 ;
                    if(part1 + part2 > Inst_getSumTPD(inst))
                    {
                        // printf("\npart1 %d (t%d-cp1 %d) + part2 %d (t2%d-cp2 %d) > Inst_getSumTPD(inst) %d  releaseDatep1 %d + criticalPathp1 %d  DJMj1 %d releaseDatep2 %d + criticalPathp %d DIM2 %d", part1, t, cp1, part2, t2, cp2, Inst_getSumTPD(inst), Project_releaseDate(project),Project_criticalPath(project) , Inst_getMaxDIM(inst,job,Mode_index(mode)), Project_releaseDate(project2),Project_criticalPath(project2),Inst_getMaxDIM(inst,job2,Mode_index(mode2)));
                        int l2 = vectorProjects[p2][elem2].idx;
                        if(varInConfs[l][l2]==1) continue;
                        //  printf("tpd p(%d,%d,%d) %d ", j2, m2, t2, l2);
                        VInt_pushBack(conflicts[l],l2); //if the sum of ends is bigger then tpd
                        varInConfs[l][l2] = 1;
                    }
                }
            }
        }
    }
    //    int weight[lp_cols(lp)];

    for(int l = 0 ; l < ncols ; l++)
    {
        int sizeConflictsJ = VInt_size(conflicts[l]);
        if(sizeConflictsJ != 0)
        {
            nConf[l] = sizeConflictsJ;
            int *confJ = VInt_getPtr(conflicts[l]);
            ALLOCATE_VECTOR(conf[l], int, sizeConflictsJ);
            for(int o = 0 ; o <  sizeConflictsJ; o++)
                conf[l][o]= confJ[o];
            // if(xf[l]<=0.00001)
            //   weight[l] = 1;
            //else
            //  weight[l] = xf[l]*1000.0;

        }//else
        //weight[l] = 0.0;
    }




    /* for(int c = 0; c < lp_cols(lp); c++) {
         if(nConf[c] > 0) {
             Res_setNSumAllVarWithConf(res,LPC_CLIQUE,1);
             Res_setNSumAllElementsConf(res,LPC_CLIQUE,nConf[c]);
             if(nConf[c] > Res_getNMaxElementsConf(res,LPC_CLIQUE))
                 Res_setNMaxElementsConf(res,LPC_CLIQUE,nConf[c]);
         }
     }
     */

    /*  printf("Conflicts\n");
       for(int ci = 0 ; ci < ncols ;ci++){
          printf("\nE %d\n", ci);
          int size = VInt_size(conflicts[ci]);
          printf("D: ");
          for(int cie = 0 ; cie < size ; cie++){
              printf(" %d ", VInt_get(conflicts[ci], cie));
          }
       }*/ // getchar();


    _time = ( (double)  (omp_get_wtime()-tinit) );
    if(VERBOSE==3) printf("\ntime finding conflicts: %f. ", _time);

    _time = ( (double) timeLeft - (double) (omp_get_wtime()-tinit) );
    CGraph *cgraph = build_cgraph_conflicts(conf, nConf, lp_cols(lp),_time);

    char recomputeDegree = 0;
    for(int j = 0 ; j < nJobs ; j++)
    {
        if(nVJobs[j]==0) continue;
        cgraph_add_clique(cgraph,vectorJobs[j],nVJobs[j]);
        recomputeDegree = 1;
    }

    if(recomputeDegree)
        cgraph_recompute_degree(cgraph);


    for(int t = 0 ; t < nTimes; t++)
    {
        free(nVTime[t]);
        for(int r = 0 ; r < nResources; r++)
            free(vectorTimes[t][r]);
        free(vectorTimes[t]);
    }
    for(int j = 0 ; j < nJobs ; j++)
        free(vectorJobs[j]);

    for(int p = 0 ; p < nProj ; p++)
        free(vectorProjects[p]);
    free(vectorJobs);
    free(vectorProjects);
    free(vectorTimes);
    free(nVJobs);
    free(nVProj);
    free(nVTime);

    for(int l = 0 ; l < lp_cols(lp) ; l++)
    {
        free(varInConfs[l]);
        VInt_free(&conflicts[l]);
        if(nConf[l] == 0) continue;
        free(conf[l]);
    }

    free(conflicts);
    free(varInConfs);
    free(conf);
    free(nConf);

    free(xIdx[0][0]);
    free(xIdx[0]);
    free(xIdx);
    free(x[0][0]);
    free(x[0]);
    free(x);

    return cgraph;
#undef MAX_IDX
}


/*create the conflict graph by searching conflicts between pair of variables all variables, also add a clique of jobs into this graph*/
CGraph *CutP_compute_conflicts_create_complete_graph( LinearProgram *lp, const int *origCols, LinearProgram *origLP, const Instance *inst, double timeLeft)
{

    double tinit = omp_get_wtime();

    assert(lp);
    assert(inst);
    assert(origLP);

    // int nround = Res_getRound(res);

    const double *xf = lp_x(lp);
    char name[256];

    /* maximum number of indexes in a var name */
#define MAX_IDX 16
    char prefix[STR_SIZE];
    int idx[MAX_IDX];
    int maxT = 0;

    int nJobs = Inst_nJobs(inst);
    int nModes = Inst_nMaxModes(inst);

    // checking maximum time  to all variables //where there is some allocation ** if
    for ( int i=0 ; (i<lp_cols(lp)) ; ++i )
    {
           lp_col_name( origLP, origCols[i], name );
            if (tolower(name[0])=='x')
            {
                parseName( name, prefix, idx );
                if (prefix[0]=='x')
                {
                    int t = idx[2];
                    maxT = MAX( maxT, t );
                    //  printf( "x(%d,%d,%d) %g\n", j, m, t, xf[i] );
                }
            }
    }
    int nTimes = maxT+1;
    // allocation of x
    double ***x;
    ALLOCATE_VECTOR( x, double **, nJobs );
    ALLOCATE_VECTOR( x[0], double *, nJobs*nModes );
    ALLOCATE_VECTOR( x[0][0], double, nJobs*nModes*nTimes );
    for ( int j=1 ; (j<nJobs) ; ++j )
        x[j] = x[j-1] + nModes;

    for ( int j=0; (j<nJobs) ; ++j )
        for ( int m=0 ; (m<nModes) ; ++m )
            x[j][m] = x[0][0] + j*nModes*nTimes + m*nTimes;

    // filling default values
    {
        double *it, *itEnd;
        it = x[0][0];
        itEnd = it + (nJobs*nModes*nTimes);
        for ( ; (it<itEnd) ; it++ )
            *it = 0.0;
    }

    // filling x from fractional solution
    for ( int i=0 ; (i<lp_cols(lp)) ; ++i )
    {
            lp_col_name( origLP, origCols[i], name );
            if (tolower(name[0]=='x'))
            {
                parseName( name, prefix, idx );
                if (prefix[0]=='x')
                    x[idx[0]][idx[1]][idx[2]] = xf[i];

            }
    }

    /* x indexes pre-processed problem*/
    int ***xIdx;
    ALLOCATE_VECTOR( xIdx, int **, nJobs );
    ALLOCATE_VECTOR( xIdx[0], int *, nJobs*nModes );
    ALLOCATE_VECTOR( xIdx[0][0], int, nJobs*nModes*nTimes );
    for ( int j=1 ; (j<nJobs) ; ++j )
        xIdx[j] = xIdx[j-1] + nModes;

    for ( int j=0; (j<nJobs) ; ++j )
        for ( int m=0 ; (m<nModes) ; ++m )
            xIdx[j][m] = xIdx[0][0] + j*nModes*nTimes + m*nTimes;
    // filling default values
    {
        int *it, *itEnd;
        it = xIdx[0][0];
        itEnd = it + (nJobs*nModes*nTimes);
        for ( ; (it<itEnd) ; it++ )
            *it = -1;
    }

    /* filling indexes from variables in this pre-processed problem */
    for ( int i=0 ; (i<lp_cols(lp)) ; ++i )
    {
            lp_col_name( origLP, origCols[i], name );
            if (tolower(name[0])=='x')
            {
                parseName( name, prefix, idx );
                if (prefix[0]=='x')
                {
                    int jlp = idx[0];
                    int mlp =idx[1];
                    int tlp = idx[2];
                    xIdx[jlp][mlp][tlp] = i;
                }
            }
    }

    int ncols = lp_cols(lp);

    int **conf;
    ALLOCATE_VECTOR(conf, int*, ncols);

    int *nConf;
    ALLOCATE_VECTOR_INI(nConf, int, ncols);

    int **varInConfs;
    ALLOCATE_VECTOR_INI(varInConfs, int*, ncols);


    int nProj = Inst_nProjects(inst);
    int nResources = Inst_nResR(inst);

    int **vectorJobs;
    ALLOCATE_VECTOR(vectorJobs,int*,nJobs);
    for(int i = 0 ; i < nJobs ; i++)
        ALLOCATE_VECTOR_INI(vectorJobs[i],int,lp_cols(lp));

    IntTriple **vectorProjects;
    ALLOCATE_VECTOR(vectorProjects,IntTriple*,nProj);
    for(int i = 0 ; i < nProj ; i++)
        ALLOCATE_VECTOR_INI(vectorProjects[i],IntTriple,ncols);

    IntTriple ***vectorTimes;
    ALLOCATE_VECTOR(vectorTimes,IntTriple**,nTimes);
    for(int t = 0 ; t < nTimes ; t++)
    {
        ALLOCATE_VECTOR(vectorTimes[t],IntTriple*,nResources);
        for(int r = 0 ; r < nResources ; r++)
            ALLOCATE_VECTOR(vectorTimes[t][r],IntTriple,ncols);
    }

    int *nVJobs;
    ALLOCATE_VECTOR_INI(nVJobs,int,nJobs);
    int *nVProj;
    ALLOCATE_VECTOR_INI(nVProj,int,nProj);
    int **nVTime;
    ALLOCATE_VECTOR(nVTime,int*,nTimes);
    for(int t = 0 ; t < nTimes; t++)
        ALLOCATE_VECTOR_INI(nVTime[t],int,nResources);

    VecInt** conflicts;
    ALLOCATE_VECTOR(conflicts,VecInt*,ncols);

    double _time;
    for(int l = 0; l < ncols ; l++)
    {
        ALLOCATE_VECTOR_INI(varInConfs[l], int, ncols);
        conflicts[l] = VInt_create();
            lp_col_name( origLP, origCols[l], name );
            if (tolower(name[0])=='x')
            {
                parseName( name, prefix, idx );
                if (prefix[0]=='x')
                {

                    int jlp = idx[0];
                    int mlp =idx[1];
                    int tlp = idx[2];

                    const Job *job = Inst_job(inst,jlp);
                    const Mode *mode = Job_mode(job,mlp);
                    if(Mode_duration(mode)==0)continue;

                    int plp = Job_project(job);

                    if(xIdx[jlp][mlp][tlp] == -1) continue;

                    vectorJobs[jlp][nVJobs[jlp]] = l;

                    for(int r = 0; r < nResources ; r++)
                    {

                        int idxROnInstD = Mode_idxResROnMode(inst,mode,r);
                        if( idxROnInstD == -1) continue;
                        int size = nVTime[tlp][r];
                        vectorTimes[tlp][r][size].idx = l;
                        vectorTimes[tlp][r][size].j = jlp;
                        vectorTimes[tlp][r][size].m = mlp;
                        vectorTimes[tlp][r][size].t = tlp;
                        vectorTimes[tlp][r][size].value = xf[l];

                        nVTime[tlp][r]++;

                    }

                    vectorProjects[plp][nVProj[plp]].idx = l;
                    vectorProjects[plp][nVProj[plp]].j = jlp;
                    vectorProjects[plp][nVProj[plp]].m = mlp;
                    vectorProjects[plp][nVProj[plp]].t = tlp;
                    vectorProjects[plp][nVProj[plp]].value = xf[l];

                    //   printf("(%d, %d, %d) %d | ", jlp,mlp,tlp,vectorJobs[jlp][nVJobs[jlp]]  );
                    nVJobs[jlp]++;
                    nVProj[plp]++;

                }
            }
    }

    for(int t = 0 ; t < nTimes ; t++)
    {
        for(int r = 0 ; r < nResources ; r++)
        {
            if(nVTime[t][r] == 0) continue;
            for(int elem = 0 ; elem < nVTime[t][r] ; elem++)
            {
                int j = vectorTimes[t][r][elem].j, m = vectorTimes[t][r][elem].m;// t1 = vectorTimes[t][r][elem].t;
                int l = vectorTimes[t][r][elem].idx;
                const Job *job = Inst_job(inst,j);
                const Mode *mode = Job_mode(job,m);
                if(Mode_duration(mode)==0)continue;
                //   printf("\nE(%d,%d,%d)", j, m, t1);
                for(int elem2 = 0 ; elem2 < nVTime[t][r] ; elem2++)
                {
                    if(elem == elem2 ) continue;
                    int j2 = vectorTimes[t][r][elem2].j, m2 = vectorTimes[t][r][elem2].m,  t2 = vectorTimes[t][r][elem2].t;
                    if(j==j2 && m==m2 && t ==t2) continue;
                    const Job *job2 = Inst_job(inst,j2);
                    const Mode *mode2 = Job_mode(job2,m2);
                    if(Mode_duration(mode2)==0)continue;

                    int res = Mode_idxResROnMode(inst,mode,r);
                    int res2 = Mode_idxResROnMode(inst,mode2,r);
                    int sumUse = Mode_useResR(mode,res)+ Mode_useResR(mode2,res2);
                    if( sumUse > Inst_capResR(inst,r))
                    {
                        int l2 = vectorTimes[t][r][elem2].idx;
                        if(varInConfs[l][l2]==1) continue;
                        //         printf("R (%d,%d,%d) %d ", j2, m2, t2, l2);
                        VInt_pushBack(conflicts[l],l2);
                        varInConfs[l][l2] = 1;
                    }
                }
            }
        }
    }

    for(int p = 0 ; p < nProj ; p++)
    {
        for(int elem = 0 ; elem < nVProj[p] ; elem++)
        {
            int l = vectorProjects[p][elem].idx;
            int j = vectorProjects[p][elem].j, m = vectorProjects[p][elem].m, t =  vectorProjects[p][elem].t;
            const Job *job = Inst_job(inst,j);
            const Mode *mode = Job_mode(job,m);
            if(Mode_duration(mode)==0)continue;
            //    printf("\nE(%d,%d,%d)", j, m, t);
            for(int elem2 = 0 ; elem2 < nVProj[p] ; elem2++)
            {
                if(elem == elem2 ) continue;
                int j2 = vectorProjects[p][elem2].j, m2 = vectorProjects[p][elem2].m, t2 = vectorProjects[p][elem2].t;
                if(j==j2 && m==m2 && t ==t2) continue;

                const Job *job2 = Inst_job(inst,j2);
                const Mode *mode2 = Job_mode(job2,m2);
                if(Mode_duration(mode2)==0)continue;

                int durationE = Mode_duration(mode);
                int timeEndPred = t+durationE;

                int over = 0;
                if( Job_hasIndSucc(job,j2) )
                {
                    int value = Inst_getMaxDIJ(inst,j,j2)- Job_minDuration(job);
                    int winTime = value ;
                    if(timeEndPred>t2 || ((t2-timeEndPred) < winTime))
                    {
                        //     printf("timeEndPred = t %d + durationE %d \n", t, durationE );
                        //   printf("\ntimeEndPred %d > t2 %d || ((t2 %d - timeEndPred %d) < winTime %d)", timeEndPred,t2,t2,timeEndPred, winTime);
                        int l2 = vectorProjects[p][elem2].idx;
                        if(varInConfs[l][l2]==1) continue;
                        //    printf("timed D(%d,%d,%d) %d ", j2, m2, t2, l2);
                        VInt_pushBack(conflicts[l],l2);
                        varInConfs[l][l2] = 1;
                        continue;
                    }

                    int part1 = 0, cp1 = 0;
                    const Project *project = Inst_project(inst,p);
                    cp1 = Project_releaseDate(project)+Project_criticalPath(project) -  Inst_getMaxDIMJM(inst, job,Mode_index(mode),j2,Mode_index(mode2));
                    part1 = t - cp1 < 0 ? 0 :
                            t - cp1 ;
                    over =  part1;
                    if(over > Inst_getSumTPD(inst))
                    {
                        int l2 = vectorProjects[p][elem2].idx;
                        if(varInConfs[l][l2]==1) continue;
                        //   printf("tpd(%d,%d,%d) %d ", j2, m2, t2, l2);
                        VInt_pushBack(conflicts[l],l2);
                        varInConfs[l][l2] = 1;
                    }
                }
            }

            for(int p2 = 0 ; p2 < nProj ; p2++)
            {
                if(p==p2) continue;
                for(int elem2 = 0 ; elem2 < nVProj[p2] ; elem2++)
                {
                    int j2 = vectorProjects[p2][elem2].j, m2 = vectorProjects[p2][elem2].m, t2 = vectorProjects[p2][elem2].t;
                    if(j==j2 && m==m2 && t ==t2) continue;

                    const Job *job2 = Inst_job(inst,j2);
                    const Mode *mode2 = Job_mode(job2,m2);
                    if(Mode_duration(mode2)==0)continue;

                    const Project *project = Inst_project(inst,p);
                    int cp1 = Project_releaseDate(project)+Project_criticalPath(project) - Inst_getMaxDIM(inst,job,Mode_index(mode));
                    const Project *project2 = Inst_project(inst,Job_project(job2));
                    int cp2 = Project_releaseDate(project2)+Project_criticalPath(project2) - Inst_getMaxDIM(inst,job2,Mode_index(mode2));

                    int part1 = t - cp1 < 0 ? 0 :
                                t - cp1 ;
                    int part2 = t2 - cp2 < 0 ? 0 :
                                t2 - cp2 ;
                    if(part1 + part2 > Inst_getSumTPD(inst))
                    {
                        // printf("\npart1 %d (t%d-cp1 %d) + part2 %d (t2%d-cp2 %d) > Inst_getSumTPD(inst) %d  releaseDatep1 %d + criticalPathp1 %d  DJMj1 %d releaseDatep2 %d + criticalPathp %d DIM2 %d", part1, t, cp1, part2, t2, cp2, Inst_getSumTPD(inst), Project_releaseDate(project),Project_criticalPath(project) , Inst_getMaxDIM(inst,job,Mode_index(mode)), Project_releaseDate(project2),Project_criticalPath(project2),Inst_getMaxDIM(inst,job2,Mode_index(mode2)));
                        int l2 = vectorProjects[p2][elem2].idx;
                        if(varInConfs[l][l2]==1) continue;
                        //  printf("tpd p(%d,%d,%d) %d ", j2, m2, t2, l2);
                        VInt_pushBack(conflicts[l],l2); //if the sum of ends is bigger then tpd
                        varInConfs[l][l2] = 1;
                    }
                }
            }
        }
    }

      //  int weight[lp_cols(lp)];

    for(int l = 0 ; l < ncols ; l++)
    {
        int sizeConflictsJ = VInt_size(conflicts[l]);
        if(sizeConflictsJ != 0)
        {
            nConf[l] = sizeConflictsJ;
            int *confJ = VInt_getPtr(conflicts[l]);
            ALLOCATE_VECTOR(conf[l], int, sizeConflictsJ);
            for(int o = 0 ; o <  sizeConflictsJ; o++)
                conf[l][o]= confJ[o];

            //if(xf[l]<=0.00001)
              //  weight[l] = 1;
            //else
                //weight[l] = xf[l]*1000.0;
        }// else
         //   weight[l] = 0.0;
    }

    _time = ( (double)  (omp_get_wtime()-tinit) );
    if(VERBOSE==3) printf("\ntime finding conflicts: %f. ", _time);

    _time = ( (double) timeLeft - (double) (omp_get_wtime()-tinit) );
    CGraph *cgraph = build_cgraph_conflicts(conf, nConf, lp_cols(lp),_time);
   // cgraph_set_weight( cgraph, weight );

    char recomputeDegree = 0;
    for(int j = 0 ; j < nJobs ; j++)
    {
        if(nVJobs[j]==0) continue;
        cgraph_add_clique(cgraph,vectorJobs[j],nVJobs[j]);
        recomputeDegree = 1;
    }

    if(recomputeDegree)
        cgraph_recompute_degree(cgraph);


    for(int t = 0 ; t < nTimes; t++)
    {
        free(nVTime[t]);
        for(int r = 0 ; r < nResources; r++)
            free(vectorTimes[t][r]);
        free(vectorTimes[t]);
    }
    for(int j = 0 ; j < nJobs ; j++)
        free(vectorJobs[j]);

    for(int p = 0 ; p < nProj ; p++)
        free(vectorProjects[p]);
    free(vectorJobs);
    free(vectorProjects);
    free(vectorTimes);
    free(nVJobs);
    free(nVProj);
    free(nVTime);

    for(int l = 0 ; l < lp_cols(lp) ; l++)
    {
        free(varInConfs[l]);
        VInt_free(&conflicts[l]);
        if(nConf[l] == 0) continue;
        free(conf[l]);
    }

    free(conflicts);
    free(varInConfs);
    free(conf);
    free(nConf);

    free(xIdx[0][0]);
    free(xIdx[0]);
    free(xIdx);
    free(x[0][0]);
    free(x[0]);
    free(x);

    return cgraph;
#undef MAX_IDX
}


/*find the cut into the hash and return if it does/does not exist or if it is inactive or dominated */
IntDblTuple cutP_findElement( CutPool *cutP, int key, int c, const int *idx, const double *coe, double rhs, int sense, int lp_cols, int type)
{

    IntDblTuple posrow;
    assert(key < lp_cols);
    assert(key >= 0);

    int size = VDbl_size(cutP->hashCuts[key]);
    int idRow = 0, idType = 1, idRhs = 2, idSense =3, idNElem = 4, idNCut = 5, init = 0, end = 0, endcoef = 0;
    double ty = 0.0, nElemCut = 0.0, nrhs = 0.0, nsense = 0.0;//  nCut = 0.0;//, row = 0;,
    //sense: 0 = L, 1 = E, 2 = G
    //    printf("0) > size %d, idRow %d, idType %d, idRhs %d, idSense %d, idNElem %d,  idNCut %d, ty %f, nrhs %f, nsesne %f, nelemcut %f,  ncut %f, type %d, rhs %f, sense %d, c %d, init %d, end %d, endcoef %d \n", size,idRow, idType, idRhs, idSense, idNElem, idNCut, ty, nrhs, nsense,  nElemCut, nCut, type, rhs, sense, c, init, end, endcoef);
    if(size > 0)
    {
        //   row = VInt_get(cutP->hashCuts[key], idRow);
        ty = VDbl_get(cutP->hashCuts[key], idType);
        nrhs = VDbl_get(cutP->hashCuts[key], idRhs);
        nsense = VDbl_get(cutP->hashCuts[key], idSense);
        nElemCut = VDbl_get(cutP->hashCuts[key], idNElem);
       // nCut = VDbl_get(cutP->hashCuts[key], idNCut);
        init = idNCut+1;
        end = init+nElemCut-1;
        endcoef = end+nElemCut;

        //    printf("1) > size %d, idRow %d, idType %d, idRhs %d, idSense %d, idNElem %d, idNCut %d,  ty %f, nrhs %f, nsesne %f, type %d, rhs %f, sense %d, nelemcut %f,  ncut %f, init %d, end %d, endcoef %d \n", size,idRow, idType, idRhs, idSense, idNElem, idNCut, ty, nrhs, nsense, type, rhs, sense, nElemCut, nCut,  init, end, endcoef);
    }


    int ih = init;
    int ihcoef = end+1;
    // printf("2) ih %d ihcoef %d c %d \n", ih, ihcoef, c);

    //  printf("\n3) Verify on key %d with size %d: \n if the cut exists: ", key, c);
    //  for(int ihce = 0 ; ihce < c ; ihce++)
    //      printf(" %f * %d ",coe[ihce], idx[ihce]);
    // printf("\n");
    // getchar();
    while( ih < size )
    {
        while( ( nElemCut < c) || (nElemCut > c) || (type != ty) || (sense != nsense ) || (rhs != nrhs))
        {
            //  printf("4) ih %d < size : c %d <> nElemCut %f . type %d != ty %f . sense %d != nsense %f . rhs %f != nrhs %f. \n", ih,c, nElemCut, type, ty, sense, nsense, rhs, nrhs);
            idNCut = endcoef+6;
            idNElem = endcoef+5;
            idSense = endcoef+4;
            idRhs = endcoef+3;
            idType = endcoef+2;
            idRow = endcoef+1;
            if(idNCut>=size)
            {
                posrow.a = 0;
                posrow.b = 0;
                posrow.c = 0;
                //    printf("*actv %f \n", 0);
                //          printf("5) idNElem %d >=size %d ih %d, init %d, ihcoef %d, nElemCut %f, nCut %f,  end %d, endcoef %d, idRow %d\n ", idNElem, size, ih, init, ihcoef, nElemCut, nCut, end, endcoef, idRow);
                return posrow;
            }
//            nCut = VDbl_get(cutP->hashCuts[key], idNCut);
            nElemCut = VDbl_get(cutP->hashCuts[key], idNElem);
            ty = VDbl_get(cutP->hashCuts[key], idType);
            nrhs = VDbl_get(cutP->hashCuts[key], idRhs);
            nsense = VDbl_get(cutP->hashCuts[key], idSense);
            init = idNCut+1;
            end = init+nElemCut-1;
            endcoef = end+nElemCut;
            ih = init;
            ihcoef = end+1;
            //     printf("6) idNElem %d still < size and <> %d ih %d, init %d, ihcoef %d, nElemCut %f, nCut %f, end %d, endcoef %d, idRow %d\n ", idNElem, size, ih, init, ihcoef, nElemCut, nCut, end, endcoef, idRow);
        }
        // printf("7) idNElem %d still < size and == %d ih %d, init %d, ihcoef %d, nElemCut %f, nCut %f, end %d, endcoef %d, idRow %d\n ", idNElem, size, ih, init, ihcoef, nElemCut, nCut, end, endcoef, idRow);
        //   getchar();
        int e, aux = 0;
        // printf("8) nElemCut %f \n ", nElemCut);
        for( e = 0 ; e < nElemCut ; e++ )
        {
            if(idx[e] != VDbl_get(cutP->hashCuts[key],ih) || coe[e] != VDbl_get(cutP->hashCuts[key],ihcoef) )
            {
                //     printf("!= idx[%d] %d != VDbl_get(cutP->hashCuts[%d],%d) %f || coe[%d] %f != VDbl_get(cutP->hashCuts[%d],%d) %f \n", e, idx[e],key, ih, VDbl_get(cutP->hashCuts[key],ih), e, coe[e], key, ihcoef, VDbl_get(cutP->hashCuts[key],ihcoef) );
                aux = 0;
                break;
            }
            // printf("== idx[%d] %d != VDbl_get(cutP->hashCuts[%d],%d) %f || coe[%d] %f != VDbl_get(cutP->hashCuts[%d],%d) %f \n", e, idx[e],key, ih, VDbl_get(cutP->hashCuts[key],ih), e, coe[e], key, ihcoef, VDbl_get(cutP->hashCuts[key],ihcoef) );
            ih++;
            ihcoef++;
            aux = 1;
        }
        if(aux)
        {
            double actv = VDbl_get(cutP->hashCuts[key],idRow); //new
            //      printf("actv %f------\n", actv);
            if(actv==-1)
            {
                posrow.a = idRow;
                posrow.b = -1;
                posrow.c = VDbl_get(cutP->hashCuts[key], idNCut);
                //   printf("\nactv %f", actv);
                //         printf("7) ACTV new idNElem %d >=size %d ih %d, init %d, ihcoef %d, nElemCut %f, nCut %f, end %d, endcoef %d, idRow %d\n ", idNElem, size, ih, init, ihcoef, nElemCut, nCut, end, endcoef, idRow);
            }
            else if(actv==-2)
            {
                posrow.a = idRow;
                posrow.b = -2;
                posrow.c = VDbl_get(cutP->hashCuts[key], idNCut);
                // printf("\nactv %f", actv);
                //        printf("7) ACTV new idNElem %d >=size %d ih %d, init %d, ihcoef %d, nElemCut %f, nCut %f, end %d, endcoef %d, idRow %d\n ", idNElem, size, ih, init, ihcoef, nElemCut, nCut, end, endcoef, idRow);
            }
            else
            {
                posrow.a = idRow;
                posrow.b = 1;
                posrow.c = VDbl_get(cutP->hashCuts[key], idNCut);
                // printf("actv %f\n", actv);
            }
            //   printf("8) Repeated new idNElem %d >=size %d ih %d, init %d, ihcoef %d, nElemCut %f, nCut %f, end %d, endcoef %d, idRow %d\n ", idNElem, size, ih, init, ihcoef, nElemCut, nCut,  end, endcoef, idRow);
            return posrow;
        }

        idNCut = endcoef+6;
        idNElem = endcoef+5;
        idSense = endcoef+4;
        idRhs = endcoef+3;
        idType = endcoef+2;
        idRow = endcoef+1;
        if(idNCut>=size)
        {
            posrow.a = 0;
            posrow.b = 0;
            posrow.c = 0;
            //            printf("**actv %f \n", 0);
            //    printf("9) new idNElem %d >=size %d ih %d, init %d, ihcoef %d, nElemCut %f, nCut %f, end %d, endcoef %d, idRow %d\n ", idNElem, size, ih, init, ihcoef, nElemCut, nCut, end, endcoef, idRow);
            return posrow;
        }
//        nCut = VDbl_get(cutP->hashCuts[key],idNCut);
        nElemCut = VDbl_get(cutP->hashCuts[key],idNElem);
        ty = VDbl_get(cutP->hashCuts[key], idType);
        nrhs = VDbl_get(cutP->hashCuts[key], idRhs);
        nsense = VDbl_get(cutP->hashCuts[key], idSense);
        init = idNCut+1;
        end = init+nElemCut-1;
        endcoef = end+nElemCut;
        ih = init;
        ihcoef = end+1;
    }

    posrow.a = 0;
    posrow.b = 0;
    posrow.c = 0;
    //  printf("***actv %f\n", 0);
    // printf("10) posrow.a %d, posrow.b %f, posrow.c %f\n ", posrow.a, posrow.b, posrow.c);

    return posrow;
}

/*applies max div common at the new cut and verify if the new cut is repeated.
the cuts are inserted into the hash at this moment*/
IntDblTuple CutP_compElemRepeatedCoef( CutPool *cutP, int key, int c, const int* idx, double* coe, double *rhs, int sense, int lp_cols, int lp_rows, int type) // cutPool
{


    assert(key < lp_cols);

    assert(key >= 0);

    int coemdc[c+1];

    //printf("key %d\n", key);
    for(int cc = 0; cc < c ; cc++)
    {
        coemdc[cc] = ROUND(coe[cc]);
        //printf (" %d * %d ", coemdc[cc], idx[cc]);
    }
    coemdc[c] = ROUND(*rhs);
    // printf (" RHS %d ", coemdc[c]);
    //printf (": end print original \n");

    int mdc = CutP_maxDivisorCommonVector(coemdc, c);

    //printf ("\n ------ MDC: %d ------ \n", mdc);
    for(int cc = 0; cc < c ; cc++)
    {
        coe[cc] = (double) coe[cc]/mdc;
        //  printf ("%f * %d ", coe[cc], idx[cc]);
    }
    *rhs = (double)*rhs/(double)mdc;


    // printf (" RHS %f ", *rhs);
    // printf (": end print reduced \n\n"); getchar();

    //double ttimefindelement = omp_get_wtime();
    IntDblTuple posrow = cutP_findElement(cutP,key,c,idx,coe,*rhs,sense,lp_cols,type);
    //double ttimefindelementfim = omp_get_wtime()-ttimefindelement;
    //printf("\n time cutP_findElement %f\n", ttimefindelementfim); fflush(stdout); getchar();
    // CutP_printHash(cutP);

    int dominado = 0;
    if(posrow.b==0)
    {
        // double ttimedominance = omp_get_wtime();
        //dominado = CutP_dominance(cutP, idx, coe, *rhs, sense,type,c);
        //double ttimedominancefim = omp_get_wtime()-ttimedominance;
        //printf("\n time CutP_dominance %f\n", ttimedominancefim);
        //fflush(stdout);

        if(dominado==1)
        {
            //      printf("\ndominated\n");
            posrow.b=3;
        }
        else
        {

            VDbl_pushBack(cutP->hashCuts[key], -1);
            VStr_pushBack(cutP->namesCuts[key], " ");
            posrow.a = VDbl_size(cutP->hashCuts[key])-1;
            VDbl_pushBack(cutP->hashCuts[key],type);
            VDbl_pushBack(cutP->hashCuts[key],*rhs);
            VDbl_pushBack(cutP->hashCuts[key],sense);
            VDbl_pushBack(cutP->hashCuts[key],c);
            VDbl_pushBack(cutP->hashCuts[key], cutP->nCuts[key] );
            posrow.c = cutP->nCuts[key];
            //printf("\n ultimo element %f ,size %d, row : %d cut %d \n adding cut to cut pool: \n", VDbl_get(cutP->hashCuts[key], posrow.a), c, posrow.a,posrow.c  );
            for(int ihce = 0 ; ihce < c ; ihce++)
            {
                VDbl_pushBack(cutP->hashCuts[key],idx[ihce]);
                //  printf(" %d ", idx[ihce]);
            }
            for(int ihce = 0 ; ihce < c ; ihce++)
            {
                VDbl_pushBack(cutP->hashCuts[key],coe[ihce]);
                //printf(" %f", coe[ihce]);
            }
            cutP->nCuts[key]++;
            //printf("\n"); fflush(stdout); getchar();
        }
    }
    return posrow;
}

/*print hash of cuts at CutPool*/
void CutP_printHash( CutPool *cutP)
{

    const double *xf = lp_x(cutP->mip);
    int nColsLP = lp_cols(cutP->mip);

    int idxSC = 4,  idxCut = 5;// idxTC = 1,  idxR = 0, idxRhs = 2, idxSense = 3;
    double sc = 0;// r = 0, nrhs = 0, nsense=0;// tc = 0;// ncut = 0; //row, cut, size cut , ,
    int init = 0, initcoef = 0, end = 0;// endcoef = 0;

    printf("\n\n Start \n\n");
    for(int n = 0 ; n < cutP->nHash ; n++)
    {

        int s = VDbl_size(cutP->hashCuts[n]); //size hash
        if(s == 0) continue;
      //  r = VDbl_get(cutP->hashCuts[n], idxR);
//        tc = VDbl_get(cutP->hashCuts[n], idxTC);
   //     nrhs = VDbl_get(cutP->hashCuts[n], idxRhs);
     //   nsense = VDbl_get(cutP->hashCuts[n], idxSense);
        sc = VDbl_get(cutP->hashCuts[n], idxSC);
//        ncut = VDbl_get(cutP->hashCuts[n], idxCut);
        init = idxCut+1;
        end = init+sc-1;
        initcoef = end+1;
//        endcoef = end+sc;
        printf("\nKey: %d Size: %d \n", n, s);

        // printf("1) idxR %d, r %f, idxTC %d, tc %f, idxRhs %d, nrhs %f, idxSense %d,  nsense %f, idxSC %d, sc %f, init %d, initcoef %d, end %d, endcoef %d, s %d \n", idxR, r, idxTC, tc, idxRhs, nrhs, idxSense, nsense, idxSC, sc, init, initcoef, end, endcoef, s);

        int *idx;
        double *coef;
        ALLOCATE_VECTOR(idx, int, nColsLP);
        ALLOCATE_VECTOR(coef, double, nColsLP);

        int o = 0;
        double value = 0;

        int ihcoef = initcoef;
        for(int ihce = init ; ihce < init+sc && ihcoef < initcoef+sc && ihcoef < s; ihce++, ihcoef++)
        {
            // printf("\nihce %d == init %d +sc %f -1 && ihcoef %d == initcoef %d +sc %f -1 %d \n", ihce, init, sc, ihcoef,initcoef,sc);

            int elem = VDbl_get(cutP->hashCuts[n],ihce);
            double coefs = VDbl_get(cutP->hashCuts[n],ihcoef);

            value += (xf[elem]*coefs);
            idx[o] = elem;
            coef[o] = coefs;
            printf(" %f * %d ,", coefs, elem);
            o++;
            if(ihce == init+sc-1 && ihcoef == initcoef+sc-1)
            {
                /*if(sc == 1){
                   printf(" somente 1 nElem %f, sense %f, rhs %f\n", sc, nsense, nrhs );
                                getchar();
                }*/
              //  printf("nElem %f, sense %f, rhs %f r(activated) %f\n", sc, nsense, nrhs, r );
               // getchar();
                o = 0;
                value = 0;
                if(ihcoef+1 == s)
                {
//                    idxR = 0;
//                    idxTC = 1;
//                    idxRhs = 2;
  //                  idxSense = 3;
                    idxSC = 4;
                    idxCut = 5;
                    init = idxCut+1;
                    end = init+sc-1;
                    initcoef = end+1;
//                    endcoef = end+sc;
                    break;
                }
//                idxR = ihcoef+1;
//                idxTC = ihcoef+2;
//                idxRhs = ihcoef+3;
//                idxSense = ihcoef+4;
                idxSC = ihcoef+5;
                idxCut = ihcoef+6;
//                r = VDbl_get(cutP->hashCuts[n], idxR);
//                tc = VDbl_get(cutP->hashCuts[n], idxTC);
  //              nrhs = VDbl_get(cutP->hashCuts[n], idxRhs);
    //            nsense = VDbl_get(cutP->hashCuts[n], idxSense);
                sc = VDbl_get(cutP->hashCuts[n], idxSC);
//                ncut = VDbl_get(cutP->hashCuts[n], idxCut);
                init = idxCut+1;
                end = init+sc-1;
                initcoef = end+1;
             //   endcoef = end+sc;
                ihce = init;
                ihcoef = initcoef;
            }
        }
        free(coef);
        free(idx);
    }


}

/* add cuts at the lp checking if cuts are non violated, repeated, dominated, new or inactive*/
int CutP_addSeparatedCuts( CutPool *cutP, Results *res, LinearProgram* lp,  VecInt **cutsElem, VecDbl **cutsCoef, VecDbl *cutsRhs, VecDbl *cutsViolation, VecStr *cutsName, VecInt *cutsSense, VecInt *cutDominated, int continuous, int nround, int cut, double timerem) // double timeseparation, int naddcut)
{
    double ttimerepfim =0,ttimeaddmodelfim=0;
    double tini = omp_get_wtime();
    if(timerem < omp_get_wtime()-tini) return 0;
    double *xf = lp_x(lp);
    int repeated = 0, dominated = 0, maxcuts = 0, newcut = 0, inactive=0;
  //  printf(" \n%d", VInt_size(cutsSense));// getchar();
    // fflush(stdout);

    //  for(int nci = 0 ; nci < VInt_size(cutsRhs); nci++)

    int nRows = 0;
//    int *idxRows;
    char **namesRows;
//    double *coefRows;
    int *startsRows;
    char *senseRows;
    double *rhsRows;


    int contelemrows =0;
  //  ALLOCATE_VECTOR(coefRows,double, lp_cols(lp)*VInt_size(cutsSense));
    //ALLOCATE_VECTOR(idxRows, int, lp_cols(lp)*VInt_size(cutsSense));
    ALLOCATE_VECTOR_INI(startsRows,int, VInt_size(cutsSense)+1);
    ALLOCATE_VECTOR(senseRows,char, VInt_size(cutsSense));
    ALLOCATE_VECTOR_INI(rhsRows,double, VInt_size(cutsSense));
    ALLOCATE_VECTOR(namesRows,char*, VInt_size(cutsSense));
    for(int in = 0 ; in <  VInt_size(cutsSense) ; in++)
        ALLOCATE_VECTOR_INI(namesRows[in], char, 512);

    VecDbl *coefRowsVec = VDbl_create();
    VecInt *idxRowsVec = VInt_create();

    for(int nci = 0 ; nci < VInt_size(cutsSense); nci++)
    {
        if( VInt_get(cutDominated,nci) == 1 )
        {
            continue;
        }
        double vi = VDbl_get(cutsViolation,nci);

        if(vi <= 0.0002)
        {
            printf("Non Violated cut %f\n", vi); //getchar();
            continue;
        }

        int sizeCut = VInt_size(cutsElem[nci]);

        if(sizeCut ==0) continue;

        int *cutIdx = VInt_getPtr(cutsElem[nci]);

        double *cutCoef = VDbl_getPtr(cutsCoef[nci]);

        int key =  hash_code_vint(  sizeCut, cutIdx,  lp_cols(lp ));

        //  printf("Key: %d \n", key);
        //    getchar();
        int sense = VInt_get(cutsSense,nci); //G

        double rhs = VDbl_get(cutsRhs,nci); //G

        const char *namecut;

        namecut = VStr_get(cutsName,nci);

        //        int aux  = 0;
        //    printf("Antes: %s\n",namecut);
        //   fflush(stdout);
        //   printf("%d %d\n", lp_cols(lp),lp_rows(lp));
        //  fflush(stdout);
        double ttimerep = omp_get_wtime();
        IntDblTuple posrow = CutP_compElemRepeatedCoef(cutP, key, sizeCut, cutIdx, cutCoef, &rhs, sense, lp_cols(lp), lp_rows(lp), cut); //erro nessa chamada
        ttimerepfim += omp_get_wtime()-ttimerep;


        //  printf("Depois: %s\n",namecut);
        //fflush(stdout);

        // printf("CutP_compElemRepeatedCoef: end\n");fflush(stdout);
        // printf("posrow.a %d, posrow.b %f \n", posrow.a, posrow.b);
        if(posrow.b==1 || posrow.b==3 )
        {
            if(posrow.b==1 ){
                repeated++;
            //    printf("repetido\n");
            }
            else{
                dominated++;
              //  printf("dominado\n");
            }
            maxcuts++;
            continue;

        }
        else if(posrow.b == 0)
        {
            CutP_setIdRow( cutP, key, posrow.a, 1);
            CutP_setNameRow(cutP,key, posrow.c, namecut);
         //   printf("novo\n"); //getchar();
            newcut++;
        }
        else
        {
            CutP_setIdRow( cutP, key, posrow.a, 1);
            CutP_setNameRow(cutP,key, posrow.c,  namecut);
         //   printf("inativo\n");
            inactive++;
        }

        char s = ' ';
        if(sense == 0)
            s = 'L';
        else if(sense == 1)
            s = 'E';
        else if(sense == 2)
            s = 'G';

        double ttimeaddmodel = omp_get_wtime();
        nRows++;
        startsRows[nRows-1] = contelemrows;
        // printf("nRows %d, startsRows[nRows]  %d,", nRows, startsRows[nRows] );getchar();
        for(int ix = 0; ix < sizeCut ; ix++)
        {
           // idxRows[contelemrows+ix] = cutIdx[ix];
            VInt_pushBack(idxRowsVec, cutIdx[ix]);
//            coefRows[contelemrows+ix] = cutCoef[ix];
            VDbl_pushBack(coefRowsVec, cutCoef[ix]);
        }

        //COMBINE_VECTOR(int, idxRows, (sizeof(idxRows)), cutIdx, sizeof(cutIdx));
        //COMBINE_VECTOR(double, coefRows, (sizeof(coefRows)), cutCoef, sizeof(cutCoef));
        contelemrows += sizeCut;
        rhsRows[nRows-1] = rhs;
        senseRows[nRows-1] = s;
        strcpy(namesRows[nRows-1],namecut);
        //  printf("namesRows[nRows] %s, namecut %s", namesRows[nRows-1],namecut);//getchar();

        //  startsRows[nRows] = startsRows[nRows-1]+sizeCut;
//        printf(" nRows %d startsRows[nRows] %d  ",  nRows , startsRows[nRows]); getchar();


        ttimeaddmodelfim += omp_get_wtime()-ttimeaddmodel;

        // printf("Add rows cut: end\n");
        //fflush(stdout);
        //  printf("Print information: start\n");fflush(stdout);
        if(VERBOSE)
        {
            Res_setNElementsCuts(res,nround, cut,  sizeCut);
            if( sizeCut > Res_getNMaxElementsCuts(res,nround,cut) )
                Res_setNMaxElementsCuts(res,nround,cut,sizeCut);
            if( sizeCut < Res_getNMinElementsCuts(res,nround,cut) )
                Res_setNMinElementsCuts(res,nround, cut,  sizeCut);

            //   double vi = valueCutIdxCO[icc].b*-1;
            //   if(cut==LPC_CLIQUE) vi = vi-1;

            Res_setSumViol(res,nround,cut,vi);
            if(vi > Res_getMaxViol(res,nround,cut) )
            {
                Res_setMaxViol(res,nround,cut,vi);
                Res_setNElementsMaxViol(res,nround,cut, sizeCut);
                Res_setIdxStartCutsMaxViol(res,cut);
                // printf("Max Viol inside: \n");
                for(int ci  = 0 ; ci <  sizeCut ; ci++)
                {
                    Res_setMaxViolCut(res, nround, cut, ci, cutIdx[ci], cutCoef[ci]);
                    //   printf(" %f*%d : %f + ", cutIdxCo[ic][ci].b,cutIdxCo[ic][ci].a, xf[ci, cutIdxCo[ic][ci].a]);
                }
            }
            if(vi < Res_getMinViol(res,nround,cut) )
            {
                Res_setMinViol(res,nround,cut,vi);
                Res_setNElementsMinViol(res,nround,cut, sizeCut);
                Res_setIdxStartCutsMinViol(res,cut);
                //   printf("Min Viol inside: \n");
                for(int ci  = 0 ; ci <  sizeCut ; ci++)
                {
                    Res_setMinViolCut(res, nround, cut, ci, cutIdx[ci], cutCoef[ci]);
                    //    printf(" %f*%d : %f + ", cutIdxCo[ic][ci].b,cutIdxCo[ic][ci].a, xf[ci, cutIdxCo[ic][ci].a]);
                }
            }

            Res_setNCutsTotal(res,nround,cut,1);
        }
        if(timerem < omp_get_wtime()-tini+1) break;
    }
    startsRows[nRows] = contelemrows;

    if (nRows > 0)
        lp_add_rows(lp, nRows, startsRows,  VInt_getPtr(idxRowsVec), VDbl_getPtr(coefRowsVec), senseRows, rhsRows, (const char**) namesRows);
        //lp_add_rows(lp, nRows, startsRows,  idxRows, coefRows, senseRows, rhsRows, (const char**) namesRows);



    // if(VInt_size(cutsSense)>0){
   // free(idxRows);
  //  free(coefRows);
    free(startsRows);
    free(senseRows);
    free(rhsRows);
    for(int in = 0 ; in <  VInt_size(cutsSense) ; in++)
        free(namesRows[in]);
    free(namesRows);

    VDbl_free(&coefRowsVec);
    VInt_free(&idxRowsVec);

    //}


    // printf("\n time CutP_compElemRepeatedCoef cut %f\n", ttimerepfim);
    //  printf("\n time add model %f\n", ttimeaddmodelfim);
    // fflush(stdout);


    //if(VERBOSE)  Res_setTCutsTotal(res,nround,cut, timeseparation);
    if(VERBOSE==2)
    {
        char namecut[256];
        if(cut==LPC_PREC) sprintf(namecut,"LPC_PREC");
        if(cut==LPC_CGCPU) sprintf(namecut,"LPC_CGCPU");
        if(cut==LPC_CG) sprintf(namecut,"LPC_CG");
        if(cut==LPC_CGGPU) sprintf(namecut,"LPC_CGGPU");
        if(cut==LPC_CGGPUR2) sprintf(namecut,"LPC_CGGPUR2");
        if(cut==LPC_RR) sprintf(namecut,"LPC_RR");
        if(cut==LPC_CLIQUE) sprintf(namecut,"LPC_CLIQUE");
        if(cut==LPC_ODDHOLES) sprintf(namecut,"LPC_ODDHOLES");
//        if(cut==LPC_JS) sprintf(namecut,"LPC_JS");

        printf("\n %s. %d cuts added. dominated %d. repeated %d. inactive %d. new %d. ", namecut, newcut+inactive, dominated, repeated, inactive, newcut);
        fflush(stdout);
    }



    // lp_write_lp(lp,"lp_original.lp");
    // getchar();
    int ncut = inactive+newcut;

    if(ncut)
    {
        if(VERBOSE==3)
        {
            //printf("\n%d new cuts  where added. nz %d.  separation time %f.  \n", ncut,  lp_nz(lp), timeseparation );
            printf("\n%d new cuts  where added. nz %d.  \n", ncut,  lp_nz(lp) );
            printf(" minViol %f, maxViol %f, avgViol %f.", Res_getMinViol(res,nround,cut), Res_getMaxViol(res,nround,cut), (double) Res_getSumViol(res,nround,cut)/  (double) Res_getNCutsTotal(res,nround,cut));
            printf("\n max violated cut: \n");
            for(int mv = 0 ; mv < Res_getNElementsMaxViol(res,nround,cut) ; mv++)
            {
                int elem = Res_getMaxViolCutA(res,nround,cut,mv);
                double co = Res_getMaxViolCutB(res,nround,cut,mv);
                printf(" %f*%d : %f + ", co,elem, xf[elem]);
            }
            printf(" \n min violated cut: \n");
            for(int mv = 0 ; mv < Res_getNElementsMinViol(res,nround,cut) ; mv++)
            {
                int elem = Res_getMinViolCutA(res,nround,cut,mv);
                double co = Res_getMinViolCutB(res,nround,cut,mv);
                printf(" %f*%d : %f + ", co,elem, xf[elem]);
            }
            printf("\n nMinElementsCuts %d, nMaxElementsCuts %d, avgElementsCuts %f.\n", Res_getNMinElementsCuts(res,nround,cut),  Res_getNMaxElementsCuts(res,nround,cut), (double) Res_getNElementsCuts(res,nround,cut)/ Res_getNCutsTotal(res,nround,cut));
        }
    }

    return ncut;
}

/* insert again the cuts of cutPool in lp in which were inactive -1 if the slack is inside the interval of the maximum allowed slack.
Cuts with -2 are dominated and they are not included again*/
void CutP_addCut( CutPool *cutP, LinearProgram* lp, int continuous, double maxslack, int nround)
{

    const double *xf = lp_x(lp);
    int nColsLP = lp_cols(lp);

    int idxR = 0, idxTC = 1, idxRhs = 2, idxSense = 3, idxSC = 4, idxCut = 5;
    double r = 0.0, tc = 0.0, ncut = 0.0, nrhs = 0.0, nsense = 0.0, sc = 0.0; //row, cut, size cut
    int init = 0, initcoef = 0, end = 0, endcoef = 0;
    //    char prefix[STR_SIZE];
    //  int nJobs = Inst_nJobs(cutP->inst);
    //    int ix[nJobs];
    //    char namecol[STR_SIZE];

    int cont=0;
    for(int n = 0 ; n < cutP->nHash ; n++)
    {

        int s = VDbl_size(cutP->hashCuts[n]); //size hash
        if(s == 0) continue;

        r = VDbl_get(cutP->hashCuts[n], idxR);
        tc = VDbl_get(cutP->hashCuts[n], idxTC);
        nrhs = VDbl_get(cutP->hashCuts[n], idxRhs);
        nsense = VDbl_get(cutP->hashCuts[n], idxSense);
        sc = VDbl_get(cutP->hashCuts[n], idxSC);
        ncut = VDbl_get(cutP->hashCuts[n], idxCut);
        init = idxCut+1;
        end = init+sc-1;
        initcoef = end+1;
        endcoef = end+sc;
        //  if(r==-1) printf("\nKey: %d Size: %d \n", n, s);
        //  if(r==-1) printf("1) idxR %d, r %f, idxTC %d, tc %f, idxRhs %d, nrhs %f, idxSense %d,  nsense %f, idxSC %d, sc %f, init %d, initcoef %d, end %d, endcoef %d, s %d \n", idxR, r, idxTC, tc, idxRhs, nrhs, idxSense, nsense, idxSC, sc, init, initcoef, end, endcoef, s);

        int endline = 0;
        // while( r != -1 && r!= -2) {
        while( r != -1 )
        {
            // printf("\nr == 1\n");
            idxCut = endcoef+6;
            idxSC = endcoef+5;
            idxSense = endcoef+4;
            idxRhs = endcoef+3;
            idxTC = endcoef+2;
            idxR = endcoef+1;
            if(idxR >= s)
            {
                idxR = 0;
                idxTC = 1;
                idxRhs = 2;
                idxSense = 3;
                idxSC = 4;
                idxCut = 5;
                init = idxCut+1;
                end = init+sc-1;
                initcoef = end+1;
                endcoef = end+sc;
                endline=1;
                break;
            }
            else
            {
                r = VDbl_get(cutP->hashCuts[n], idxR);
                tc = VDbl_get(cutP->hashCuts[n], idxTC);
                nrhs = VDbl_get(cutP->hashCuts[n], idxRhs);
                nsense = VDbl_get(cutP->hashCuts[n], idxSense);
                sc = VDbl_get(cutP->hashCuts[n], idxSC);
                ncut = VDbl_get(cutP->hashCuts[n], idxCut);
                init = idxCut+1;
                end = init+sc-1;
                initcoef = end+1;
                endcoef = end+sc;
                //    printf("3) idxR %d, r %f, idxTC %d, tc %f, idxRhs %d, nrhs %f, idxSense %d, nsense %f, idxSC %d, sc %f,  init %d, initcoef %d, end %d endcoef %d\n", idxR, r, idxTC, tc, idxRhs, nrhs, idxSense, nsense,idxSC, sc, init, initcoef, end, endcoef);
            }
        }

        if(endline) continue;

        int *idx;
        double *coef;
        ALLOCATE_VECTOR(idx, int, nColsLP);
        ALLOCATE_VECTOR(coef, double, nColsLP);

        int o = 0;
        double value = 0;

        int ihcoef = initcoef;
        for(int ihce = init ; ihce < init+sc && ihcoef < initcoef+sc && ihcoef < s; ihce++, ihcoef++)
        {
            //   printf("\n ihce %d == init %d +sc %f -1 && ihcoef %d == initcoef %d +sc %f -1\n", ihce, init, sc, ihcoef,initcoef,sc);

            int elem = VDbl_get(cutP->hashCuts[n],ihce);
            double coefs = VDbl_get(cutP->hashCuts[n],ihcoef);

            value += (xf[elem]*coefs);
            idx[o] = elem;
            coef[o] = coefs;
            //  if(r==-1)
            //      printf(" %f * %d ,", coefs, elem);
            //  fflush(stdout);
            o++;
            if(ihce == init+sc-1 && ihcoef == initcoef+sc-1)
            {
                //    if(r==-1) printf(" nElem %f, sense %f, rhs %f\n", sc, nsense, nrhs );
                //   printf(" -- \n\n");
                char name[STR_SIZE];
                end = ihce;
                endcoef = ihcoef;

                double slack = fabs(nrhs-value);
                if(slack < maxslack)
                {

                    char sens = (nsense== 0 ? 'L' : (nsense == 1 ? 'E' : 'G'));

/*                    if(tc==LPC_JS)
                        sprintf( name, "cutJS#%d_%d", lp_rows(lp)+1, nround);*/
                    if(tc==LPC_PREC)
                        sprintf( name, "cutPrec#%d_%d", lp_rows(lp)+1, nround);
                    if(tc==LPC_CGCPU)
                        sprintf( name, "cutCGCPU#%d_%d", lp_rows(lp)+1, nround);
                    if(tc== LPC_CLIQUE)
                        sprintf( name, "cutClique#%d_%d", lp_rows(lp)+1, nround );
                    if(tc== LPC_ODDHOLES)
                        sprintf( name, "cutODDHOLES#%d_%d", lp_rows(lp)+1, nround );
                    if(tc== LPC_RR)
                        sprintf( name, "cutRR#%d_%d", lp_rows(lp)+1, nround);
                    if(tc== LPC_CG)
                        sprintf( name, "cutCG#%d_%d", lp_rows(lp)+1, nround );
                    if(tc== LPC_CGGPU)
                        sprintf( name, "CGGPU#%d_%d", lp_rows(lp)+1, nround );
                    if(tc== LPC_CGGPUR2)
                        sprintf( name, "CGGPUR2#%d_%d", lp_rows(lp)+1, nround );
                    if(continuous)
                        lp_add_row(  lp, sc, idx, coef, name, sens, nrhs );
                    else
                        lp_add_cut( lp, sc, idx, coef, name, sens, nrhs );
                    VDbl_set(cutP->hashCuts[n],idxR,1);
                    VStr_set(cutP->namesCuts[n],ncut,name);
                    cont++;
                    //   printf("Add slack %f, maxslack %f\n", slack, maxslack);

                    //    printf("n %d, idxR %d, 1\n", n,idxR);
                }


                CLEAR_VECTOR(idx,int,nColsLP);
                CLEAR_VECTOR(coef,int,nColsLP);

                o = 0;
                value = 0;
                if(endcoef+1 == s)
                {
                    idxR = 0;
                    idxTC = 1;
                    idxRhs = 2;
                    idxSense = 3;
                    idxSC = 4;
                    idxCut = 5;
                    init = idxCut+1;
                    end = init+sc-1;
                    initcoef = end+1;
                    endcoef = end+sc;
                    endline=1;
                    ihcoef = initcoef;
                    ihce = init;
                    break;
                }


                idxR = endcoef+1;
                idxTC = endcoef+2;
                idxRhs = endcoef+3;
                idxSense = endcoef+4;
                idxSC = endcoef+5;
                idxCut = endcoef+6;
                r = VDbl_get(cutP->hashCuts[n], idxR);
                tc = VDbl_get(cutP->hashCuts[n], idxTC);
                nrhs = VDbl_get(cutP->hashCuts[n], idxRhs);
                nsense = VDbl_get(cutP->hashCuts[n], idxSense);
                sc = VDbl_get(cutP->hashCuts[n], idxSC);
                ncut = VDbl_get(cutP->hashCuts[n], idxCut);
                init = idxCut+1;
                end = init+sc-1;
                initcoef = end+1;
                endcoef = end+sc;
                ihce = init;
                ihcoef = initcoef;
                // printf("6.1) idxR %d, r %f, idxTC %d, tc %f, idxRhs %d, nrhs %f, idxSense %d, nsense %f, idxSC %d, sc %f, init %d, ihce %d, initcoef %d, ihcoef %d, end %d endcoef %d\n", idxR, r, idxTC, tc, idxRhs, nrhs, idxSense, nsense,idxSC, sc, init, ihce, initcoef, ihcoef, end, endcoef);

                // while( r != -1 && r!= -2) {
                while( r != -1 )
                {
                    //     printf("\nSec r == 1\n");
                    idxCut = endcoef+6;
                    idxSC = endcoef+5;
                    idxSense = endcoef+4;
                    idxRhs = endcoef+3;
                    idxTC = endcoef+2;
                    idxR = endcoef+1;
                    if(idxR >= s)
                    {
                        //      printf("\nendline\n");
                        idxR = 0;
                        idxTC = 1;
                        idxRhs = 2;
                        idxSense = 3;
                        idxSC = 4;
                        idxCut = 5;
                        r = 0.0;
                        tc = 0.0;
                        nrhs = 0.0;
                        nsense = 0.0;
                        sc = 0.0;
                        ncut = 0.0;
                        init = idxCut+1;
                        end = init+sc-1;
                        initcoef = end+1;
                        endcoef = end+sc;
                        ihcoef = initcoef;
                        ihce = init;
                        endline=1;
                        //    printf("7) idxR %d, r %f, idxTC %d, tc %f, idxRhs %d, nrhs %f, idxSense %d, nsense %f, idxSC %d, sc %f,  init %d, initcoef %d, end %d endcoef %d\n", idxR, r, idxTC, tc, idxRhs, nrhs, idxSense, nsense,idxSC, sc, init, initcoef, end, endcoef);
                        break;
                    }
                    else
                    {
                        r = VDbl_get(cutP->hashCuts[n], idxR);
                        tc = VDbl_get(cutP->hashCuts[n], idxTC);
                        nrhs = VDbl_get(cutP->hashCuts[n], idxRhs);
                        nsense = VDbl_get(cutP->hashCuts[n], idxSense);
                        sc = VDbl_get(cutP->hashCuts[n], idxSC);
                        ncut = VDbl_get(cutP->hashCuts[n], idxCut);
                        init = idxCut+1;
                        end = init+sc-1;
                        initcoef = end+1;
                        endcoef = end+sc;
                        ihcoef = initcoef;
                        ihce = init;
                        //    printf("8) idxR %d, r %f, idxTC %d, tc %f, idxRhs %d, nrhs %f, idxSense %d, nsense %f, idxSC %d, sc %f,  init %d, initcoef %d, end %d endcoef %d\n", idxR, r, idxTC, tc, idxRhs, nrhs, idxSense, nsense,idxSC, sc, init, initcoef, end, endcoef);
                    }
                }

                // getchar();
                ihcoef = initcoef-1;
                ihce = init-1;
                if(endline) break;
            }
        }
        free(coef);
        free(idx);
    }
    if(VERBOSE==2)  printf("\nnumbers of rows added %d, lp_rows after add cuts %d. \n", cont, lp_rows(lp));
}

/*remove cut from lp if it has extrapolated the established maximum slack allowed to the next rounds
  the cut remain at the CutPool with stats inactive (-1)
*/
void CutP_removeCut( CutPool *cutP, LinearProgram* lp,  double maxslack)
{
    //    double removecuttime = omp_get_wtime();
    double findSlack = 0;
    double findElem = 0;
    double findRow = 0;

    //    const double *xf = lp_x(lp);
    int nColsLP = lp_cols(lp);
    int nRows = lp_rows(lp);
    int *rowstoremove;
    ALLOCATE_VECTOR(rowstoremove, int, nRows);
    int cont = 0 ;

    double tt = omp_get_wtime();
    double *slacks = lp_row_slack(lp);
    findSlack +=  (omp_get_wtime()-tt);
    //  printf("\nnRows %d: ", nRows);
    for(int n = 0 ; n < nRows ; n++)
    {

        int aux = 0;

        //  printf("\nSlack %f < maxslack %f ", slacks[n], maxslack);
        if( slacks[n] < maxslack) continue;

        char name[STR_SIZE];
        lp_row_name(lp, n, name);
        if(strncmp(name,"cutR",4)== 0)
            aux = LPC_RR;
        else if(strncmp(name,"cutCl",5)== 0)
            aux = LPC_CLIQUE;
        else if(strncmp(name,"cutODDHOLES",11)== 0)
            aux = LPC_ODDHOLES;
        else if(strncmp(name,"cutP",4)== 0)
            aux = LPC_PREC;
        else if(strncmp(name,"cutCG",5)== 0)
            aux = LPC_CG;
        else if(strncmp(name,"CGGPUR2",7)== 0)
            aux = LPC_CGGPUR2;
        else if(strncmp(name,"CGGPU",5)== 0)
            aux = LPC_CGGPU;
        else if(strncmp(name,"cutCGCPU",4)== 0)
            aux = LPC_CGCPU;
/*        else if(strncmp(name,"cutJS",5)== 0)
            aux = LPC_JS;*/
        else
            continue;

        // printf("\n name res to be removed %s ", name);
        //  fflush(stdout);

        //printf("Slacks %f, maxslack %f\n", slacks[n], maxslack);

        int *idx;
        double *coef;
        ALLOCATE_VECTOR(idx, int, nColsLP);
        ALLOCATE_VECTOR(coef, double, nColsLP);


        double tttt = omp_get_wtime();
        int c = lp_row(lp,n,idx,coef);

        //   printf("Cut to be removed: \n");
        //  for(int nn = 0 ; nn < c ; nn++){
        //       printf(" %f * %d ", coef[nn], idx[nn]);
        // }
        // printf("\n");

        char sen = lp_sense(lp,n);
        int sense = (sen == 'L' ? 0 : (sen == 'E' ? 1 : 2));
        double rhs = lp_rhs(lp,n);
        int key = hash_code_vint(  c, idx,  nColsLP );
        /*
                 printf("Remover cut name %s: \n", name);
                  for(int sel = 0 ; sel<c ; sel++) {
                      printf(" %f * %d ", coef[sel], idx[sel]);
                 }
                 printf("\n");
                  printf(" key removing %d \n", key );
                           getchar();
        */

        findRow +=  (omp_get_wtime()-tttt);

        double ttt = omp_get_wtime();
        rowstoremove[cont] = n;
        IntDblTuple posrow = cutP_findElement(cutP,key,c,idx,coef,rhs,sense,nColsLP,aux);
        /*  if(posrow.b==0){
            printf(" removing on key %d idRow %d value %f ncut %f ", key, posrow.a, posrow.b, posrow.c);
           // printf("Remover cut name %s: \n", name);
            for(int sel = 0 ; sel<c ; sel++) {
                printf(" %f * %d ", coef[sel], idx[sel]);
            }
            printf("\n");
            printf(" rhs %f , key removing %d \n", rhs, key );
            getchar();
             CutP_printHash(cutP); getchar();
              fflush(stdout);
          }*/
        //if(posrow.b==0) getchar();
        //CutP_printHash(cutP,lp); getchar();
        // if(posrow.b==1)
        CutP_setIdRow( cutP, key, posrow.a, -1);
        //  posrow = cutP_findElement(cutP,key,c,idx,coef,rhs,sense,nColsLP,aux);
        //  int idCut = CutP_getCut(cutP,key,posrow.a);
        // if(posrow.b==0){
        //    printf(" new removing on key %d idRow %d name %s, value %f ", key, posrow.a, VStr_get(cutP->namesCuts[key], posrow.c), VDbl_get(cutP->hashCuts[key], posrow.a));
        //  fflush(stdout);
        // }


        cont++;
        //   printf("time find elements %f. \n",  (omp_get_wtime()-ttt));
        findElem += (omp_get_wtime()-ttt);

        free(coef);
        free(idx);
    }
    //   CutP_printHash(cutP,lp);
    //  printf("time find row %f, time find slacks %f time find elem %f total time find rows to remove %f. \n", findRow,  findSlack, findElem, (omp_get_wtime()-removecuttime));
    //    double rrr = omp_get_wtime();
    lp_remove_rows(lp,cont,rowstoremove);

    free(rowstoremove);
    if(VERBOSE==2) printf("\nnumber of rows removed %d, lp_rows after remove cuts %d. ", cont, lp_rows(lp));
    //   getchar();

}

/*starts the procedure to separate in parallel cuts to the cutting plane*/
int CutP_separationParallel(CutPool *cutP, FILE *fp2, const Parameters *par,  Results *res, int cr, int cp, int cc, int co, int cg, int cggpur2, int cgcpu,  int** maxTJM, int* maxTJ, double startT, double timeLeft)
{

    double initT = omp_get_wtime();
    int newCuts = 0;
    int nround = Res_getRound(res);
    double _time = 0.0;
    //    double *xf = lp_x(cutP->mip);

    CutC *ccr= CutC_create(cutP->inst);
    CutPR *cpr = CutPR_create(cutP->inst);
//    CutJS *ccjs = CutJS_create(cutP->inst);
    CutCL *ccl = CutCL_create(cutP->inst);
    CutOH *coh = CutOH_create(cutP->inst);
    CutCG *ccgcpu;
    CutCG *ccggpur2;
    CutCG *ccg;
    int thr = cr+cp+cc+co+cg+cgcpu+cggpur2;//+cjs;
    double cutTimeCR = 0.0,  cutTimeCP = 0.0,   cutTimeCC = 0.0, cutTimeCO = 0.0, cutTimeCG = 0.0; // cutTimeJS = 0.0,

    int origCols[lp_cols( cutP->mip )];
    for ( int i=0 ; (i<lp_cols( cutP->mip )) ; ++i )
        origCols[i] = i;

    _time = ( (double) timeLeft - (omp_get_wtime()-initT) );
    CGraph *cgraph;

    if(cc || co     )
        cgraph = CutP_compute_conflicts_create_graph(cutP->mip, origCols, cutP->mip, par->inst,_time);

    _time = ( (double) timeLeft - (omp_get_wtime()-initT) );

    #pragma omp parallel sections shared(_time,cutP,res,par,cgraph,origCols) num_threads(thr)
    {
        #pragma omp section
        {

            if(cr)
            {
               // printf("cr ID %d\n", omp_get_thread_num()); fflush(stdout);//getchar();
                double initTCR = omp_get_wtime();
                //CutC_add_cuts_cover_parallel( ccr, cutP->mip, origCols, cutP->mip, par->inst,_time, par->continuous, par->lifting, Res_getRound(res));
                if(par->lifting)
                    CutC_add_cuts_cover_model_parallel( ccr, cutP->mip, origCols, cutP->mip, par->inst,_time, par->continuous, par->lifting, Res_getRound(res));
                else
                    CutC_add_cuts_cover_parallel( ccr, cutP->mip, origCols, cutP->mip, par->inst,_time, par->continuous, par->lifting, Res_getRound(res));
                //printf("end cr ID %d\n", omp_get_thread_num());fflush(stdout);//getchar();
                cutTimeCR = (omp_get_wtime()-initTCR);
            }
        }
        #pragma omp section
        {

            if(cp)
            {
                //printf("cp ID %d\n", omp_get_thread_num());fflush(stdout);//getchar();
                double initTCP = omp_get_wtime();
                //CutPR_add_cuts_precedence_parallel( cpr, cutP->mip, origCols, cutP->mip, par->inst,_time, par->continuous, par->lifting, maxTJM, maxTJ, par->maxcuts,  Res_getRound(res) );
                CutPR_add_cuts_precedence_parallel( cpr, cutP->mip,  origCols, cutP->mip, par->inst,_time, par->continuous, par->lifting, par->maxcuts,  Res_getRound(res) );
                cutTimeCP = (omp_get_wtime()-initTCP);
               // printf("end cp ID %d\n", omp_get_thread_num());
              //  fflush(stdout);//getchar();
            }
        }
       /*  #pragma omp section
        {

            if(cjs)
            {
                //printf("cp ID %d\n", omp_get_thread_num());fflush(stdout);//getchar();
                double initTJS = omp_get_wtime();
                CutJS_add_cuts_jobset_parallel( ccjs, cutP->mip,  origCols, cutP->mip, par->inst,_time, par->continuous, par->lifting, par->maxcuts,  Res_getRound(res) );
                cutTimeJS = (omp_get_wtime()-initTJS);
               // printf("end cp ID %d\n", omp_get_thread_num());
              //  fflush(stdout);//getchar();
            }
        }*/
        #pragma omp section
        {
            if(cc)
            {
               // printf("cc ID %d\n", omp_get_thread_num());fflush(stdout);//getchar();
                double initTCC = omp_get_wtime();
                CutCL_add_cuts_conflicts_clique_parallel( cgraph, ccl, cutP->mip, par->inst,_time, par->maxcuts,  Res_getRound(res));
                cutTimeCC = (omp_get_wtime()-initTCC);
             //  printf("end cc ID %d\n", omp_get_thread_num());
//                fflush(stdout);//getchar();
            }
        }
        #pragma omp section
        {

            if(co)
            {
               // printf("co ID %d\n", omp_get_thread_num());fflush(stdout);//getchar();
                double initTCO = omp_get_wtime();
                CutOH_add_cuts_conflicts_odd_parallel( cgraph, coh, cutP->mip, par->inst,_time, par->maxcuts,   Res_getRound(res));
                cutTimeCO = (omp_get_wtime()-initTCO);
               // printf("end co ID %d\n", omp_get_thread_num());fflush(stdout);//getchar();
            }
        }
        #pragma omp section
        {
            if(cg)
            {
                int ID = omp_get_thread_num();
                //  printf("MipC_cutCG: start on thread %d\n", ID); fflush(stdout);
                double initTCG = omp_get_wtime();
                double timeconst = 0;

                //Selected rows to GPU

                //ccg = CutCG_creat_and_identify_rows( cutP->mip, origCols, cutP->mip, par->inst, _time, par->continuous, maxTJM, par->maxcuts,par->mininstant, par->maxinstant, par->jump,  cutP);
                int horizon=0;
                //ccg = CutCG_create_and_identify_rows_all( cutP->mip, origCols, cutP->mip, par->inst, _time,par->mininstant, par->maxinstant, par->jump,  Res_getRound(res));
                ccg = CutCG_create_and_identify_rows_all_conflitcts( cutP->mip, origCols, cutP->mip, par->inst, _time,par->mininstant, par->maxinstant, par->jump,  Res_getRound(res), &horizon);

                //CutCG_print(ccg);
                timeconst = ( (double) _time - (omp_get_wtime()-initTCG));

                CutCG_add_cuts_model_CG_parallel( ccg, cgraph, LPC_CG, cutP->mip, par->inst, timeconst, par->maxcuts,  Res_getRound(res), horizon);

                cutTimeCG = (omp_get_wtime()-initTCG);

            }
        }
        #pragma omp section
        {
            if(cgcpu)
            {
              //  printf("co CGCPU %d\n", omp_get_thread_num());fflush(stdout);//getchar();
                double initTCG = omp_get_wtime();
                parameters_ccg *parCCG =(parameters_ccg*)malloc(sizeof(parameters_ccg));
                double timeconst = 0;

                //Selected rows to CG
                int horizon=0;
                ccgcpu = CutCG_create_and_identify_rows_all_conflitcts( cutP->mip, origCols, cutP->mip, par->inst, _time, par->mininstant, par->maxinstant, par->jump,  Res_getRound(res), &horizon);
                //CutCG_create_and_identify_rows_all
                //CutCG_create_and_identify_rows_all_conflitcts
                setParameters_ccg(parCCG, 0); //0 - normal mode, 1- intermediate mode, 2- hard mode

                int cont=0;

                Cut_gpu *ccg;
                int nC1,nV1;

                nV1 = VDbl_size(ccgcpu->xfElemPP);
                nC1 = VDbl_size(ccgcpu->rowrhs);
                for( int r = 0 ; r < VDbl_size(ccgcpu->rowrhs) ; r++)
                    cont+= VDbl_size(ccgcpu->rowCoef[r]);
                ccg = fillStruct(ccgcpu,parCCG->precision, nV1, nC1,cont);

                timeconst = ( (double) _time - (omp_get_wtime()-initTCG));
            //    printf("\n N ccg %d ", ccg->numberConstrains); //getchar();



                ccg = initial_runCPU(ccg,parCCG->maxDenominator,parCCG->precision,timeconst);
                timeconst = ( (double) _time - (omp_get_wtime()-initTCG));
                if(timeconst>=1){
                    ccg  = second_phase_runCPU(ccg,parCCG->numberMaxConst,parCCG->nRuns,parCCG->maxDenominator,parCCG->precision,parCCG->nSizeR1, timeconst);
                }
                returnCutCG(ccg,ccgcpu,parCCG->precision);
                free(parCCG);
                free(ccg);

                //end CG
                timeconst = ( (double) _time - (omp_get_wtime()-initTCG) );
                CutCG_add_cuts_CG_parallel( ccgcpu, cgraph, LPC_CGCPU, cutP->mip, par->inst, timeconst, par->continuous,  par->maxcuts,  Res_getRound(res), horizon);

                cutTimeCG = (omp_get_wtime()-initTCG);
              //  printf("end cg ID %d\n", omp_get_thread_num()); fflush(stdout);//getchar();
            }
        }

         #pragma omp section
        {
            if(cggpur2)
            {

#ifdef __NVCC__
                //setGpuThread(0);
#endif // __NVCC__
               // int ID = omp_get_thread_num();
               // printf("cg %d MipC_cutCGGPU Rank 2: start on thread %d\n", cggpur2, ID);
               // fflush(stdout);
                double initTCG = omp_get_wtime();
                parameters_ccg *parCCG =(parameters_ccg*)malloc(sizeof(parameters_ccg));

                double timeconst = 0;
                //Selected rows to GPU
                int horizon=0;
                ccggpur2 = CutCG_create_and_identify_rows_all_conflitcts( cutP->mip, origCols, cutP->mip, par->inst, _time,par->mininstant, par->maxinstant, par->jump,  Res_getRound(res), &horizon);
                timeconst = ( (double) _time - (omp_get_wtime()-initTCG) );
                setParameters_ccg(parCCG, 0); //0 - normal mode, 1- intermediate mode, 2- hard mode


#ifdef __NVCC__
                //start CG
                //printf("AQUI!");
                //getchar();
                int cont=0;
                int nRR_cggpu = Inst_nResR(par->inst);
                int nC1,nV1;

                for( int r = 0 ; r < VDbl_size(ccggpur2->rowrhs) ; r++)
                    cont+= VDbl_size(ccggpur2->rowCoef[r]);

                nV1 = VDbl_size(ccggpur2->xfElemPP);
                nC1 = VDbl_size(ccggpur2->rowrhs);

                Cut_gpu *ccg;
                Cut_gpu_aux *ccg_aux;
                ccg = fillStruct(ccggpur2,parCCG->precision,nV1,nC1,cont);

                int n_gpu = verifyGpu();
                //printf("verify gpu %d\n", n_gpu);
                fflush(stdout);
                int nBlocks, nThreads, nRuns_temp;
                int i;
                int pos_R1 = 0;
                double aux;
                int n_cuts,nRepeat = 1;

                ccg_aux = AllocationStructCutAux(nC1,nV1,cont);
                returnDimension(&nBlocks,&nThreads, parCCG->nRuns,ccg->numberConstrains);
                n_cuts = ccg->numberConstrains;

                ccg = zeroHalf_runGPU(ccg, ccg_aux,20,parCCG->precision,nThreads, nBlocks);
//                getchar();
                //ccg = initial_runGPU(ccg, ccg_aux,parCCG->maxDenominator,parCCG->precision,1,nThreads,nBlocks, nRR_cggpu);
 //               if(n_cuts!=ccg->numberConstrains)
//                    ccg_aux = reallocCut(ccg,ccg_aux);
                n_cuts = ccg->numberConstrains;
                if(parCCG->nRuns < 0.7*(nBlocks*nThreads))
                {
                    nThreads = (parCCG->nRuns/nBlocks) ;
                }
                else
                {
                    aux = (float)(parCCG->nRuns - 0.7*(nThreads*nBlocks))/(float)(0.7*(nThreads*nBlocks));
                    nRepeat += ceil(aux);
                    nThreads = (0.7*(nThreads*nBlocks))/nBlocks;
                }
                for(i=1; i<=nRepeat; i++)
                {
                    timeconst = ( (double) _time - (omp_get_wtime()-initTCG));
                    fflush(stdout);
                    if((nRepeat>1)&&(i == nRepeat))
                    {
                        nThreads = (parCCG->nRuns - (nThreads*nBlocks)*(i-1))/nBlocks;
                    }
                    nRuns_temp = nThreads*nBlocks;
                    if(timeconst>= 1){
                        /*ATUALIZAR TEMPO DENTRO DO METODO timeGPU*/
                        ccg = second_phase_runGPU(ccg, ccg_aux, parCCG->numberMaxConst,nRuns_temp,parCCG->maxDenominator,parCCG->precision, nBlocks, nThreads, &pos_R1, parCCG->nSizeR1,timeconst);
                    }
                }
                returnCutCG(ccg,ccggpur2,parCCG->precision);

                free(parCCG);
                free(ccg);
                free(ccg_aux);
                //end CG
#endif // __NVCC__

                double timegpu = ( (double) timeLeft - (omp_get_wtime()-initT) );
                CutCG_add_cuts_CG_parallel( ccggpur2, cgraph, LPC_CGGPUR2, cutP->mip, par->inst, timegpu, par->continuous,  par->maxcuts,  Res_getRound(res), horizon);

                cutTimeCG = (omp_get_wtime()-initTCG);
             //printf("end cg ID %d\n", omp_get_thread_num()); fflush(stdout);//getchar();
            }
        }
    }



    int cutcr = 0;
    int cutprec = 0;
//    int cutjs = 0;
    int cutcgcpu =0;
    int cutclique = 0;
    int cutoddholes = 0;
    int cutcggpur2 =0;
    int cutcg=0;

    if(cr)
    {
        _time = ( (double) timeLeft - (omp_get_wtime()-initT) );
        double initTCR = omp_get_wtime();
        for(int r = 0; r < Inst_nResR(par->inst); r++)
        {
            //   cutcr += CutP_addSeparatedCuts(cutP,res,cutP->mip,VInt_getPtr(ccr->cutElem[r]),VDbl_getPtr(ccr->cutCoef[r]),VDbl_getPtr(ccr->cutrhs[r]),VDbl_getPtr(ccr->cutviolation[r]),VStr_ptr(ccr->cutname[r]),VInt_getPtr(ccr->cutsense[r]),VInt_getPtr(ccr->cutdominated[r]),par->continuous,nround,LPC_RR,_time);//,ccr->timeseparation,ccr->nAddCuts);
            cutcr += CutP_addSeparatedCuts(cutP,res,cutP->mip,ccr->cutElem[r],ccr->cutCoef[r],ccr->cutrhs[r],ccr->cutviolation[r],ccr->cutname[r],ccr->cutsense[r],ccr->cutdominated[r],par->continuous,nround,LPC_RR,_time);//,ccr->timeseparation,ccr->nAddCuts);
        }
        newCuts += cutcr;
        cutTimeCR += (omp_get_wtime()-initTCR);
        if(VERBOSE) Res_setTCutsTotal(res,Res_getRound(res),LPC_RR,cutTimeCR);
       // printf("\nCR added: %d\n", cutcr); fflush(stdout);

    }

    if(cp)
    {
        _time = ( (double) timeLeft - (omp_get_wtime()-initT) );
        double initTCP = omp_get_wtime();
        cutprec = CutP_addSeparatedCuts(cutP,res,cutP->mip,cpr->cutElem,cpr->cutCoef,cpr->cutrhs,cpr->cutviolation,cpr->cutname,cpr->cutsense,cpr->cutdominated, par->continuous,nround,LPC_PREC,_time);//,cpr->timeseparation,cpr->nAddCuts);
        newCuts += cutprec;
        cutTimeCP += (omp_get_wtime()-initTCP);
        if(VERBOSE) Res_setTCutsTotal(res,Res_getRound(res),LPC_PREC,cutTimeCP);
      //  printf("\nCP added: %d\n", cutprec); fflush(stdout);
    }

  /*  if(cjs)
    {
        _time = ( (double) timeLeft - (omp_get_wtime()-initT) );
        double initTJS = omp_get_wtime();
        cutjs = CutP_addSeparatedCuts(cutP,res,cutP->mip,ccjs->cutElem,ccjs->cutCoef,ccjs->cutrhs,ccjs->cutviolation,ccjs->cutname,ccjs->cutsense,ccjs->cutdominated, par->continuous,nround,LPC_JS,_time);//,cpr->timeseparation,cpr->nAddCuts);
        newCuts += cutjs;
        cutTimeJS += (omp_get_wtime()-initTJS);
        if(VERBOSE) Res_setTCutsTotal(res,Res_getRound(res),LPC_JS,cutTimeJS);
      //  printf("\nCP added: %d\n", cutprec); fflush(stdout);
    }*/

    if(cc)
    {
        _time = ( (double) timeLeft - (omp_get_wtime()-initT) );
        double initTCC = omp_get_wtime();
        cutclique = CutP_addSeparatedCuts(cutP,res,cutP->mip,ccl->cutElem,ccl->cutCoef,ccl->cutrhs,ccl->cutviolation,ccl->cutname,ccl->cutsense,ccl->cutdominated,par->continuous,nround,LPC_CLIQUE,_time);//,ccl->timeseparation,ccl->nAddCuts);ĵ
        newCuts += cutclique;
        cutTimeCC += (omp_get_wtime()-initTCC);
        if(VERBOSE) Res_setTCutsTotal(res,Res_getRound(res),LPC_CLIQUE,cutTimeCC);
        //printf("\nCC added: %d\n", cutclique);  fflush(stdout);
    }

    if(co)
    {
        _time = ( (double) timeLeft - (omp_get_wtime()-initT) );
        double initTCO = omp_get_wtime();
        cutoddholes = CutP_addSeparatedCuts(cutP,res,cutP->mip,coh->cutElem,coh->cutCoef,coh->cutrhs,coh->cutviolation,coh->cutname,coh->cutsense,coh->cutdominated,par->continuous,nround,LPC_ODDHOLES,_time);//,coh->timeseparation,coh->nAddCuts);
        newCuts += cutoddholes;
        cutTimeCO += (omp_get_wtime()-initTCO);
        if(VERBOSE) Res_setTCutsTotal(res,Res_getRound(res),LPC_ODDHOLES,cutTimeCO);
        //printf("\nCO added: %d\n", cutoddholes); fflush(stdout);
    }

    if(cgcpu)
    {

        _time = ( (double) timeLeft - (omp_get_wtime()-initT) );
        double initTCG = omp_get_wtime();
        if(VDbl_size(ccgcpu->cutrhs)>0)
        {
            cutcgcpu = CutP_addSeparatedCuts(cutP,res,cutP->mip,ccgcpu->cutElem,ccgcpu->cutCoef,ccgcpu->cutrhs,ccgcpu->cutviolation,ccgcpu->cutname,ccgcpu->cutsense,ccgcpu->cutdominated,par->continuous,nround,LPC_CGCPU,_time);//,ccgcpu->timeseparation,ccgcpu->nAddCuts);
            newCuts += cutcgcpu;
        }
        CutCG_free(&ccgcpu);

        cutTimeCG += (omp_get_wtime()-initTCG);
        if(VERBOSE) Res_setTCutsTotal(res,Res_getRound(res),LPC_CGCPU,cutTimeCG);
       // printf(" CG added: %d ", cutcgcpu);fflush(stdout);
    }
    if(cggpur2)
    {

        _time = ( (double) timeLeft - (omp_get_wtime()-initT) );
        double initTCG = omp_get_wtime();
        if(VDbl_size(ccggpur2->cutrhs)>0)
        {
            cutcggpur2 = CutP_addSeparatedCuts(cutP,res,cutP->mip,ccggpur2->cutElem,ccggpur2->cutCoef,ccggpur2->cutrhs,ccggpur2->cutviolation,ccggpur2->cutname,ccggpur2->cutsense,ccggpur2->cutdominated,par->continuous,nround,LPC_CGGPUR2,_time);//,ccggpur2->timeseparation,ccggpur2->nAddCuts);
            newCuts += cutcggpur2;
        }
        CutCG_free(&ccggpur2);

        cutTimeCG += (omp_get_wtime()-initTCG);
        if(VERBOSE) Res_setTCutsTotal(res,Res_getRound(res),LPC_CGGPUR2,cutTimeCG);
    }

    if(cg)
    {
        _time = ( (double) timeLeft - (omp_get_wtime()-initT) );
        double initTCG = omp_get_wtime();
        //  printf("add CG\n");fflush(stdout);
        cutcg = CutP_addSeparatedCuts(cutP,res,cutP->mip,ccg->cutElem,ccg->cutCoef,ccg->cutrhs,ccg->cutviolation,ccg->cutname,ccg->cutsense,ccg->cutdominated,par->continuous,nround,LPC_CG,_time);//,ccg->timeseparation,ccg->nAddCuts);
        Res_setUsoUR(res,nround,ccg->contUR);
        Res_setUsoURN(res,nround,ccg->contURN);
        Res_setUsoUM(res,nround,ccg->contUM);
        Res_setUsoUCI(res,nround,ccg->contUCI);
        Res_setUsoUCP(res,nround,ccg->contUCP);
        Res_setTotalUsoUR(res,ccg->contUR);
        Res_setTotalUsoURN(res,ccg->contURN);
        Res_setTotalUsoUM(res,ccg->contUM);
        Res_setTotalUsoUCI(res,ccg->contUCI);
        Res_setTotalUsoUCP(res,ccg->contUCP);
        Res_setNJumps(res,ccg->nJumps);
        newCuts += cutcg;
        cutTimeCG += (omp_get_wtime()-initTCG);
        if(VERBOSE) Res_setTCutsTotal(res,Res_getRound(res),LPC_CG,cutTimeCG);
        CutCG_free(&ccg);
    }

    if(cg)
        fprintf( fp2, "%d %d %d %d %d ", cutcr, cutprec, cutclique, cutoddholes,cutcg);
    if(cggpur2)
        fprintf( fp2, "%d %d %d %d %d ", cutcr, cutprec, cutclique, cutoddholes,cutcggpur2);
    if(cgcpu)
        fprintf( fp2, "%d %d %d %d %d ", cutcr, cutprec, cutclique, cutoddholes,cutcgcpu);

    cutP->totalTimeSeparation += (double) (omp_get_wtime()-initT);

    CutCL_free(&ccl);
    CutC_free(&ccr);
    CutPR_free(&cpr);
//    CutJS_free(&ccjs);
    CutOH_free(&coh);
    if(cc || co)
        cgraph_free( &cgraph );

    return newCuts;

}


/*starts the procedure to separate in parallel cuts to the cutting plane it already receives the cgraph created by all variables. Separated strategies that uses cgraph creates a induced subgraph*/
int CutP_separationParallelCgraph(CutPool *cutP, FILE *fp2, CGraph *cgraph, const Parameters *par,  Results *res, int cr, int cp, int cc, int co, int cg, int cggpur2, int cgcpu,  int** maxTJM, int* maxTJ, double startT, double timeLeft)
{


    double initT = omp_get_wtime();
    int newCuts = 0;
    int nround = Res_getRound(res);
    double _time = 0.0;

    CutC *ccr= CutC_create(cutP->inst);
    CutPR *cpr = CutPR_create(cutP->inst);
//    CutJS *ccjs = CutJS_create(cutP->inst);
    CutCL *ccl = CutCL_create(cutP->inst);
    CutOH *coh = CutOH_create(cutP->inst);
    CutCG *ccgcpu;
    CutCG *ccggpur2;
    CutCG *ccg;
    int thr = cr+cp+cc+co+cg+cgcpu+cggpur2;//+cjs;
    double cutTimeCR = 0.0,  cutTimeCP = 0.0,  cutTimeCC = 0.0, cutTimeCO = 0.0, cutTimeCG = 0.0; // cutTimeJS = 0.0,

    int origCols[lp_cols( cutP->mip )];
    for ( int i=0 ; (i<lp_cols( cutP->mip )) ; ++i )
        origCols[i] = i;

    _time = ( (double) timeLeft - (omp_get_wtime()-initT) );

    #pragma omp parallel sections shared(_time,cutP,res,par,cgraph,origCols) num_threads(thr)
    {
        #pragma omp section
        {

            if(cr)
            {
               // printf("cr ID %d\n", omp_get_thread_num()); fflush(stdout);//getchar();
                double initTCR = omp_get_wtime();
                //CutC_add_cuts_cover_parallel( ccr, cutP->mip, origCols, cutP->mip, par->inst,_time, par->continuous, par->lifting, Res_getRound(res));
                if(par->lifting)
                    CutC_add_cuts_cover_model_parallel( ccr, cutP->mip, origCols, cutP->mip, par->inst,_time, par->continuous, par->lifting, Res_getRound(res));
                else
                    CutC_add_cuts_cover_parallel( ccr, cutP->mip, origCols, cutP->mip, par->inst,_time, par->continuous, par->lifting, Res_getRound(res));
                //printf("end cr ID %d\n", omp_get_thread_num());fflush(stdout);//getchar();
                cutTimeCR = (omp_get_wtime()-initTCR);
            }
        }
        #pragma omp section
        {

            if(cp)
            {
                //printf("cp ID %d\n", omp_get_thread_num());fflush(stdout);//getchar();
                double initTCP = omp_get_wtime();
                //CutPR_add_cuts_precedence_parallel( cpr, cutP->mip, origCols, cutP->mip, par->inst,_time, par->continuous, par->lifting, maxTJM, maxTJ, par->maxcuts,  Res_getRound(res) );
                CutPR_add_cuts_precedence_parallel( cpr, cutP->mip,  origCols, cutP->mip, par->inst,_time, par->continuous, par->lifting, par->maxcuts,  Res_getRound(res) );
                cutTimeCP = (omp_get_wtime()-initTCP);
               // printf("end cp ID %d\n", omp_get_thread_num());
              //  fflush(stdout);//getchar();
            }
        }
      /*   #pragma omp section
        {

            if(cjs)
            {
                //printf("cp ID %d\n", omp_get_thread_num());fflush(stdout);//getchar();
                double initTJS = omp_get_wtime();
                CutJS_add_cuts_jobset_parallel( ccjs, cutP->mip,  origCols, cutP->mip, par->inst,_time, par->continuous, par->lifting, par->maxcuts,  Res_getRound(res) );
                cutTimeJS = (omp_get_wtime()-initTJS);
               // printf("end cp ID %d\n", omp_get_thread_num());
              //  fflush(stdout);//getchar();
            }
        }*/
        #pragma omp section
        {
            if(cc)
            {
               // printf("cc ID %d\n", omp_get_thread_num());fflush(stdout);//getchar();
                double initTCC = omp_get_wtime();
                CutCL_add_cuts_conflicts_clique_parallel( cgraph, ccl, cutP->mip,par->inst,_time, par->maxcuts,  Res_getRound(res));
                cutTimeCC = (omp_get_wtime()-initTCC);
             //  printf("end cc ID %d\n", omp_get_thread_num());
//                fflush(stdout);//getchar();
            }
        }
        #pragma omp section
        {

            if(co)
            {
               // printf("co ID %d\n", omp_get_thread_num());fflush(stdout);//getchar();
                double initTCO = omp_get_wtime();
                CutOH_add_cuts_conflicts_odd_parallel( cgraph, coh, cutP->mip, par->inst,_time, par->maxcuts,   Res_getRound(res));
                cutTimeCO = (omp_get_wtime()-initTCO);
               // printf("end co ID %d\n", omp_get_thread_num());fflush(stdout);//getchar();
            }
        }
        #pragma omp section
        {
            if(cg)
            {
                int ID = omp_get_thread_num();
                //  printf("MipC_cutCG: start on thread %d\n", ID);
                fflush(stdout);
                double initTCG = omp_get_wtime();
                double timeconst = 0;

                //Selected rows to GPU

                //ccg = CutCG_creat_and_identify_rows( cutP->mip, origCols, cutP->mip, par->inst, _time, par->continuous, maxTJM, par->maxcuts,par->mininstant, par->maxinstant, par->jump,  cutP);
                int horizon=0;
                ccg = CutCG_create_and_identify_rows_all_conflitcts_cgraph( cutP->mip, cgraph, origCols, cutP->mip, par->inst, _time,par->mininstant, par->maxinstant, par->jump,  Res_getRound(res), &horizon);

                //CutCG_print(ccg);
                timeconst = ( (double) _time - (omp_get_wtime()-initTCG));

                CutCG_add_cuts_model_CG_parallel( ccg, cgraph, LPC_CG, cutP->mip, par->inst, timeconst,  par->maxcuts,  Res_getRound(res), horizon);

                cutTimeCG = (omp_get_wtime()-initTCG);

            }
        }
        #pragma omp section
        {
            if(cgcpu)
            {
//                printf("co CGCPU %d\n", omp_get_thread_num());fflush(stdout);//getchar();
                double initTCG = omp_get_wtime();
                parameters_ccg *parCCG =(parameters_ccg*)malloc(sizeof(parameters_ccg));
                double timeconst = 0;

                //Selected rows to CG
                int horizon=0;
                ccgcpu = CutCG_create_and_identify_rows_all_conflitcts_cgraph( cutP->mip, cgraph, origCols, cutP->mip, par->inst, _time, par->mininstant, par->maxinstant, par->jump,  Res_getRound(res), &horizon);
                //CutCG_create_and_identify_rows_all
                //CutCG_create_and_identify_rows_all_conflitcts
                setParameters_ccg(parCCG, 0); //0 - normal mode, 1- intermediate mode, 2- hard mode

                int cont=0;

                Cut_gpu *ccg;
                int nC1,nV1;

                nV1 = VDbl_size(ccgcpu->xfElemPP);
                nC1 = VDbl_size(ccgcpu->rowrhs);
                for( int r = 0 ; r < VDbl_size(ccgcpu->rowrhs) ; r++)
                    cont+= VDbl_size(ccgcpu->rowCoef[r]);
                ccg = fillStruct(ccgcpu,parCCG->precision, nV1, nC1,cont);

                timeconst = ( (double) _time - (omp_get_wtime()-initTCG));
            //    printf("\n N ccg %d ", ccg->numberConstrains); //getchar();
                ccg = initial_runCPU(ccg,parCCG->maxDenominator,parCCG->precision,timeconst);
                timeconst = ( (double) _time - (omp_get_wtime()-initTCG));
                if(timeconst>=1){
                    ccg  = second_phase_runCPU(ccg,parCCG->numberMaxConst,parCCG->nRuns,parCCG->maxDenominator,parCCG->precision,parCCG->nSizeR1, timeconst);
                }
                returnCutCG(ccg,ccgcpu,parCCG->precision);
                free(parCCG);
                free(ccg);

                //end CG
                timeconst = ( (double) _time - (omp_get_wtime()-initTCG) );
                CutCG_add_cuts_CG_parallel( ccgcpu, cgraph, LPC_CGCPU, cutP->mip, par->inst, timeconst, par->continuous,  par->maxcuts,  Res_getRound(res), horizon);

                cutTimeCG = (omp_get_wtime()-initTCG);
            //    printf("end cg ID %d\n", omp_get_thread_num()); fflush(stdout);//getchar();
            }
        }

         #pragma omp section
        {
            if(cggpur2)
            {

#ifdef __NVCC__
                //setGpuThread(0);
#endif // __NVCC__
               // int ID = omp_get_thread_num();
               // printf("cg %d MipC_cutCGGPU Rank 2: start on thread %d\n", cggpur2, ID);
               // fflush(stdout);
                double initTCG = omp_get_wtime();
                parameters_ccg *parCCG =(parameters_ccg*)malloc(sizeof(parameters_ccg));

                double timeconst = 0;
                int nRR_cggpu = Inst_nResR(par->inst);

                //Selected rows to GPU
                int horizon=0;
                ccggpur2 = CutCG_create_and_identify_rows_all_conflitcts( cutP->mip, origCols, cutP->mip, par->inst, _time, par->mininstant, par->maxinstant, par->jump,  Res_getRound(res), &horizon);

                timeconst = ( (double) _time - (omp_get_wtime()-initTCG) );
                setParameters_ccg(parCCG, 0); //0 - normal mode, 1- intermediate mode, 2- hard mode


#ifdef __NVCC__
                //start CG
                int cont=0;
                int nC1,nV1;

                nV1 = VDbl_size(ccggpur2->xfElemPP);
                nC1 = VDbl_size(ccggpur2->rowrhs);
                for( int r = 0 ; r < VDbl_size(ccggpur2->rowrhs) ; r++)
                    cont+= VDbl_size(ccggpur2->rowCoef[r]);

                Cut_gpu *ccg;
                Cut_gpu_aux *ccg_aux;
                ccg = fillStruct(ccggpur2,parCCG->precision,nV1,nC1,cont);

                int n_gpu = verifyGpu();
                //printf("verify gpu %d\n", n_gpu);
                fflush(stdout);
                int nBlocks, nThreads, nRuns_temp;
                int i;
                int pos_R1 = 0;
                double aux;
                int n_cuts,nRepeat = 1;

                ccg_aux = AllocationStructCutAux(nC1,nV1,cont);
                returnDimension(&nBlocks,&nThreads, parCCG->nRuns,ccg->numberConstrains);
                n_cuts = ccg->numberConstrains;
                ccg = initial_runGPU(ccg, ccg_aux,parCCG->maxDenominator,parCCG->precision,1,nThreads,nBlocks, nRR_cggpu);
 //               if(n_cuts!=ccg->numberConstrains)
//                    ccg_aux = reallocCut(ccg,ccg_aux);
                n_cuts = ccg->numberConstrains;
                if(parCCG->nRuns < 0.7*(nBlocks*nThreads))
                {
                    nThreads = (parCCG->nRuns/nBlocks) ;
                }
                else
                {
                    aux = (float)(parCCG->nRuns - 0.7*(nThreads*nBlocks))/(float)(0.7*(nThreads*nBlocks));
                    nRepeat += ceil(aux);
                    nThreads = (0.7*(nThreads*nBlocks))/nBlocks;
                }
                for(i=1; i<=nRepeat; i++)
                {
                    timeconst = ( (double) _time - (omp_get_wtime()-initTCG));
                    fflush(stdout);
                    if((nRepeat>1)&&(i == nRepeat))
                    {
                        nThreads = (parCCG->nRuns - (nThreads*nBlocks)*(i-1))/nBlocks;
                    }
                    nRuns_temp = nThreads*nBlocks;
                    if(timeconst>= 1){
                        /*ATUALIZAR TEMPO DENTRO DO METODO timeGPU*/
                        ccg = second_phase_runGPU(ccg, ccg_aux, parCCG->numberMaxConst,nRuns_temp,parCCG->maxDenominator,parCCG->precision, nBlocks, nThreads, &pos_R1, parCCG->nSizeR1,timeconst);
                    }
                }
                returnCutCG(ccg,ccggpur2,parCCG->precision);

                free(parCCG);
                free(ccg);
                free(ccg_aux);
                //end CG
#endif // __NVCC__

                double timegpu = ( (double) timeLeft - (omp_get_wtime()-initT) );
                CutCG_add_cuts_CG_parallel( ccggpur2, cgraph, LPC_CGGPUR2, cutP->mip, par->inst, timegpu, par->continuous,  par->maxcuts,  Res_getRound(res),horizon);

                cutTimeCG = (omp_get_wtime()-initTCG);
             //printf("end cg ID %d\n", omp_get_thread_num()); fflush(stdout);//getchar();
            }
        }
    }



    int cutcr = 0;
    int cutprec = 0;
//    int cutjs = 0;
    int cutcgcpu =0;
    int cutclique = 0;
    int cutoddholes = 0;
    int cutcggpur2 =0;
    int cutcg=0;

    if(cr)
    {
        _time = ( (double) timeLeft - (omp_get_wtime()-initT) );
        double initTCR = omp_get_wtime();
        for(int r = 0; r < Inst_nResR(par->inst); r++)
        {
            //   cutcr += CutP_addSeparatedCuts(cutP,res,cutP->mip,VInt_getPtr(ccr->cutElem[r]),VDbl_getPtr(ccr->cutCoef[r]),VDbl_getPtr(ccr->cutrhs[r]),VDbl_getPtr(ccr->cutviolation[r]),VStr_ptr(ccr->cutname[r]),VInt_getPtr(ccr->cutsense[r]),VInt_getPtr(ccr->cutdominated[r]),par->continuous,nround,LPC_RR,_time);//,ccr->timeseparation,ccr->nAddCuts);
            cutcr += CutP_addSeparatedCuts(cutP,res,cutP->mip,ccr->cutElem[r],ccr->cutCoef[r],ccr->cutrhs[r],ccr->cutviolation[r],ccr->cutname[r],ccr->cutsense[r],ccr->cutdominated[r],par->continuous,nround,LPC_RR,_time);//,ccr->timeseparation,ccr->nAddCuts);
        }
        newCuts += cutcr;
        cutTimeCR += (omp_get_wtime()-initTCR);
        if(VERBOSE) Res_setTCutsTotal(res,Res_getRound(res),LPC_RR,cutTimeCR);
       // printf("\nCR added: %d\n", cutcr); fflush(stdout);

    }

    if(cp)
    {
        _time = ( (double) timeLeft - (omp_get_wtime()-initT) );
        double initTCP = omp_get_wtime();
        cutprec = CutP_addSeparatedCuts(cutP,res,cutP->mip,cpr->cutElem,cpr->cutCoef,cpr->cutrhs,cpr->cutviolation,cpr->cutname,cpr->cutsense,cpr->cutdominated, par->continuous,nround,LPC_PREC,_time);//,cpr->timeseparation,cpr->nAddCuts);
        newCuts += cutprec;
        cutTimeCP += (omp_get_wtime()-initTCP);
        if(VERBOSE) Res_setTCutsTotal(res,Res_getRound(res),LPC_PREC,cutTimeCP);
      //  printf("\nCP added: %d\n", cutprec); fflush(stdout);
    }

  /*  if(cjs)
    {
        _time = ( (double) timeLeft - (omp_get_wtime()-initT) );
        double initTJS = omp_get_wtime();
        cutjs = CutP_addSeparatedCuts(cutP,res,cutP->mip,ccjs->cutElem,ccjs->cutCoef,ccjs->cutrhs,ccjs->cutviolation,ccjs->cutname,ccjs->cutsense,ccjs->cutdominated, par->continuous,nround,LPC_JS,_time);//,cpr->timeseparation,cpr->nAddCuts);
        newCuts += cutjs;
        cutTimeJS += (omp_get_wtime()-initTJS);
        if(VERBOSE) Res_setTCutsTotal(res,Res_getRound(res),LPC_JS,cutTimeJS);
      //  printf("\nCP added: %d\n", cutprec); fflush(stdout);
    }*/

    if(cc)
    {
        _time = ( (double) timeLeft - (omp_get_wtime()-initT) );
        double initTCC = omp_get_wtime();
        cutclique = CutP_addSeparatedCuts(cutP,res,cutP->mip,ccl->cutElem,ccl->cutCoef,ccl->cutrhs,ccl->cutviolation,ccl->cutname,ccl->cutsense,ccl->cutdominated,par->continuous,nround,LPC_CLIQUE,_time);//,ccl->timeseparation,ccl->nAddCuts);ĵ
        newCuts += cutclique;
        cutTimeCC += (omp_get_wtime()-initTCC);
        if(VERBOSE) Res_setTCutsTotal(res,Res_getRound(res),LPC_CLIQUE,cutTimeCC);
        //printf("\nCC added: %d\n", cutclique);  fflush(stdout);
    }

    if(co)
    {
        _time = ( (double) timeLeft - (omp_get_wtime()-initT) );
        double initTCO = omp_get_wtime();
        cutoddholes = CutP_addSeparatedCuts(cutP,res,cutP->mip,coh->cutElem,coh->cutCoef,coh->cutrhs,coh->cutviolation,coh->cutname,coh->cutsense,coh->cutdominated,par->continuous,nround,LPC_ODDHOLES,_time);//,coh->timeseparation,coh->nAddCuts);
        newCuts += cutoddholes;
        cutTimeCO += (omp_get_wtime()-initTCO);
        if(VERBOSE) Res_setTCutsTotal(res,Res_getRound(res),LPC_ODDHOLES,cutTimeCO);
        //printf("\nCO added: %d\n", cutoddholes); fflush(stdout);
    }

    if(cgcpu)
    {

        _time = ( (double) timeLeft - (omp_get_wtime()-initT) );
        double initTCG = omp_get_wtime();
        if(VDbl_size(ccgcpu->cutrhs)>0)
        {
            cutcgcpu = CutP_addSeparatedCuts(cutP,res,cutP->mip,ccgcpu->cutElem,ccgcpu->cutCoef,ccgcpu->cutrhs,ccgcpu->cutviolation,ccgcpu->cutname,ccgcpu->cutsense,ccgcpu->cutdominated,par->continuous,nround,LPC_CGCPU,_time);//,ccgcpu->timeseparation,ccgcpu->nAddCuts);
            newCuts += cutcgcpu;
        }
        CutCG_free(&ccgcpu);

        cutTimeCG += (omp_get_wtime()-initTCG);
        if(VERBOSE) Res_setTCutsTotal(res,Res_getRound(res),LPC_CGCPU,cutTimeCG);
       // printf(" CG added: %d ", cutcgcpu);fflush(stdout);
    }
    if(cggpur2)
    {

        _time = ( (double) timeLeft - (omp_get_wtime()-initT) );
        double initTCG = omp_get_wtime();
        if(VDbl_size(ccggpur2->cutrhs)>0)
        {
            cutcggpur2 = CutP_addSeparatedCuts(cutP,res,cutP->mip,ccggpur2->cutElem,ccggpur2->cutCoef,ccggpur2->cutrhs,ccggpur2->cutviolation,ccggpur2->cutname,ccggpur2->cutsense,ccggpur2->cutdominated,par->continuous,nround,LPC_CGGPUR2,_time);//,ccggpur2->timeseparation,ccggpur2->nAddCuts);
            newCuts += cutcggpur2;
        }
        CutCG_free(&ccggpur2);

        cutTimeCG += (omp_get_wtime()-initTCG);
        if(VERBOSE) Res_setTCutsTotal(res,Res_getRound(res),LPC_CGGPUR2,cutTimeCG);
    }

    if(cg)
    {
        _time = ( (double) timeLeft - (omp_get_wtime()-initT) );
        double initTCG = omp_get_wtime();
        //  printf("add CG\n");fflush(stdout);
        cutcg = CutP_addSeparatedCuts(cutP,res,cutP->mip,ccg->cutElem,ccg->cutCoef,ccg->cutrhs,ccg->cutviolation,ccg->cutname,ccg->cutsense,ccg->cutdominated,par->continuous,nround,LPC_CG,_time);//,ccg->timeseparation,ccg->nAddCuts);
        Res_setUsoUR(res,nround,ccg->contUR);
        Res_setUsoURN(res,nround,ccg->contURN);
        Res_setUsoUM(res,nround,ccg->contUM);
        Res_setUsoUCI(res,nround,ccg->contUCI);
        Res_setUsoUCP(res,nround,ccg->contUCP);
        Res_setTotalUsoUR(res,ccg->contUR);
        Res_setTotalUsoURN(res,ccg->contURN);
        Res_setTotalUsoUM(res,ccg->contUM);
        Res_setTotalUsoUCI(res,ccg->contUCI);
        Res_setTotalUsoUCP(res,ccg->contUCP);
        Res_setNJumps(res,ccg->nJumps);
        newCuts += cutcg;
        cutTimeCG += (omp_get_wtime()-initTCG);
        if(VERBOSE) Res_setTCutsTotal(res,Res_getRound(res),LPC_CG,cutTimeCG);
        CutCG_free(&ccg);
    }

    if(cg)
        fprintf( fp2, "%d %d %d %d %d ", cutcr, cutprec, cutclique, cutoddholes,cutcg);
    if(cggpur2)
        fprintf( fp2, "%d %d %d %d %d ", cutcr, cutprec, cutclique, cutoddholes,cutcggpur2);
    if(cgcpu)
        fprintf( fp2, "%d %d %d %d %d ", cutcr, cutprec, cutclique, cutoddholes,cutcgcpu);

    cutP->totalTimeSeparation += (double) (omp_get_wtime()-initT);

    CutCL_free(&ccl);
    CutC_free(&ccr);
    CutPR_free(&cpr);
//    CutJS_free(&ccjs);
    CutOH_free(&coh);
    return newCuts;

}



/*verify dominance of a value into a vector*/
int CutP_dominance(CutPool *cutP, const int *idxA, double *coefA, double rhs, int sense, int typeA, int sizeA)
{


    /*   double* xf;
       ALLOCATE_VECTOR_INI(xf,double,lp_cols(cutP->mip));
       xf = lp_x(cutP->mip);

       double violA =0;
       printf ("\nCut analyzed A \n");
       for(int cccc = 0; cccc < sizeA ;cccc++){
           char varname1[250];
           lp_col_name(cutP->mip,idxA[cccc],varname1);
           printf ("  %f * %d (%s) x(%f) ", coefA[cccc], idxA[cccc], varname1, xf[idxA[cccc]]);
           violA +=  xf[idxA[cccc]];
       }
       printf (": sense %d rhs %f violation %f\n ", sense, rhs,violA-rhs); //getchar();
    */


    int idxR = 0, dominance = 0;
    if(typeA==LPC_CG ||typeA==LPC_PREC || typeA==LPC_ODDHOLES) return dominance; // typeA==LPC_JS ||

    double r = 0;

    int cont  = lp_cols(cutP->mip);
    int *incElements;
    ALLOCATE_VECTOR_INI(incElements,int,cont);
    int *elements;
    ALLOCATE_VECTOR_INI(elements,int,cont);
    int *posElemts;
    ALLOCATE_VECTOR_INI(posElemts,int,cont);
    double *coefElemA;
    ALLOCATE_VECTOR_INI(coefElemA,double,cont);
    double *coefElemB;
    ALLOCATE_VECTOR_INI(coefElemB,double,cont);

    for(int n = 0 ; n < cutP->nHash ; n++)
    {

        int s = VDbl_size(cutP->hashCuts[n]); //size hash
        idxR = 0;

        // printf("passou -0 n %d idxR %d s %d\n ", n, idxR, s);
        if(s == 0) continue;


        while(idxR < s)
        {

            double typeB = VDbl_get(cutP->hashCuts[n], idxR+1);
            double rhsB = VDbl_get(cutP->hashCuts[n], idxR+2);
//            double senseB = VDbl_get(cutP->hashCuts[n], idxR+3);
            double sizeB = VDbl_get(cutP->hashCuts[n], idxR+4);
            double nCutB = VDbl_get(cutP->hashCuts[n], idxR+5);
            int init =  idxR+5+1;
            int end = init+sizeB-1;
            int initcoef = end+1;
            int endcoef = end+sizeB;

            r = VDbl_get(cutP->hashCuts[n], idxR);

            if(r==-2 )
            {
                idxR = endcoef+1;
                continue; //cut dominated removed
            }

            if( (((typeB == LPC_CLIQUE) && (typeA == LPC_RR)) || ((typeB == LPC_RR) && (typeA == LPC_CLIQUE))) && (rhs <= 1 || rhsB <=1))
            {
                //double viol = 0;


                int nElem = 0;

                //  printf ("\nCut B %f typeB %f \n", r,typeB);
                for(int i =  0 ; i < cont ; i ++)
                {
                    if(i < sizeA)
                    {
                        if(incElements[idxA[i]]==0)
                        {
                            incElements[idxA[i]] = 1;
                            elements[nElem] = idxA[i];
                            posElemts[idxA[i]] = nElem;
                            coefElemA[nElem] = coefA[i];
                            nElem++;
                        }
                        else
                            coefElemA[posElemts[idxA[i]]]= coefA[i];
                    }
                    if(i < sizeB)
                    {
                        double eleme = VDbl_get(cutP->hashCuts[n],init);
                        int elem = (int)eleme;
                        double coefs = VDbl_get(cutP->hashCuts[n],initcoef);
                        char varname2[250];
                        lp_col_name(cutP->mip,elem,varname2);
                        //  printf (" %f * %d (%s) x(%f) ", coefs, elem, varname2, xf[elem]);
                        //   viol += xf[elem];
                        init+=1;
                        initcoef+=1;
                        if(incElements[elem]==0)
                        {
                            incElements[elem] = 1;
                            elements[nElem] = elem;
                            posElemts[elem] = nElem;
                            coefElemB[nElem] = coefs;
                            nElem++;
                        }
                        else
                            coefElemB[posElemts[elem]] = coefs;
                    }

                    if((i>sizeA) && (i>sizeB)) break;
                }
                // printf ("Sense: %f rhs %f violation %f \n ", senseB, rhsB,viol-rhsB); //getchar();

                int contDomA = 0;
                int contDomB = 0;
                dominance = 0;

                for(int cc = 0; cc < nElem ; cc++)
                {
                    // printf("A %f*%d B %f*%d\n", coefElemA[posElemts[elements[cc]]], elements[cc], coefElemB[posElemts[elements[cc]]], elements[cc] );
                    if(coefElemA[posElemts[elements[cc]]] >= coefElemB[posElemts[elements[cc]]]+1e-8)
                        contDomA++;
                    else
                    {
                        if((coefElemB[posElemts[elements[cc]]] >= coefElemA[posElemts[elements[cc]]]+1e-8))
                            contDomB++;
                    }

                    if((contDomA > 0) && (contDomB > 0))
                    {
                        // printf("\nBoth not dominated\n"); //getchar();
                        dominance = -1;
                        break;
                    }
                }

                /*free(incElements);
                free(elements);
                free(posElemts);
                free(coefElemA);
                free(coefElemB);*/
                for(int ie = 0 ; ie <nElem; ie++)
                {
                    incElements[elements[ie]] = 0;
                    coefElemB[posElemts[elements[ie]]] = 0;
                    coefElemA[posElemts[elements[ie]]] = 0;
                    posElemts[elements[ie]] = 0;
                    elements[ie] = 0;
                }

                if(dominance==-1)
                {
                    idxR = endcoef+1;
                    dominance = 0;
                    continue;
                }
                if((contDomA) && (rhs <= rhsB+1e-8))
                {
                    dominance = 2;
                    if( r == 1)
                    {

                        // printf(" INSIDE A domina B "); fflush(stdout);
                        // getchar();
                        // lp_write_lp(cutP->mip,"cutslpdominanceAB.lp"); getchar();
                        const char * namecut = VStr_get(cutP->namesCuts[n],nCutB);
                        int idRemoveR = lp_row_index(cutP->mip,namecut);

                        lp_remove_row(cutP->mip,idRemoveR);
                    }

                    CutP_setIdRow(cutP,n,idxR,-2);
                }
                else
                {
                    if((contDomB) && (rhs+1e-8 >= rhsB))
                    {
                        //printf(" INSIDE B domina A cont A %d contB %d ", contDomA, contDomB);
                        //fflush(stdout);
                        free(incElements);
                        free(elements);
                        free(posElemts);
                        free(coefElemA);
                        free(coefElemB);
                        //getchar();
                        // lp_write_lp(cutP->mip,"cutslpdominanceBA.lp"); getchar();
                        return 1;
                    }
                }
            }
            idxR = endcoef+1;
        }

    }

    free(incElements);
    free(elements);
    free(posElemts);
    free(coefElemA);
    free(coefElemB);
    return dominance; // is not  dominated
}

/*verify dominance between two vector*/
int CutP_dominanceBetweenTwo( LinearProgram *lp, int idxcutA,  int *idxA, double *coefA, double rhs, int sense, int sizeA,  int idxcutB,  int *idxB, double *coefB, double rhsB, int senseB, int sizeB, VecInt *dominatedcuts)
{


    /*  double* xf;
      ALLOCATE_VECTOR_INI(xf,double,lp_cols(lp));
      xf = lp_x(lp);


      double violA =0;
      double violB =0;
      printf ("\nCut analyzed A \n");
      for(int cc = 0; cc < sizeA ;cc++){
          char varname1[250];
          lp_col_name(lp,idxA[cc],varname1);
          printf ("  %f * %d (%s) x(%f) ", coefA[cc], idxA[cc], varname1, xf[idxA[cc]]);
          violA +=  xf[idxA[cc]];
      }
      printf (": sense %d rhs %f violation %f\n ", sense, rhs,violA-rhs); //getchar();

       printf ("\nCut B \n");
      for(int cc = 0; cc < sizeB ;cc++){
          char varname1[250];
          lp_col_name(lp,idxB[cc],varname1);
          printf ("  %f * %d (%s) x(%f) ", coefB[cc], idxB[cc], varname1, xf[idxB[cc]]);
          violB +=  xf[idxB[cc]];
      }
      printf (": sense %d rhs %f violation %f\n ", senseB, rhsB,violB-rhs); //getchar();
    */

    int dominance = -1;

    int cont  = lp_cols(lp);
    int *incElements;
    ALLOCATE_VECTOR_INI(incElements,int,cont);
    int *elements;
    ALLOCATE_VECTOR_INI(elements,int,cont);
    int *posElemts;
    ALLOCATE_VECTOR_INI(posElemts,int,cont);

    int nElem = 0;
    double *coefElemA;
    ALLOCATE_VECTOR_INI(coefElemA,double,cont);
    double *coefElemB;
    ALLOCATE_VECTOR_INI(coefElemB,double,cont);

    for(int i =  0 ; i < cont ; i ++)
    {
        if(i < sizeA)
        {
            if(incElements[idxA[i]]==0)
            {
                incElements[idxA[i]] = 1;
                elements[nElem] = idxA[i];
                posElemts[idxA[i]] = nElem;
                coefElemA[nElem] = coefA[i];
                nElem++;
            }
            else
                coefElemA[posElemts[idxA[i]]]= coefA[i];
        }
        if(i < sizeB)
        {
            if(incElements[idxB[i]]==0)
            {
                incElements[idxB[i]] = 1;
                elements[nElem] = idxB[i];
                posElemts[idxB[i]] = nElem;
                coefElemB[nElem] = coefB[i];
                nElem++;
            }
            else
                coefElemB[posElemts[idxB[i]]]= coefB[i];
        }

        if((i>sizeA) && (i>sizeB)) break;
    }


    int contDomA = 0;
    int contDomB = 0;

    for(int cc = 0; cc < nElem ; cc++)
    {
        if(coefElemA[cc] >= coefElemB[cc]+1e-8)
            contDomA++;
        else
        {
            if((coefElemB[cc] >= coefElemA[cc]+1e-8))
                contDomB++;
        }

        if((contDomA > 0) && (contDomB > 0))
        {
            // printf("\nBoth not dominated\n"); //getchar();
            dominance = 0;
            break;
        }
    }


    if(dominance==0)
    {
        free(incElements);
        free(elements);
        free(posElemts);
        free(coefElemA);
        free(coefElemB);
//         printf("\nBoth not dominated\n"); fflush(stdout); //getchar();
        return dominance; // is not  dominated
    }

    if((contDomA) && (rhs <= rhsB+1e-8))
    {
        dominance = 2;
         // printf("\n A domina B ");
        // fflush(stdout);// getchar();
        //  lp_write_lp(cutP->mip,"cutslpdominanceAB.lp"); getchar();
        VInt_set(dominatedcuts,idxcutB,1);
    }
    else if((contDomB) && (rhs+1e-8 >= rhsB))
    {
        dominance = 1;
      //    printf("\n B domina A ");
        //fflush(stdout); //getchar();
        // lp_write_lp(cutP->mip,"cutslpdominanceBA.lp"); getchar();
        VInt_set(dominatedcuts,idxcutA,1);
    }


    free(incElements);
    free(elements);
    free(posElemts);
    free(coefElemA);
    free(coefElemB);

    return dominance;
}

/*binary search function*/
int binary_search_cut (int elem, const int *idxB,  int l, int r)
{
    assert(idxB);
    //  printf("l %d r %d\n",l,r);fflush(stdout);
    int i = 0;
    i = (l + r)/2;
    if (idxB[i] == elem)
        return i;
    if (l >= r)
        return -1; // Não foi encontrado
    else if (idxB[i] < elem)
        return binary_search_cut(elem, idxB, i+1, r);
    else
        return binary_search_cut(elem, idxB, l, i-1);
}


/*calculates the maximum divisor common  for integers numbers in a vector.
this function uses CutP_maxDivisorCommonRec*/
int CutP_maxDivisorCommonVector(int coefs[], int nElem)
{

    int n = coefs[nElem];
    int mdc = 1;
    while(nElem > 0)
    {
        int m = coefs[nElem-1];
        mdc = CutP_maxDivisorCommonRec( m, n);
        //printf("%d %d: MDC: %d\n", m, n, mdc); fflush(stdout);
        n = mdc;
        nElem--;
    }

    int m = coefs[nElem];
    mdc = CutP_maxDivisorCommonRec( m, n);
    //  printf("%d %d: MDC: %d\n", m, n, mdc); fflush(stdout);
    n = mdc;
//    printf("mdc %d\n",mdc);
    return n;

}

/*calculates the maximum common divisor for two integers.*/
int CutP_maxDivisorCommonRec(int m, int n)
{

    int t = 0;
    m = m < 0 ? -m : m; /* abs(u) */
    n = n < 0 ? -n : n;
    if (m < n)
    {
        t = m;
        m = n;
        n = t;
    }

    if(n==0)
        return m;
    else
    {
        int resto = m % n;
        return CutP_maxDivisorCommonRec( n, resto );
    }

}


/*functions to sort vectors and by int or double*/
void CutP_quick_sort (IntDblPair *a, int n)
{
    if (n < 2)
        return;
    int i = 0, j = 0, p = 0;
    p = a[n / 2].a;
    for (i = 0, j = n - 1;; i++, j--)
    {
        int t = 0;
        while (a[i].a < p)
            i++;
        while (p < a[j].a)
            j--;
        if (i >= j)
            break;
        t = a[i].a;
        a[i].a = a[j].a;
        a[j].a = t;
    }
    CutP_quick_sort(a, i);
    CutP_quick_sort(a + i, n - i);
}

void CutP_quick_sort_vec ( int *idx, double *coe, int n)
{
    if (n < 2)
        return;
    int i = 0, j = 0, p = 0;
    p = idx[n/2];
    for (i = 0, j = n - 1;; i++, j--)
    {
        int t = 0;
        double co = 0;
        while (idx[i] < p)
            i++;
        while (p < idx[j])
            j--;
        if (i >= j)
            break;
        t = idx[i];
        idx[i]=idx[j];
        idx[j]=t;

        co = coe[i];
        coe[i]=coe[j];
        coe[j]=co;
    }
    CutP_quick_sort_vec(idx, coe, i);
    CutP_quick_sort_vec(idx+i, coe+i, n - i);
}

void CutP_quick_sort_vec_by_double ( int *idx, double *coe, int n)
{
    if (n < 2)
        return;
    int i = 0, j = 0;
    double p = 0;
    p = coe[n/2];
    for (i = 0, j = n - 1;; i++, j--)
    {
        double va = 0.0;
        int id = 0;
        while (coe[i] < p)
            i++;
        while (p < coe[j])
            j--;
        if (i >= j)
            break;
        va = coe[i];
        coe[i]=coe[j];
        coe[j]=va;

        id = idx[i];
        idx[i]=idx[j];
        idx[j]=id;
    }
    CutP_quick_sort_vec_by_double(idx, coe, i);
    CutP_quick_sort_vec_by_double(idx+i, coe+i, n - i);
}


 /*Lift a cut by tightening the value on the right side based on a MIP with constraints about RR, NR, Job, Conflicts to all variables of the cut */
double CutP_model_lift_cgraph( const Instance *inst, const CGraph *cgraph, const char* cutname, int cutnelem, int *cutElem, double *cutCoef, double rhs, LinearProgram *lp, double timeLeft, int horizon){

    assert(lp);
    assert(inst);

    double tinit = omp_get_wtime();



    char prefix[STR_SIZE];
    int id[5];
    char name[STR_SIZE];

    int *idx;
    ALLOCATE_VECTOR_INI(idx,int,lp_cols(lp));
    double *objbin;
    ALLOCATE_VECTOR_INI(objbin,double,cutnelem);
    int **diffzerorr;
    ALLOCATE_VECTOR_INI(diffzerorr,int*,Inst_nResR(inst));
    VecInt ***idxR;
    ALLOCATE_VECTOR_INI(idxR,VecInt**,Inst_nResR(inst));
    VecInt **idxNR;
    ALLOCATE_VECTOR_INI(idxNR,VecInt*,Inst_nResN(inst));

    VecDbl ***coefR;
    ALLOCATE_VECTOR_INI(coefR,VecDbl**,Inst_nResR(inst));
    for(int res = 0 ; res < Inst_nResR(inst); res++){
        ALLOCATE_VECTOR_INI(coefR[res],VecDbl*,horizon);
        ALLOCATE_VECTOR_INI(idxR[res],VecInt*,horizon);
        ALLOCATE_VECTOR_INI(diffzerorr[res],int,horizon);
        for(int td = 0; td < horizon ;td++){
            coefR[res][td] = VDbl_create();
            idxR[res][td] = VInt_create();
        }
    }


    VecDbl **coefNR;
    ALLOCATE_VECTOR_INI(coefNR,VecDbl*,Inst_nResN(inst));
    for(int res = 0 ; res < Inst_nResN(inst); res++){
        coefNR[res] = VDbl_create();
        idxNR[res] = VInt_create();
    }
    int *diffzeronr;
    ALLOCATE_VECTOR_INI(diffzeronr,int,Inst_nResN(inst));

    parseName( cutname, prefix, id );
   // printf("%s:", cutname);
//    int time = id[1];

    VecStr *namesbin = VStr_create(STR_SIZE);
    double *x = lp_x(lp);
    double *xf;
    ALLOCATE_VECTOR_INI(xf,double,cutnelem);
    double *xfcghaph;
    ALLOCATE_VECTOR_INI(xfcghaph, double, lp_cols(lp));

    for(int e = 0 ; e < cutnelem ;e++){
        int j,m,t;
        lp_col_name(lp,cutElem[e],name);
        parseName( name, prefix, id );
        j = id[0];
        m = id[1];
        const Job *job = Inst_job(inst,j);
        const Mode *mode = Job_mode(job,m);
        t = id[2];
        //xf[e]= x[j][m][t];
        xf[e]= x[cutElem[e]];

        char vnamex[STR_SIZE];
        sprintf( vnamex, "x(%d,%d,%d)", j, m,t);
        idx[cutElem[e]] = e;
        VStr_pushBack( namesbin, vnamex );
        objbin[e] = cutCoef[e]*-1;
        //printf(" %f * %s + ", objbin[e], name);

       for(int res = 0 ; res < Mode_nResR(mode); res++){
            for(int ttd = t ; ttd< t+Mode_duration(mode)+1 ; ttd++){
                int r = Mode_idxResR(mode,res);
                VDbl_pushBack(coefR[r][ttd],Mode_useResR(mode,res));
                VInt_pushBack(idxR[r][ttd],e);
                diffzerorr[r][ttd]++;
            }
        }
        for(int res = 0 ; res < Mode_nResN(mode); res++){
            int nr = Mode_idxResN(mode,res);
            VDbl_pushBack(coefNR[nr],Mode_useResN(mode,res));
            VInt_pushBack(idxNR[nr],e);
           diffzeronr[nr]++;
        }

    }
  //  printf(" rhs %f \n", rhs); fflush(stdout);

    LinearProgram *mipCutm;
    mipCutm = lp_create();
    double _time = (double)timeLeft - (omp_get_wtime()-tinit);
    lp_set_max_seconds(mipCutm, _time);
    lp_set_print_messages(mipCutm,0);

    lp_add_bin_cols( mipCutm, cutnelem, objbin , VStr_ptr(namesbin) );

    for(int res = 0 ; res < Inst_nResR(inst); res++){
        for(int ttd = 0 ; ttd< horizon ; ttd++){
            if(diffzerorr[res][ttd]==0) continue;
            char nameR[STR_SIZE];
            sprintf( nameR, "cResR(%d,%d)", res,ttd);
            lp_add_row(mipCutm,diffzerorr[res][ttd],VInt_getPtr(idxR[res][ttd]),VDbl_getPtr(coefR[res][ttd]),nameR,'L',Inst_capResR(inst,res));
            //lp_add_row(mipCutm,cutnelem,idx,coefR[res],nameR,'L',Inst_capResR(inst,res));
        }
    }
    for(int res = 0 ; res < Inst_nResN(inst); res++){
        if(diffzeronr[res]==0) continue;
       char nameR[STR_SIZE];
       sprintf( nameR, "cResNR(%d)", res);
       lp_add_row(mipCutm,diffzeronr[res],VInt_getPtr(idxNR[res]),VDbl_getPtr(coefNR[res]),nameR,'L',Inst_capResN(inst,res));
       //lp_add_row(mipCutm,cutnelem,idx,coefNR[res],nameR,'L',Inst_capResN(inst,res));
    }


    int **varInConfs;
    ALLOCATE_VECTOR_INI(varInConfs, int*, lp_cols(lp));
    for(int il = 0 ; il < lp_cols(lp) ; il++)
        ALLOCATE_VECTOR_INI(varInConfs[il], int, lp_cols(lp));

    int nProj = Inst_nProjects(inst);
    int *nVProj;
    ALLOCATE_VECTOR_INI(nVProj,int,nProj);
    IntTriple **vectorProjects;
    ALLOCATE_VECTOR(vectorProjects,IntTriple*,nProj);
    for(int i = 0 ; i < nProj ; i++)
        ALLOCATE_VECTOR_INI(vectorProjects[i],IntTriple,cutnelem);


    VecInt **sameJob;
    ALLOCATE_VECTOR_INI(sameJob,VecInt*,Inst_nJobs(inst));
    for(int i = 0; i < Inst_nJobs(inst) ; i++)
        sameJob[i] = VInt_create();

   // int *jobsvar;
   // ALLOCATE_VECTOR_INI(jobsvar,int,Inst_nJobs(inst));

    for(int l = 0; l < cutnelem; l++) {

        int j,m,t;
        lp_col_name(lp,cutElem[l],name);
        parseName( name, prefix, id );
        j = id[0];
        m = id[1];
        const Job *job = Inst_job(inst,j);
//        const Mode *mode = Job_mode(job,m);
        t = id[2];

       // jobsvar[j] = 1;
        VInt_pushBack(sameJob[j], l);
        int plp = Job_project(job);

        vectorProjects[plp][nVProj[plp]].idx = cutElem[l];
        vectorProjects[plp][nVProj[plp]].j = j;
        vectorProjects[plp][nVProj[plp]].m = m;
        vectorProjects[plp][nVProj[plp]].t = t;
        vectorProjects[plp][nVProj[plp]].value = x[cutElem[l]];
        nVProj[plp]++;

    }



    for(int i = 0; i < Inst_nJobs(inst) ; i++) {
        int size = VInt_size(sameJob[i]);
        if(size<=0) continue;
        //if(jobsvar[] != 1) continue;
//        int nme=0;
        double coefsj[size];
        // printf("\nElem for job %d\n", i);
        for(int e = 0 ; e < size ; e++ )
            coefsj[e] = 1.0;
        double rhsm = 1;
        char nameR[STR_SIZE];
        sprintf( nameR, "cSelmode(%d)", i);
        lp_add_row(mipCutm,size,VInt_getPtr(sameJob[i]),coefsj,nameR,'L',rhsm);

    }

    for(int p = 0 ; p < nProj; p++) {
        //  if(cont >= lp_rows(lp)) break;
        for(int elem = 0 ; elem < nVProj[p] ; elem++) {
            //if(cont >= lp_rows(lp)) break;
            int l = vectorProjects[p][elem].idx;
            int j = vectorProjects[p][elem].j, m = vectorProjects[p][elem].m, t =  vectorProjects[p][elem].t;
        //        printf("\n E (%d,%d,%d) %d ", j, m, t, l);
            const Job *job = Inst_job(inst,j);
            const Mode *mode = Job_mode(job,m);
            if(Mode_duration(mode)==0)continue;

            for(int elem2 = 0 ; elem2 < nVProj[p] ; elem2++) {
                //  if(cont >= lp_rows(lp)) break;
                if(elem == elem2 ) continue;
                int j2 = vectorProjects[p][elem2].j, m2 = vectorProjects[p][elem2].m, t2 = vectorProjects[p][elem2].t;
                if(j==j2 && m==m2 && t ==t2) continue;

                if(j==j2) {
                    int l2 = vectorProjects[p][elem2].idx;
                    if(varInConfs[l][l2]==1) continue;
                  //  VInt_pushBack(conflicts[l],l2);
                    xfcghaph[l] = vectorProjects[p][elem].value;
                    xfcghaph[l2] = vectorProjects[p][elem2].value;
              //          printf("Same Job (%d,%d,%d) %d ", j2, m2, t2, l2);
                    varInConfs[l][l2] = 1;
                    continue;
                }


                const Job *job2 = Inst_job(inst,j2);
                const Mode *mode2 = Job_mode(job2,m2);
                if(Mode_duration(mode2)==0)continue;

                int durationE = Mode_duration(mode);
                int timeEndPred = t+durationE;

                if( Job_hasIndSucc(job,j2) ) {
                    int value = Inst_getMaxDIJ(inst,j,j2)- Job_minDuration(job);
                    int winTime = value ;
                    if(timeEndPred>t2 || ((t2-timeEndPred) < winTime)) {
                        //     printf("timeEndPred = t %d + durationE %d \n", t, durationE );
                   //     printf("\ntimeEndPred %d > t2 %d || ((t2 %d - timeEndPred %d) < winTime %d)", timeEndPred,t2,t2,timeEndPred, winTime);
                        int l2 = vectorProjects[p][elem2].idx;
                        if(varInConfs[l][l2]==1) continue;
                        //     printf("prec D(%d,%d,%d) %d ", j2, m2, t2, l2);
                        xfcghaph[l] = vectorProjects[p][elem].value;
                        xfcghaph[l2] = vectorProjects[p][elem2].value;
                       // VInt_pushBack(conflicts[l],l2);
                        varInConfs[l][l2] = 1;
                        continue;
                    }
                }
            }

            for(int p2 = 0 ; p2 < nProj ; p2++) {
                //  if(cont >= lp_rows(lp)) break;
                if(p==p2) continue;
                for(int elem2 = 0 ; elem2 < nVProj[p2] ; elem2++) {
                    //  if(cont >= lp_rows(lp)) break;
                    int j2 = vectorProjects[p2][elem2].j, m2 = vectorProjects[p2][elem2].m, t2 = vectorProjects[p2][elem2].t;
                    if(j==j2 && m==m2 && t ==t2) continue;

                    const Job *job2 = Inst_job(inst,j2);
                    const Mode *mode2 = Job_mode(job2,m2);
                    if(Mode_duration(mode2)==0)continue;

                    const Project *project = Inst_project(inst,p);
                    int cp1 = Project_releaseDate(project)+Project_criticalPath(project) - Inst_getMaxDIM(inst,job,Mode_index(mode));
                    const Project *project2 = Inst_project(inst,Job_project(job2));
                    int cp2 = Project_releaseDate(project2)+Project_criticalPath(project2) - Inst_getMaxDIM(inst,job2,Mode_index(mode2));

                    int part1 = t - cp1 < 0 ? 0 :
                                t - cp1 ;
                    int part2 = t2 - cp2 < 0 ? 0 :
                                t2 - cp2 ;
                    if(part1 + part2 > Inst_getSumTPD(inst)) {
                        //printf("\npart1 %d (t%d-cp1 %d) + part2 %d (t2%d-cp2 %d) > Inst_getSumTPD(inst) %d  releaseDatep1 %d + criticalPathp1 %d  DJMj1 %d releaseDatep2 %d + criticalPathp %d DIM2 %d", part1, t, cp1, part2, t2, cp2, Inst_getSumTPD(inst), Project_releaseDate(project),Project_criticalPath(project) , Inst_getMaxDIM(inst,job,Mode_index(mode)), Project_releaseDate(project2),Project_criticalPath(project2),Inst_getMaxDIM(inst,job2,Mode_index(mode2)));
                        int l2 = vectorProjects[p2][elem2].idx;
                        if(varInConfs[l][l2]==1) continue;
                       // VInt_pushBack(conflicts[l],l2);
                        xfcghaph[l] = vectorProjects[p][elem].value;
                        xfcghaph[l2] = vectorProjects[p2][elem2].value;
                         //        printf("tpd p(%d,%d,%d) %d ", j2, m2, t2, l2);
                        varInConfs[l][l2] = 1;

                    }
                }
            }
        }
    }


    CliqueSeparation *sep = clq_sep_create(cgraph);
    clq_sep_set_maxcliques(sep,  100);
    clq_sep_set_max_it_bk(sep, 1000);
    clq_sep_set_extend_method(sep,0);

    _time = ( (double) timeLeft - (omp_get_wtime()-tinit) );
    if(_time <= clq_sep_get_max_time_bk(sep))
        clq_sep_set_max_time_bk(sep,_time);

    clq_sep_separate( sep, xfcghaph);

    const CliqueSet *clqSet = clq_sep_get_cliques( sep );
    int nCliques = clq_get_number_of_cliques(clqSet);


    int totalConf = 0;
    for(int nc = 0 ; nc < nCliques ; nc++) {
        int * isetelements = clq_set_clique_elements(clqSet,nc);
        int size = clq_set_clique_size(clqSet,nc);
        int idxconv[size];
        for(int nci = 0 ; nci < size ; nci++){
            idxconv[nci] = idx[isetelements[nci]];
        }
        double coef[size];
        FILL(coef,0,size,1);
        char nameR[STR_SIZE];
        sprintf( nameR, "cConf(%d)", nc);
        lp_add_row(mipCutm,size,idxconv,coef,nameR,'L',1.0);
        totalConf++;
    }
    clq_sep_free(&sep);


//lp_set_concurrentMIP(mipCutm,1);
    //lp_set_method(mipCutm,4);
        //lp_set_seed(mipCutm,100000);
    int st = lp_optimize(mipCutm);

    if(st==LP_OPTIMAL)
    {


//        for(int s = 0 ; s <)
        if (lp_obj_value(mipCutm)*-1 < rhs){
            double newrhs =lp_obj_value(mipCutm)*-1;
      //  printf("creating mipCutm RR: end\n");
      //  fflush(stdout);
       // lp_write_sol(mipCutm,"mipCutm.sol");
      //  lp_write_lp(mipCutm,"mipCutm.lp");// getchar();

           // FILE *fp = fopen( "liftcglog.txt", "a" );
           // if ( fp == NULL ) {
             //   printf( "File was not opened. : liftcglog.txt\n");
            //    exit( 0 );
            //}
            //fprintf( fp, " %s ; %d ; %f ; %f ; %d ; \n", cutname, cutnelem, rhs, newrhs, totalConf);
           // printf( " %s ; %d ; %f ; %f ; %d ; \n", cutname, cutnelem, rhs, newrhs, totalConf);
           // fclose( fp );

            for(int i = 0; i < Inst_nJobs(inst) ; i++)
                VInt_free(&sameJob[i]);
            free(sameJob);

            for(int il = 0 ; il < lp_cols(lp) ; il++)
                free(varInConfs[il]);
            free(varInConfs);

            free(idx);
            free(objbin);
            for(int res = 0 ; res < Inst_nResR(inst); res++){
                for(int td = 0; td < horizon ;td++){
                    VDbl_free(&coefR[res][td]);
                    VInt_free(&idxR[res][td]);
                }
                free(coefR[res]);
                free(idxR[res]);
            }

            for(int res = 0 ; res < Inst_nResN(inst); res++){
                VDbl_free(&coefNR[res]);
                VInt_free(&idxNR[res]);
            }
            free(idxR);
            free(idxNR);
            free(coefR);
            free(coefNR);
            free(diffzeronr);
            free(diffzerorr);
            free(xf);
            free(xfcghaph);

            VStr_free(&namesbin);
            lp_free(&mipCutm);
//            getchar();
            return newrhs;
        }
    }


    for(int i = 0; i < Inst_nJobs(inst) ; i++)
        VInt_free(&sameJob[i]);
    free(sameJob);

            for(int il = 0 ; il < lp_cols(lp) ; il++)
                free(varInConfs[il]);
            free(varInConfs);

    free(idx);
    free(objbin);
    for(int res = 0 ; res < Inst_nResR(inst); res++){
        for(int td = 0; td < horizon ;td++){
            VDbl_free(&coefR[res][td]);
            VInt_free(&idxR[res][td]);
        }
        free(coefR[res]);
        free(idxR[res]);
    }


    for(int res = 0 ; res < Inst_nResN(inst); res++){
        VDbl_free(&coefNR[res]);
        VInt_free(&idxNR[res]);
    }
    free(idxR);
    free(idxNR);
    free(coefR);
    free(coefNR);
    free(diffzeronr);
    free(diffzerorr);
    free(xf);
    free(xfcghaph);

    VStr_free(&namesbin);
    lp_free(&mipCutm);
    return -1;
}

/*Lift a cut by tightening the value on the right side based on a MIP with constraints about RR, NR, Job, Conflicts to all variables of the cut */
double CutP_model_lift( const Instance *inst, const char* cutname, int cutnelem, int *cutElem, double *cutCoef, double rhs, LinearProgram *lp, double timeLeft, int horizon){

    assert(lp);
    assert(inst);

    double tinit = omp_get_wtime();



    char prefix[STR_SIZE];
    int id[5];
    char name[STR_SIZE];

    int *idx;
    ALLOCATE_VECTOR_INI(idx,int,cutnelem);
    double *objbin;
    ALLOCATE_VECTOR_INI(objbin,double,cutnelem);
    int **diffzerorr;
    ALLOCATE_VECTOR_INI(diffzerorr,int*,Inst_nResR(inst));
    VecInt ***idxR;
    ALLOCATE_VECTOR_INI(idxR,VecInt**,Inst_nResR(inst));
    VecInt **idxNR;
    ALLOCATE_VECTOR_INI(idxNR,VecInt*,Inst_nResN(inst));

    VecDbl ***coefR;
    ALLOCATE_VECTOR_INI(coefR,VecDbl**,Inst_nResR(inst));
    for(int res = 0 ; res < Inst_nResR(inst); res++){
        ALLOCATE_VECTOR_INI(coefR[res],VecDbl*,horizon);
        ALLOCATE_VECTOR_INI(idxR[res],VecInt*,horizon);
        ALLOCATE_VECTOR_INI(diffzerorr[res],int,horizon);
        for(int td = 0; td < horizon ;td++){
            coefR[res][td] = VDbl_create();
            idxR[res][td] = VInt_create();
        }
    }


    VecDbl **coefNR;
    ALLOCATE_VECTOR_INI(coefNR,VecDbl*,Inst_nResN(inst));
    for(int res = 0 ; res < Inst_nResN(inst); res++){
        coefNR[res] = VDbl_create();
        idxNR[res] = VInt_create();
    }
    int *diffzeronr;
    ALLOCATE_VECTOR_INI(diffzeronr,int,Inst_nResN(inst));

    parseName( cutname, prefix, id );
   // printf("%s:", cutname);
//    int time = id[1];

    VecStr *namesbin = VStr_create(STR_SIZE);
    double *x = lp_x(lp);
    double *xf;
    ALLOCATE_VECTOR_INI(xf,double,cutnelem);

    for(int e = 0 ; e < cutnelem ;e++){
        int j,m,t;
        lp_col_name(lp,cutElem[e],name);
        parseName( name, prefix, id );
        j = id[0];
        m = id[1];
        const Job *job = Inst_job(inst,j);
        const Mode *mode = Job_mode(job,m);
        t = id[2];
        //xf[e]= x[j][m][t];
        xf[e]= x[cutElem[e]];

        char vnamex[STR_SIZE];
        sprintf( vnamex, "x(%d,%d,%d)", j, m,t);
        idx[e] = e;
        VStr_pushBack( namesbin, vnamex );
        objbin[e] = cutCoef[e]*-1;
        //printf(" %f * %s + ", objbin[e], name);

       for(int res = 0 ; res < Mode_nResR(mode); res++){
            for(int ttd = t ; ttd< t+Mode_duration(mode)+1 ; ttd++){
                int r = Mode_idxResR(mode,res);
                VDbl_pushBack(coefR[r][ttd],Mode_useResR(mode,res));
                VInt_pushBack(idxR[r][ttd],e);
                diffzerorr[r][ttd]++;
            }
        }
        for(int res = 0 ; res < Mode_nResN(mode); res++){
            int nr = Mode_idxResN(mode,res);
            VDbl_pushBack(coefNR[nr],Mode_useResN(mode,res));
            VInt_pushBack(idxNR[nr],e);
           diffzeronr[nr]++;
        }

    }
  //  printf(" rhs %f \n", rhs); fflush(stdout);

    LinearProgram *mipCutm;
    mipCutm = lp_create();
    double _time = (double)timeLeft - (omp_get_wtime()-tinit);
    lp_set_max_seconds(mipCutm, _time);
    lp_set_print_messages(mipCutm,0);

    lp_add_bin_cols( mipCutm, cutnelem, objbin , VStr_ptr(namesbin) );

    for(int res = 0 ; res < Inst_nResR(inst); res++){
        for(int ttd = 0 ; ttd< horizon ; ttd++){
            if(diffzerorr[res][ttd]==0) continue;
            char nameR[STR_SIZE];
            sprintf( nameR, "cResR(%d,%d)", res,ttd);
            lp_add_row(mipCutm,diffzerorr[res][ttd],VInt_getPtr(idxR[res][ttd]),VDbl_getPtr(coefR[res][ttd]),nameR,'L',Inst_capResR(inst,res));
            //lp_add_row(mipCutm,cutnelem,idx,coefR[res],nameR,'L',Inst_capResR(inst,res));
        }
    }
    for(int res = 0 ; res < Inst_nResN(inst); res++){
        if(diffzeronr[res]==0) continue;
       char nameR[STR_SIZE];
       sprintf( nameR, "cResNR(%d)", res);
       lp_add_row(mipCutm,diffzeronr[res],VInt_getPtr(idxNR[res]),VDbl_getPtr(coefNR[res]),nameR,'L',Inst_capResN(inst,res));
       //lp_add_row(mipCutm,cutnelem,idx,coefNR[res],nameR,'L',Inst_capResN(inst,res));
    }


    VecInt** conflicts;
    ALLOCATE_VECTOR(conflicts,VecInt*,cutnelem);

    int **varInConfs;
    ALLOCATE_VECTOR_INI(varInConfs, int*, cutnelem);
    int nProj = Inst_nProjects(inst);
    int *nVProj;
    ALLOCATE_VECTOR_INI(nVProj,int,nProj);
    IntTriple **vectorProjects;
    ALLOCATE_VECTOR(vectorProjects,IntTriple*,nProj);
    for(int i = 0 ; i < nProj ; i++)
        ALLOCATE_VECTOR_INI(vectorProjects[i],IntTriple,cutnelem);


    VecInt **sameJob;
    ALLOCATE_VECTOR_INI(sameJob,VecInt*,Inst_nJobs(inst));
    for(int i = 0; i < Inst_nJobs(inst) ; i++)
        sameJob[i] = VInt_create();

   // int *jobsvar;
   // ALLOCATE_VECTOR_INI(jobsvar,int,Inst_nJobs(inst));

    for(int l = 0; l < cutnelem; l++) {
        conflicts[l] = VInt_create();
        ALLOCATE_VECTOR_INI(varInConfs[l], int, cutnelem);

        int j,m,t;
        lp_col_name(lp,cutElem[l],name);
        parseName( name, prefix, id );
        j = id[0];
        m = id[1];
        const Job *job = Inst_job(inst,j);
//        const Mode *mode = Job_mode(job,m);
        t = id[2];

       // jobsvar[j] = 1;
        VInt_pushBack(sameJob[j], l);
        int plp = Job_project(job);

        vectorProjects[plp][nVProj[plp]].idx = l;
        vectorProjects[plp][nVProj[plp]].j = j;
        vectorProjects[plp][nVProj[plp]].m = m;
        vectorProjects[plp][nVProj[plp]].t = t;
        vectorProjects[plp][nVProj[plp]].value = xf[l];
        nVProj[plp]++;

    }



    for(int i = 0; i < Inst_nJobs(inst) ; i++) {
        int size = VInt_size(sameJob[i]);
        if(size<=0) continue;
        //if(jobsvar[] != 1) continue;
//        int nme=0;
        double coefsj[size];
        // printf("\nElem for job %d\n", i);
        for(int e = 0 ; e < size ; e++ )
            coefsj[e] = 1.0;
        double rhsm = 1;
        char nameR[STR_SIZE];
        sprintf( nameR, "cSelmode(%d)", i);
        lp_add_row(mipCutm,size,VInt_getPtr(sameJob[i]),coefsj,nameR,'L',rhsm);

    }

    for(int p = 0 ; p < nProj; p++) {
        //  if(cont >= lp_rows(lp)) break;
        for(int elem = 0 ; elem < nVProj[p] ; elem++) {
            //if(cont >= lp_rows(lp)) break;
            int l = vectorProjects[p][elem].idx;
            int j = vectorProjects[p][elem].j, m = vectorProjects[p][elem].m, t =  vectorProjects[p][elem].t;
        //        printf("\n E (%d,%d,%d) %d ", j, m, t, l);
            const Job *job = Inst_job(inst,j);
            const Mode *mode = Job_mode(job,m);
            if(Mode_duration(mode)==0)continue;

            for(int elem2 = 0 ; elem2 < nVProj[p] ; elem2++) {
                //  if(cont >= lp_rows(lp)) break;
                if(elem == elem2 ) continue;
                int j2 = vectorProjects[p][elem2].j, m2 = vectorProjects[p][elem2].m, t2 = vectorProjects[p][elem2].t;
                if(j==j2 && m==m2 && t ==t2) continue;

                if(j==j2) {
                    int l2 = vectorProjects[p][elem2].idx;
                    if(varInConfs[l][l2]==1) continue;
                    VInt_pushBack(conflicts[l],l2);
              //          printf("Same Job (%d,%d,%d) %d ", j2, m2, t2, l2);
                    varInConfs[l][l2] = 1;
                    continue;
                }


                const Job *job2 = Inst_job(inst,j2);
                const Mode *mode2 = Job_mode(job2,m2);
                if(Mode_duration(mode2)==0)continue;

                int durationE = Mode_duration(mode);
                int timeEndPred = t+durationE;

                if( Job_hasIndSucc(job,j2) ) {
                    int value = Inst_getMaxDIJ(inst,j,j2)- Job_minDuration(job);
                    int winTime = value ;
                    if(timeEndPred>t2 || ((t2-timeEndPred) < winTime)) {
                        //     printf("timeEndPred = t %d + durationE %d \n", t, durationE );
                   //     printf("\ntimeEndPred %d > t2 %d || ((t2 %d - timeEndPred %d) < winTime %d)", timeEndPred,t2,t2,timeEndPred, winTime);
                        int l2 = vectorProjects[p][elem2].idx;
                        if(varInConfs[l][l2]==1) continue;
                        //     printf("prec D(%d,%d,%d) %d ", j2, m2, t2, l2);
                        VInt_pushBack(conflicts[l],l2);
                        varInConfs[l][l2] = 1;
                        continue;
                    }
                }
            }

            for(int p2 = 0 ; p2 < nProj ; p2++) {
                //  if(cont >= lp_rows(lp)) break;
                if(p==p2) continue;
                for(int elem2 = 0 ; elem2 < nVProj[p2] ; elem2++) {
                    //  if(cont >= lp_rows(lp)) break;
                    int j2 = vectorProjects[p2][elem2].j, m2 = vectorProjects[p2][elem2].m, t2 = vectorProjects[p2][elem2].t;
                    if(j==j2 && m==m2 && t ==t2) continue;

                    const Job *job2 = Inst_job(inst,j2);
                    const Mode *mode2 = Job_mode(job2,m2);
                    if(Mode_duration(mode2)==0)continue;

                    const Project *project = Inst_project(inst,p);
                    int cp1 = Project_releaseDate(project)+Project_criticalPath(project) - Inst_getMaxDIM(inst,job,Mode_index(mode));
                    const Project *project2 = Inst_project(inst,Job_project(job2));
                    int cp2 = Project_releaseDate(project2)+Project_criticalPath(project2) - Inst_getMaxDIM(inst,job2,Mode_index(mode2));

                    int part1 = t - cp1 < 0 ? 0 :
                                t - cp1 ;
                    int part2 = t2 - cp2 < 0 ? 0 :
                                t2 - cp2 ;
                    if(part1 + part2 > Inst_getSumTPD(inst)) {
                        //printf("\npart1 %d (t%d-cp1 %d) + part2 %d (t2%d-cp2 %d) > Inst_getSumTPD(inst) %d  releaseDatep1 %d + criticalPathp1 %d  DJMj1 %d releaseDatep2 %d + criticalPathp %d DIM2 %d", part1, t, cp1, part2, t2, cp2, Inst_getSumTPD(inst), Project_releaseDate(project),Project_criticalPath(project) , Inst_getMaxDIM(inst,job,Mode_index(mode)), Project_releaseDate(project2),Project_criticalPath(project2),Inst_getMaxDIM(inst,job2,Mode_index(mode2)));
                        int l2 = vectorProjects[p2][elem2].idx;
                        if(varInConfs[l][l2]==1) continue;
                        VInt_pushBack(conflicts[l],l2);
                         //        printf("tpd p(%d,%d,%d) %d ", j2, m2, t2, l2);
                        varInConfs[l][l2] = 1;

                    }
                }
            }
        }
    }


    int **conf;
    ALLOCATE_VECTOR(conf, int*, cutnelem);

    int *nConf;
    ALLOCATE_VECTOR_INI(nConf, int, cutnelem);
    int weight[cutnelem];
    //printf("find conflicts clique: end\n");fflush(stdout);
    for(int l = 0 ; l < cutnelem ; l++) {
        int sizeConflictsJ = VInt_size(conflicts[l]);
        nConf[l] = sizeConflictsJ;
        // printf("\nconflicts[%d] size %d: ", l, sizeConflictsJ);
        if(sizeConflictsJ != 0) {
            int *confJ = VInt_getPtr(conflicts[l]);
            double *coef;
            ALLOCATE_VECTOR_INI(coef, double, sizeConflictsJ);
            FILL(coef,0,sizeConflictsJ,1.0);
            ALLOCATE_VECTOR(conf[l], int, sizeConflictsJ);
            for(int o = 0 ; o <  sizeConflictsJ; o++) {

                //  printf("confJ[o]%d", confJ[o]);
                conf[l][o]= confJ[o];
                //printf("conf[%d][%d] %d ", l,o, conf[l][o]);
                //getchar();
            }
            CutP_quick_sort_vec(conf[l],coef,sizeConflictsJ);
            if(xf[l]<=0.00001)
                weight[l] = 1;
            else
                weight[l] = xf[l]*1000.0;
            free(coef);
        } else
            weight[l] = 0.0;
    }

    _time = ( (double) timeLeft - (omp_get_wtime()-tinit) );
    CGraph *cgraph = build_cgraph_conflicts(conf, nConf, cutnelem, _time);

    //cgraph_create_induced_subgraph(cgraph,)

    cgraph_set_weight( cgraph, weight );
    BronKerbosch *bk =  bk_create(cgraph);
    _time = ( (double) timeLeft - (omp_get_wtime()-tinit) );
//    bk_set_max_it(bk, 100);
    //bk_set_max_cliques(bk, 5);
  //  if(_time<1.5)
    //    bk_set_timelimit(bk, _time);
   // else
     //   bk_set_timelimit(bk, 1.0);
  //  bk_set_min_weight(bk, weight);
    //bk_set_timelimit(bk, _time);
    bk_run(bk);



    CliqueSet *clqSet = bk_get_clq_set(bk);
    int nCliqueSet =  clq_get_number_of_cliques(clqSet);

    int totalConf = 0;
    for(int nc = 0 ; nc < nCliqueSet ; nc++) {
        int* iset = clq_set_clique_elements(clqSet, nc);
        int size = clq_set_clique_size(clqSet, nc);
        double coef[size];
        FILL(coef,0,size,1);
        char nameR[STR_SIZE];
        sprintf( nameR, "cConf(%d)", nc);
        lp_add_row(mipCutm,size,iset,coef,nameR,'L',1.0);
        totalConf++;
    }


//lp_set_concurrentMIP(mipCutm,1);
    //lp_set_method(mipCutm,4);
        //lp_set_seed(mipCutm,100000);
    int st = lp_optimize(mipCutm);

    if(st==LP_OPTIMAL)
    {


//        for(int s = 0 ; s <)
        if (lp_obj_value(mipCutm)*-1 < rhs){
            double newrhs =lp_obj_value(mipCutm)*-1;
      //  printf("creating mipCutm RR: end\n");
      //  fflush(stdout);
       // lp_write_sol(mipCutm,"mipCutm.sol");
      //  lp_write_lp(mipCutm,"mipCutm.lp");// getchar();

           // FILE *fp = fopen( "liftcglog.txt", "a" );
           // if ( fp == NULL ) {
             //   printf( "File was not opened. : liftcglog.txt\n");
            //    exit( 0 );
            //}
            //fprintf( fp, " %s ; %d ; %f ; %f ; %d ; \n", cutname, cutnelem, rhs, newrhs, totalConf);
           // printf( " %s ; %d ; %f ; %f ; %d ; \n", cutname, cutnelem, rhs, newrhs, totalConf);
           // fclose( fp );

            for(int i = 0; i < Inst_nJobs(inst) ; i++)
                VInt_free(&sameJob[i]);
            free(sameJob);

            for(int l = 0 ; l < cutnelem ; l++) {
                if( nConf[l]!=0)  free(conf[l]);
                VInt_free(&conflicts[l]);
                free(varInConfs[l]);
            }

            free(conf);
            free(nConf);
            free(conflicts);
            free(varInConfs);

            free(idx);
            free(objbin);
            for(int res = 0 ; res < Inst_nResR(inst); res++){
                for(int td = 0; td < horizon ;td++){
                    VDbl_free(&coefR[res][td]);
                    VInt_free(&idxR[res][td]);
                }
                free(coefR[res]);
                free(idxR[res]);
                free(diffzerorr[res]);
            }

            for(int res = 0 ; res < Inst_nResN(inst); res++){
                VDbl_free(&coefNR[res]);
                VInt_free(&idxNR[res]);
            }
            free(idxR);
            free(idxNR);
            free(coefR);
            free(coefNR);
            free(diffzeronr);
            free(diffzerorr);
            free(xf);


            for(int p = 0 ; p < nProj ; p++)
                free(vectorProjects[p]);
            free(vectorProjects);
            free(nVProj);
            cgraph_free(&cgraph);
            VStr_free(&namesbin);
            bk_free(bk);
            lp_free(&mipCutm);
//            getchar();
            return newrhs;
        }
    }


    for(int i = 0; i < Inst_nJobs(inst) ; i++)
        VInt_free(&sameJob[i]);
    free(sameJob);

    for(int l = 0 ; l < cutnelem ; l++) {
        if( nConf[l]!=0)  free(conf[l]);
        VInt_free(&conflicts[l]);
        free(varInConfs[l]);
    }

    free(conf);
    free(nConf);
    free(conflicts);
    free(varInConfs);

    free(idx);
    free(objbin);
    for(int res = 0 ; res < Inst_nResR(inst); res++){
        for(int td = 0; td < horizon ;td++){
            VDbl_free(&coefR[res][td]);
            VInt_free(&idxR[res][td]);
        }
        free(diffzerorr[res]);
        free(coefR[res]);
        free(idxR[res]);
    }


    for(int res = 0 ; res < Inst_nResN(inst); res++){
        VDbl_free(&coefNR[res]);
        VInt_free(&idxNR[res]);
    }
    free(idxR);
    free(idxNR);
    free(coefR);
    free(coefNR);
    free(diffzeronr);
    free(diffzerorr);
    free(xf);


    for(int p = 0 ; p < nProj ; p++)
        free(vectorProjects[p]);
    free(vectorProjects);
    free(nVProj);
    cgraph_free(&cgraph);
    VStr_free(&namesbin);
    bk_free(bk);
    lp_free(&mipCutm);
    return -1;
}

/*get and set functions*/
VecDbl** CutP_getHC( CutPool *cutP )
{
    assert(cutP!=NULL);
    return cutP->hashCuts;
}

double CutP_getTotalTimeSeparation( CutPool *cutP )
{
    assert(cutP!=NULL);
    return cutP->totalTimeSeparation;
}

void CutP_setNameRow( CutPool *cutP, int key, int idCut, const char *value)
{
    VStr_set(cutP->namesCuts[key], idCut, value);
}

void CutP_setIdRow( CutPool *cutP, int key, int idRow, int value)
{
    VDbl_set(cutP->hashCuts[key], idRow, value);
}

IdxCoefCut *CutP_getCut(CutPool *cutP, int key, int idxR)
{
    assert(cutP);

    IdxCoefCut *cut;
    ALLOCATE_INI(cut,IdxCoefCut)

    double s = VDbl_size(cutP->hashCuts[key]); //size hash

    double r = VDbl_get(cutP->hashCuts[key], idxR);
    //    int tc = VDbl_get(cutP->hashCuts[key], idxR+1);
    double nrhs = VDbl_get(cutP->hashCuts[key], idxR+2);
    double nsense = VDbl_get(cutP->hashCuts[key], idxR+3);
    double sc = VDbl_get(cutP->hashCuts[key], idxR+4);
    double ncut = VDbl_get(cutP->hashCuts[key], idxR+5);
    int init =  idxR+5+1;
    int end = init+sc-1;
    int initcoef = end+1;
    int endcoef = end+sc;
    //   printf("\nKey: %d Cut Size: %d \n", key, s);

    // printf("1) idxR %d, r %f, idxTC %d, tc %f, idxRhs %d, nrhs %f, idxSense %d,  nsense %f, idxSC %d, sc %f, init %d, initcoef %d, end %d, endcoef %d, s %d \n", idxR, r, idxTC, tc, idxRhs, nrhs, idxSense, nsense, idxSC, sc, init, initcoef, end, endcoef, s);

    ALLOCATE_VECTOR(cut->idx, int, lp_cols(cutP->mip));
    ALLOCATE_VECTOR(cut->coef, double, lp_cols(cutP->mip));

    int o = 0;
    //    double value = 0;

    int ihcoef = initcoef;

    for(int ihce = init ; ihce < init+sc && ihcoef < initcoef+sc && ihcoef < s; ihce++, ihcoef++)
    {
        // printf("\nihce %d == init %d +sc %f -1 && ihcoef %d == initcoef %d +sc %f -1 %d \n", ihce, init, sc, ihcoef,initcoef,sc);

        double elem = VDbl_get(cutP->hashCuts[key],ihce);
        double coefs = VDbl_get(cutP->hashCuts[key],ihcoef);

        cut->idx[o] = (int) elem;
        cut->coef[o] = coefs;
        //   printf(" %f * %d ,", coefs, elem);
        o++;

        if(ihce == init+sc-1 && ihcoef == initcoef+sc-1)
        {
            cut->rhs = nrhs;
            cut->sense = (int) nsense;
            cut->sizec = (int) sc;
            cut->row = (int) r;
            cut->endc =  endcoef;
            cut->ncut = (int) ncut;
            // printf("nElem %f, sense %f, rhs %f r(activated) %f, ncut %f\n", sc, nsense, nrhs, r, ncut );
        }
    }

    return cut;

}


/*free memory*/
void CutP_free( CutPool **_cutP )
{
    CutPool *cutP = *_cutP;

    for ( int i=0 ; i<cutP->nHash ; ++i )
    {
        VDbl_free( &cutP->hashCuts[i] );
        VStr_free( &cutP->namesCuts[i] );
    }
    free( cutP->hashCuts );
    free( cutP->namesCuts );
    free( cutP->nCuts );
    //    free( cutP->staticValues );
    free( cutP );
    *_cutP = NULL;
}

