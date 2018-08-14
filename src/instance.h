/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef INSTANCE_H_DEFINED
#define INSTANCE_H_DEFINED

typedef struct _Instance Instance;
typedef struct _Job Job;
typedef struct _Project Project;
typedef struct _Mode Mode;

char Mode_isFeasible( const Instance *inst, const Mode *mode );

/* loads instance from a file */
Instance *Inst_create( const char *dir, const char *instance );
Instance *Inst_read( char ** argv, int argc );
/* number of jobs */
int Inst_nJobs( const Instance *inst );

/* maximum number of modes per job */
int Inst_nMaxModes( const Instance *inst );

/* number of projects */
int Inst_nProjects( const Instance *inst );

/* returns job idxJob */
const Job *Inst_job( const Instance *inst, int idxJob );

const Job *Project_job( const Project *p, int idxJob );

/* total number of renewable resources, i.e. global and local,
 * considering that local resources are have unique indexes */
int Inst_nResR( const Instance *inst );

/* number of global renewable resources */
int Inst_nResRGlobal( const Instance *inst );

/* capacity of a renewable resource */
int Inst_capResR( const Instance *inst, int i );

int Inst_getSizePath(const Instance *inst, int i);
int Inst_getValuePosPath(const Instance *inst, int i, int pos);
int Inst_getMaxDIM(const Instance *inst, const Job *job,int m);
int Inst_getMaxDIMJM(const Instance *inst, const Job *job, int m,  int j2, int m2);
/* capacity of a non-renewable resource */
int Inst_capResN( const Instance *inst, int i );

/* total number of non-renewable resources (all local),
 * all have unique indexes */
int Inst_nResN( const Instance *inst );

/* returns the index where the non renewable resource of project starts*/
int Inst_idxResNProj( const Instance *inst, int idxProj );

/* returns project idxProject */
const Project *Inst_project( const Instance *inst, int idxProject );

/* returns the index of project*/
int Project_index( const Project *p );

/* returns the release date of project */
int Project_releaseDate( const Project *p );

/* returns the dueDate of project */
int Project_dueDate( const Project *p );

/* returns the tardcost of project */
int Project_tardCost( const Project *p );

/* returns the MPM-Time of project */
int Project_mpmTime( const Project *p );

/* returns the criticalPath of project */
int Project_criticalPath( const Project *p );

/* calculate and return the criticalPath of project */
int Project_calc_criticalPath(Instance *inst);

/* returns the number of jobs on project */
int Project_nJobs( const Project *p );

/* returns the index of the first project job in the instance */
int Project_idxFirstJob( const Project *p );

/*returns the number of successors of job*/
int Job_nSucc( const Job *job );

/* to access the successor idxSucc of job */
int Job_succ( const Job *job, int idxSucc );

/*returns the number of predecessors of job*/
int Job_nPred( const Job *job );

/* to access the predecessor idxPred of job */
int Job_pred( const Job *job, int idxPred );

/*returns the min duration of job*/
int Job_minDuration( const Job *job );

/*returns the max duration of job*/
int Job_maxDuration( const Job *job );

/*returns the early start time of job*/
int Job_est( const Job *job );

/*returns the index of job on instance*/
int Job_index( const Job *job );

/*returns the index of job on project*/
int Job_idxOnProject( const Job *job );

/* returns the number of modes for a job */
int Job_nModes( const Job *job );

/* returns the index of modes on instance*/
int Job_idxMode( const Job *job, int idMode );

/* returns true if job have the successor idxJob*/
int Job_hasSucc( const Instance *inst, const Job *job, int idxJob );

/* returns true if job have the predecessor idxJob*/
int Job_hasPred( const Instance *inst, const Job *job, int idxJob );

/* returns true if job have the indirect predecessor idxJob*/
int Job_hasIndPred(  const Job *job, int idxJob );

/* returns true if job have the indirect successor idxJob*/
int Job_hasIndSucc(  const Job *job, int idxJob );

/* returns the index of project of job */
int Job_project( const Job *job );

/* returns the number of infeasible modes for a job */
int Job_nInfeasModes( const Job *job );

/* returns mode information for a job */
const Mode *Job_mode( const Job *job, int idxMode );

/* duration of a mode */
int Mode_duration( const Mode *mode );

/* number of renewable resources used by this mode */
int Mode_nResR( const Mode *mode );

/* indexes of renewable resources used
 * with this mode */
int Mode_idxResR( const Mode *mode, int i );

/* use of renewable resources */
int Mode_useResR( const Mode *mode, int i );

/* number of non renewable resources used by a mode */
int Mode_nResN( const Mode *mode );

/* to access indexes of local non-renewable resources used in a mode */
int Mode_idxResN( const Mode *mode, int i );

int Mode_idxResROnMode( const Instance *inst, const Mode *mode, int i );

/* to access the use of local non-renewable resources used in a mode */
int Mode_useResN( const Mode *mode, int i );

int somaMinDurationPathsByJob( const Instance *inst, int j, int nSuccs);
int somaMinDurationPathsByJobAndInter(const Instance *inst, int j, int nSuccs, int inter, int durInter);
int Inst_computeCompPathsByJobAndInter( const Instance *inst, int i, int m, int d, int s, int m2, int ds, int dsmin );
int Inst_computeCompPathsByJob( const Instance *inst, int i, int m,  int d);
//int Job_minDurationPath( const Job *job, int m );
void floydWarshallMax(const Instance *inst, int **D, int nJobs);
void maxDistanceByModes(const Instance *inst, int **D, int *** DJM);
int Inst_getMaxDIJ(const Instance *inst, int i, int j);

/* returns index of mode*/
int Mode_index( const Mode *mode );

/* set sumTPD: max Delay*/
void Inst_setSumTPD(Instance *inst, int value);

/* returns sumTPD: max Delay*/
int Inst_getSumTPD(const Instance *inst);
int cmp_int_pair_b( const void *v1, const void *v2 );


void Inst_setSumTMS(Instance *inst, const int value);
int Inst_getSumTMS(const Instance *inst);

//int Job_durationModeSort(const Job *job, int m);
//int Job_idxModeSort( const Job *job, int m);

int Inst_getMaxDIJM(const Instance *inst, int i, int j, int m);
/* frees memory used by a instance */
void Inst_free( Instance **_inst );

/*print the entire contents of an instance*/
void Inst_print( const Instance *inst );

/* checks if the file exists*/
int Inst_fileExists(char fileName[]);

/* check if there are any job infeasible with all others to any resource*/
void Inst_jobsInfeasible(const Instance *inst);

/* prints the contents of a project*/
void Project_print( const Instance *inst, int p );

/* print the entire contents of a job*/
void Job_print( const Instance *inst, int i );

/* print the entire contents of a mode*/
void Mode_print( const Instance *inst, const Job *job, int i, int nHid );

#endif /* INSTANCE_H_DEFINED */

