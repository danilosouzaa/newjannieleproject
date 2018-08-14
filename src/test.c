#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "macros.h"
#include "neighborhood.h"
#include "test.h"

struct _Test {

    int nNeigh;
    int *improveNeigh;
    int *visitNeigh;
    Cost *improveFO;
    double *timeNeigh;
    double totalTime;
    Cost currentFO;

    double currentTime;
    int currentNeigh;

    int T;
    int SAmax;
    double alpha;

    int nStages;
    int lastN;
    int it;

    const struct _Instance *inst;

    long double **I;
    long double **EQ;
    long double **TI;
    long double **TE;
    long double **TIV;
    long double **FI;
    long double **FE;

    long double **normFI;
    long double **normFE;

    long double **PFI;
    long double **PFE;

    long double *maxFI;

    long double *minFI;

    long double *maxFE;
    long double *minFE;

    long double *intervalI;
    long double *intervalE;
};

void Test_free( Test **_test )
{

    Test *test = *_test;

    free( test->improveNeigh );
    free( test->visitNeigh );
    free( test->timeNeigh );
    free( test->improveFO );

    for(int i = 0; i < test->nNeigh ; i++)   {
        free( test->EQ[i] );
        free( test->I[i] );
        free( test->TI[i] );
        free( test->TE[i] );
        free( test->TIV[i] );
        free( test->FI[i] );
        free( test->FE[i] );
        free( test->normFI[i] );
        free( test->normFE[i] );
        free( test->PFE[i] );
        free( test->PFI[i] );
    }

    free( test->EQ );
    free( test->I );
    free( test->TI );
    free( test->TE );
    free( test->TIV );
    free( test->FI );
    free( test->FE );
    free( test->normFI );
    free( test->normFE );
    free( test->PFE );
    free( test->PFI );

    free( test->intervalE);
    free( test->intervalI);
    free( test->maxFI );
    free( test->maxFE );
    free( test->minFI );
    free( test->minFE );

    free( test );

    *_test = NULL;

}

Test *Test_create(int nNeigh, const Instance* inst)
{

    assert( nNeigh > 0 );

    Test *test;



    ALLOCATE_INI( test, Test );
    ALLOCATE_VECTOR_INI( test->improveNeigh, int, nNeigh );
    ALLOCATE_VECTOR_INI( test->visitNeigh, int, nNeigh );
    ALLOCATE_VECTOR_INI( test->improveFO, Cost, nNeigh );
    ALLOCATE_VECTOR_INI( test->timeNeigh, double, nNeigh );

    test->inst = inst;

    test->totalTime=0;
    test->nNeigh=nNeigh;
    test->currentNeigh = 1;
    test->currentFO=0;

    test->nStages = 2;
    test->it = 10000;

    long double **I = (long double **) calloc (test->nNeigh, sizeof(long double *));
    for(int i = 0; i < test->nNeigh ; i++)
        I[i] = (long double*) calloc (test->nStages, sizeof(long double));
    test->I = I;

    long double **EQ = (long double **) calloc (test->nNeigh, sizeof(long double *));
    for(int i = 0; i < test->nNeigh ; i++)
        EQ[i] = (long double*) calloc (test->nStages, sizeof(long double));
    test->EQ = EQ;

    long double **TI = (long double **) calloc (test->nNeigh, sizeof(long double *));
    for(int i = 0; i < test->nNeigh ; i++)
        TI[i] = (long double*) calloc (test->nStages, sizeof(long double));
    test->TI = TI;


    long double **TE = (long double **) calloc (test->nNeigh, sizeof(long double *));
    for(int i = 0; i < test->nNeigh ; i++)
        TE[i] = (long double*) calloc (test->nStages, sizeof(long double));
    test->TE = TE;

    long double **TIV = (long double **) calloc (test->nNeigh, sizeof(long double *));
    for(int i = 0; i < test->nNeigh ; i++)
        TIV[i] = (long double*) calloc (test->nStages, sizeof(long double));
    test->TIV = TIV;


    long double **FI= (long double **) calloc (test->nNeigh, sizeof(long double *));
    for(int i = 0; i < test->nNeigh ; i++)
        FI[i] = (long double*) calloc (test->nStages, sizeof(long double));
    test->FI = FI;

    long double **FE= (long double **) calloc (test->nNeigh, sizeof(long double *));
    for(int i = 0; i < test->nNeigh ; i++)
        FE[i] = (long double*) calloc (test->nStages, sizeof(long double));
    test->FE = FE;


    long double **normFI= (long double **) calloc (test->nNeigh, sizeof(long double *));
    for(int i = 0; i < test->nNeigh ; i++)
        normFI[i] = (long double*) calloc (test->nStages, sizeof(long double));
    test->normFI = normFI;

    long double **normFE= (long double **) calloc (test->nNeigh, sizeof(long double *));
    for(int i = 0; i < test->nNeigh ; i++)
        normFE[i] = (long double*) calloc (test->nStages, sizeof(long double));
    test->normFE = normFE;

    long double **PFI= (long double **) calloc (test->nNeigh, sizeof(long double *));
    for(int i = 0; i < test->nNeigh ; i++)
        PFI[i] = (long double*) calloc (test->nStages, sizeof(long double));
    test->PFI = PFI;

    long double **PFE= (long double **) calloc (test->nNeigh, sizeof(long double *));
    for(int i = 0; i < test->nNeigh ; i++)
        PFE[i] = (long double*) calloc (test->nStages, sizeof(long double));
    test->PFE = PFE;

    long double *maxFI = (long double *) calloc (test->nStages, sizeof(long double));
    test->maxFI =maxFI;

    long double *intervalI = ( long double *) calloc (test->nStages, sizeof(long double));
    test->intervalI = intervalI;

    long double *maxFE = (long double *) calloc (test->nStages, sizeof(long double));
    test->maxFE =maxFE;

    long double *intervalE = (long double *) calloc (test->nStages, sizeof(long double));
    test->intervalE = intervalE;


    long double *minFI= (long double *) calloc (test->nStages, sizeof(long double));
    long double *minFE= (long double *) calloc (test->nStages, sizeof(long double));

    for(int i = 0 ; i < test->nStages; i++) {
        minFI[i] = (long double)INT_MAX_M;
        minFE[i] = (long double)INT_MAX_M;
    }


    test->minFI =minFI;
    test->minFE =minFE;



    return test;
}

void Test_readAnalisysNeigh(Test *test, int sw, int it, double psw, char *argv)
{

    char name[256];
    int timeP;
    float** sumProb;
    int** cont;
    float prob;



    ALLOCATE_VECTOR(sumProb, float*, test->nNeigh );
    for(int i = 0; i < test->nNeigh ; i++)
        ALLOCATE_VECTOR_INI( sumProb[i], float, 31);


    ALLOCATE_VECTOR(cont, int*, test->nNeigh );
    for(int i = 0; i < test->nNeigh ; i++)
        ALLOCATE_VECTOR_INI( cont[i], int, 31);



    /*
        for(int c = 1 ; c <=3 ; c++){
        for(int iiii = 1 ; iiii<=10 ; iiii++ ){

        char insta1[256]="";
        char number[256]="";

        sprintf(number, "%d", iiii);

        if(c==1){
        strcat (insta1, "A-");
        strcat (insta1, number);
        strcat (insta1,".txt");
        }if(c==2){
        strcat (insta1, "B-");
        strcat (insta1, number);
        strcat (insta1,".txt");
        }if(c==3){
        strcat (insta1, "X-");
        strcat (insta1, number);
        strcat (insta1,".txt");
        }

        printf("insta1 %s ******\n", insta1);
    */
    char fil[256] = "probabilitiesOnline_";
    char bufferII[256];
    char side[256];
    char penalty[256];
    char instance[256];
    sprintf(instance, "%s", argv);//insta1);
    strcat (fil,instance);
    strcat (fil,"_");
    sprintf(bufferII, "%d", it);
    strcat (fil,bufferII);
    strcat (fil,"_");
    sprintf(side, "%d", sw);
    strcat (fil,side);
    strcat (fil,"_");
    sprintf(penalty, "%f", psw);
    strcat (fil,penalty);
    strcat (fil,".txt");

    printf("fil %s ******\n", fil);
    getchar();

    FILE* ap = fopen(fil, "r");

    int idx=0, idxT=0;

    if(ap == NULL)
        printf("Erro, nao foi possivel abrir o arquivo\n");
    else {
        while( (fscanf(ap,"%d %s %f\n", &timeP, name, &prob))!=EOF ) {

            if(strcmp("ISJ", name)==0)
                idx = seqInvert-1;
            if(strcmp("OJ", name)==0)
                idx = seqShiftJob-1;
            if(strcmp("STJ", name)==0)
                idx = seqSwapJob-1;
            if(strcmp("OP", name)==0)
                idx = seqShiftProj-1;
            if(strcmp("SCTP", name)==0)
                idx = seqSwapProj-1;
            if(strcmp("CPP", name)==0)
                idx = seqCompactProj-1;
            if(strcmp("COM", name)==0)
                idx = changeOneMode-1;
            if(strcmp("CTM", name)==0)
                idx = changeTwoMode-1;
            if(strcmp("CThM", name)==0)
                idx = changeThreeMode-1;
            if(strcmp("CFM", name)==0)
                idx = changeFourMode-1;
            if(strcmp("SSJW", name)==0)
                idx = seqSwapJobFILS-1;
            if(strcmp("SSIW", name)==0)
                idx = seqInsertJobFILS-1;
            if(strcmp("SPE", name)==0)
                idx = seqCompOnExtrem-1;
            if(strcmp("CSP", name)==0)
                idx = seqMoveProj-1;


            if(timeP == 0)   {
                idxT = 0;
                cont[idx][idxT]++;
            }
            if(timeP == 10)  {
                idxT = 1;
                cont[idx][idxT]++;
            }
            if(timeP == 20)  {
                idxT = 2;
                cont[idx][idxT]++;
            }
            if(timeP == 30)  {
                idxT = 3;
                cont[idx][idxT]++;
            }
            if(timeP == 40)  {
                idxT = 4;
                cont[idx][idxT]++;
            }
            if(timeP == 50)  {
                idxT = 5;
                cont[idx][idxT]++;
            }
            if(timeP == 60)  {
                idxT = 6;
                cont[idx][idxT]++;
            }
            if(timeP == 70)  {
                idxT = 7;
                cont[idx][idxT]++;
            }
            if(timeP == 80)  {
                idxT = 8;
                cont[idx][idxT]++;
            }
            if(timeP == 90)  {
                idxT = 9;
                cont[idx][idxT]++;
            }
            if(timeP == 100) {
                idxT = 10;
                cont[idx][idxT]++;
            }
            if(timeP == 110) {
                idxT = 11;
                cont[idx][idxT]++;
            }
            if(timeP == 120) {
                idxT = 12;
                cont[idx][idxT]++;
            }
            if(timeP == 130) {
                idxT = 13;
                cont[idx][idxT]++;
            }
            if(timeP == 140) {
                idxT = 14;
                cont[idx][idxT]++;
            }
            if(timeP == 150) {
                idxT = 15;
                cont[idx][idxT]++;
            }
            if(timeP == 160) {
                idxT = 16;
                cont[idx][idxT]++;
            }
            if(timeP == 170) {
                idxT = 17;
                cont[idx][idxT]++;
            }
            if(timeP == 180) {
                idxT = 18;
                cont[idx][idxT]++;
            }
            if(timeP == 190) {
                idxT = 19;
                cont[idx][idxT]++;
            }
            if(timeP == 200) {
                idxT = 20;
                cont[idx][idxT]++;
            }
            if(timeP == 210) {
                idxT = 21;
                cont[idx][idxT]++;
            }
            if(timeP == 220) {
                idxT = 22;
                cont[idx][idxT]++;
            }
            if(timeP == 230) {
                idxT = 23;
                cont[idx][idxT]++;
            }
            if(timeP == 240) {
                idxT = 24;
                cont[idx][idxT]++;
            }
            if(timeP == 250) {
                idxT = 25;
                cont[idx][idxT]++;
            }
            if(timeP == 260) {
                idxT = 26;
                cont[idx][idxT]++;
            }
            if(timeP == 270) {
                idxT = 27;
                cont[idx][idxT]++;
            }
            if(timeP == 280) {
                idxT = 28;
                cont[idx][idxT]++;
            }
            if(timeP == 290) {
                idxT = 29;
                cont[idx][idxT]++;
            }
            if(timeP == 300) {
                idxT = 30;
                cont[idx][idxT]++;
            }

            sumProb[idx][idxT] += prob ;
            //cont[idx][timeP]++;
            //sumProb[idx][timeP] += prob ;
            //printf("%d %s %f\n", timeP, name, prob);

        }
    }
    fclose(ap);
    //}
    //}

    for(int n = 0 ; n <test->nNeigh ; n++) {
        for(int i = 0 ; i <=31 ; i++) {
            char* sigla;
            switch(n) {
                case seqInvert-1:
                    sigla = "ISJ";
                    break;
                case seqShiftJob-1:
                    sigla = "OJ";
                    break;
                case seqSwapJob-1:
                    sigla = "STJ";
                    break;
                case seqShiftProj-1:
                    sigla = "OP";
                    break;
                case seqSwapProj-1:
                    sigla = "SCTP";
                    break;
                case seqCompactProj-1:
                    sigla = "CPP";
                    break;
                case changeOneMode-1:
                    sigla = "COM";
                    break;
                case changeTwoMode-1:
                    sigla = "CTM";
                    break;
                case changeThreeMode-1:
                    sigla = "CThM";
                    break;
                case changeFourMode-1:
                    sigla = "CFM";
                    break;
                case seqSwapJobFILS-1:
                    sigla = "SSJW";
                    break;
                case seqInsertJobFILS-1:
                    sigla = "SSIW";
                    break;
                case seqCompOnExtrem-1:
                    sigla = "SPE";
                    break;
                case seqMoveProj-1:
                    sigla = "CSP";
                    break;
            }
            int iT;
            switch(i) {
                case 0:
                    iT = 0;
                    break;
                case 1:
                    iT = 10;
                    break;
                case 2:
                    iT = 20;
                    break;
                case 3:
                    iT = 30;
                    break;
                case 4:
                    iT = 40;
                    break;
                case 5:
                    iT = 50;
                    break;
                case 6:
                    iT = 60;
                    break;
                case 7:
                    iT = 70;
                    break;
                case 8:
                    iT = 80;
                    break;
                case 9:
                    iT = 90;
                    break;
                case 10:
                    iT = 100;
                    break;
                case 11:
                    iT = 110;
                    break;
                case 12:
                    iT = 120;
                    break;
                case 13:
                    iT = 130;
                    break;
                case 14:
                    iT = 140;
                    break;
                case 15:
                    iT = 150;
                    break;
                case 16:
                    iT = 160;
                    break;
                case 17:
                    iT = 170;
                    break;
                case 18:
                    iT = 180;
                    break;
                case 19:
                    iT = 190;
                    break;
                case 20:
                    iT = 200;
                    break;
                case 21:
                    iT = 210;
                    break;
                case 22:
                    iT = 220;
                    break;
                case 23:
                    iT = 230;
                    break;
                case 24:
                    iT = 240;
                    break;
                case 25:
                    iT = 250;
                    break;
                case 26:
                    iT = 260;
                    break;
                case 27:
                    iT = 270;
                    break;
                case 28:
                    iT = 280;
                    break;
                case 29:
                    iT = 290;
                    break;
                case 30:
                    iT = 300;
                    break;

            }

            //            iT=i;
            char fil[256] = "chartProbabilities_";
            char bufferII[5];
            sprintf(bufferII, "%s", sigla);
            strcat (fil,bufferII);
            strcat (fil,".dat");

            FILE* app = fopen(fil, "a");
            float probabilidades =  (float) sumProb[n][i] / (float)cont[n][i];
            printf("\n%d %.8f", iT, probabilidades );
            fprintf(app,"\n%d %.8f", iT, probabilidades );
            fclose(app);
        }
    }




}

void Test_readAnalisysNeighToInstance(Test *test, int sw, int it, double psw, char *argv)
{

    char name[256];
    float timeP;
    float prob;

    char fil[256] = "probabilitiesOnlineInstance_";
    char bufferII[256];
    char side[256];
    char penalty[256];
    char instance[256];
    sprintf(instance, "%s", argv);//insta1);
    strcat (fil,instance);
    strcat (fil,"_");
    sprintf(bufferII, "%d", it);
    strcat (fil,bufferII);
    strcat (fil,"_");
    sprintf(side, "%d", sw);
    strcat (fil,side);
    strcat (fil,"_");
    sprintf(penalty, "%f", psw);
    strcat (fil,penalty);
    strcat (fil,".txt");

    printf("fil %s ******\n", fil);
    getchar();

    FILE* ap = fopen(fil, "r");

    while( (fscanf(ap,"%f %s %f\n", &timeP, name, &prob))!=EOF ) {
        char* sigla;

        if(strcmp("ISJ", name)==0) {
            sigla = "ISJ";
            char fil[256] = "chartProbabilitiesInstance_";
            char bufferII[5];
            sprintf(bufferII, "%s", sigla);
            strcat (fil,bufferII);
            strcat (fil,".dat");

            FILE* app = fopen(fil, "a");
            printf("\n%f %.8f", timeP, prob );
            fprintf(app,"\n%f %.8f", timeP, prob );
            fclose(app);

        }
        if(strcmp("OJ", name)==0) {
            sigla = "OJ";
            char fil[256] = "chartProbabilitiesInstance_";
            char bufferII[5];
            sprintf(bufferII, "%s", sigla);
            strcat (fil,bufferII);
            strcat (fil,".dat");

            FILE* app = fopen(fil, "a");
            printf("\n%f %.8f", timeP, prob );
            fprintf(app,"\n%f %.8f", timeP, prob );
            fclose(app);

        }
        if(strcmp("STJ", name)==0) {
            sigla = "STJ";
            char fil[256] = "chartProbabilitiesInstance_";
            char bufferII[5];
            sprintf(bufferII, "%s", sigla);
            strcat (fil,bufferII);
            strcat (fil,".dat");

            FILE* app = fopen(fil, "a");
            printf("\n%f %.8f", timeP, prob );
            fprintf(app,"\n%f %.8f", timeP, prob );
            fclose(app);

        }
        if(strcmp("OP", name)==0) {
            sigla = "OP";
            char fil[256] = "chartProbabilitiesInstance_";
            char bufferII[5];
            sprintf(bufferII, "%s", sigla);
            strcat (fil,bufferII);
            strcat (fil,".dat");

            FILE* app = fopen(fil, "a");
            printf("\n%f %.8f", timeP, prob );
            fprintf(app,"\n%f %.8f", timeP, prob );
            fclose(app);

        }
        if(strcmp("SCTP", name)==0) {
            sigla = "SCTP";
            char fil[256] = "chartProbabilitiesInstance_";
            char bufferII[5];
            sprintf(bufferII, "%s", sigla);
            strcat (fil,bufferII);
            strcat (fil,".dat");

            FILE* app = fopen(fil, "a");
            printf("\n%f %.8f", timeP, prob );
            fprintf(app,"\n%f %.8f", timeP, prob );
            fclose(app);

        }
        if(strcmp("CPP", name)==0) {
            sigla = "CPP";
            char fil[256] = "chartProbabilitiesInstance_";
            char bufferII[5];
            sprintf(bufferII, "%s", sigla);
            strcat (fil,bufferII);
            strcat (fil,".dat");

            FILE* app = fopen(fil, "a");
            printf("\n%f %.8f", timeP, prob );
            fprintf(app,"\n%f %.8f", timeP, prob );
            fclose(app);

        }
        if(strcmp("COM", name)==0) {
            sigla = "COM";
            char fil[256] = "chartProbabilitiesInstance_";
            char bufferII[5];
            sprintf(bufferII, "%s", sigla);
            strcat (fil,bufferII);
            strcat (fil,".dat");

            FILE* app = fopen(fil, "a");
            printf("\n%f %.8f", timeP, prob );
            fprintf(app,"\n%f %.8f", timeP, prob );
            fclose(app);

        }
        if(strcmp("CTM", name)==0) {
            sigla = "CTM";
            char fil[256] = "chartProbabilitiesInstance_";
            char bufferII[5];
            sprintf(bufferII, "%s", sigla);
            strcat (fil,bufferII);
            strcat (fil,".dat");

            FILE* app = fopen(fil, "a");
            printf("\n%f %.8f", timeP, prob );
            fprintf(app,"\n%f %.8f", timeP, prob );
            fclose(app);
        }
        if(strcmp("CThM", name)==0) {
            sigla = "CThM";
            char fil[256] = "chartProbabilitiesInstance_";
            char bufferII[5];
            sprintf(bufferII, "%s", sigla);
            strcat (fil,bufferII);
            strcat (fil,".dat");

            FILE* app = fopen(fil, "a");
            printf("\n%f %.8f", timeP, prob );
            fprintf(app,"\n%f %.8f", timeP, prob );
            fclose(app);
        }
        if(strcmp("CFM", name)==0) {
            sigla = "CFM";
            char fil[256] = "chartProbabilitiesInstance_";
            char bufferII[5];
            sprintf(bufferII, "%s", sigla);
            strcat (fil,bufferII);
            strcat (fil,".dat");

            FILE* app = fopen(fil, "a");
            printf("\n%f %.8f", timeP, prob );
            fprintf(app,"\n%f %.8f", timeP, prob );
            fclose(app);
        }
        if(strcmp("SSJW", name)==0) {
            sigla = "SSJW";
            char fil[256] = "chartProbabilitiesInstance_";
            char bufferII[5];
            sprintf(bufferII, "%s", sigla);
            strcat (fil,bufferII);
            strcat (fil,".dat");

            FILE* app = fopen(fil, "a");
            printf("\n%f %.8f", timeP, prob );
            fprintf(app,"\n%f %.8f", timeP, prob );
            fclose(app);
        }
        if(strcmp("SSIW", name)==0) {
            sigla = "SSIW";
            char fil[256] = "chartProbabilitiesInstance_";
            char bufferII[5];
            sprintf(bufferII, "%s", sigla);
            strcat (fil,bufferII);
            strcat (fil,".dat");

            FILE* app = fopen(fil, "a");
            printf("\n%f %.8f", timeP, prob );
            fprintf(app,"\n%f %.8f", timeP, prob );
            fclose(app);
        }
        if(strcmp("SPE", name)==0) {
            sigla = "SPE";
            char fil[256] = "chartProbabilitiesInstance_";
            char bufferII[5];
            sprintf(bufferII, "%s", sigla);
            strcat (fil,bufferII);
            strcat (fil,".dat");

            FILE* app = fopen(fil, "a");
            printf("\n%f %.8f", timeP, prob );
            fprintf(app,"\n%f %.8f", timeP, prob );
            fclose(app);
        }
        if(strcmp("CSP", name)==0) {
            sigla = "CSP";
            char fil[256] = "chartProbabilitiesInstance_";
            char bufferII[5];
            sprintf(bufferII, "%s", sigla);
            strcat (fil,bufferII);
            strcat (fil,".dat");

            FILE* app = fopen(fil, "a");
            printf("\n%f %.8f", timeP, prob );
            fprintf(app,"\n%f %.8f", timeP, prob );
            fclose(app);
        }


    }

    fclose(ap);
}

void Test_writeResultVND(Test *test, char **argv, Cost initFO, int *assortment, Solution* sol, int fi)
{
    FILE* ap;

    if(fi)
        ap = fopen("resultsVND_1.txt", "a");
    else
        ap = fopen("resultsVND_0.txt", "a");

    fprintf(ap, "%s & %ld & %ld & %ld & %ld & ", argv[2], initFO, Sol_getTPD(sol), Sol_getTMS(sol), Sol_getCost(sol));

    for(int i=0; i<test->nNeigh; ++i)
        fprintf(ap, "%d & ", assortment[i]);

    for(int i=1; i<=test->nNeigh; ++i)
        fprintf(ap, "%d & ", Test_getImproveNeigh(test, i));

    for(int i=1; i<=test->nNeigh; ++i)
        fprintf(ap, "%d & ", Test_getVisitNeigh(test, i));

    for(int i=1; i<=test->nNeigh; ++i)
        fprintf(ap, "%ld & ", Test_getImproveFO(test, i));

    for(int i=1; i<=test->nNeigh; ++i)
        fprintf(ap, "%.4f & ", Test_getTimeNeigh(test, i));

    fprintf(ap, "%.4f ", Test_getTotalTime(test));

    fprintf(ap, "\n");

    fclose(ap);

}

void Test_writeResultVNS(Test *test, char **argv, Cost initFO, int *assortment, int first, int type, Solution* sol,  int lfa, int itRNA, int itLAHC, int nMoves, int nSizeSamplingShake, double perc, double percRS, int tm, int tj, int rm, int rj, int nSol)
{

    char fil[256] = "resultsVNS.txt";

    FILE* ap = fopen(fil, "a");

    fprintf(ap, "%s & %ld & %ld & %ld & %ld & ", argv[2], initFO, Sol_getTPD(sol), Sol_getTMS(sol), Sol_getCost(sol));

    for(int i=0; i<test->nNeigh; ++i)
        fprintf(ap, "%d & ", assortment[i]);

    for(int i=1; i<=test->nNeigh; ++i)
        fprintf(ap, "%d & ", Test_getImproveNeigh(test, i));

    for(int i=1; i<=test->nNeigh; ++i)
        fprintf(ap, "%d & ", Test_getVisitNeigh(test, i));

    for(int i=1; i<=test->nNeigh; ++i)
        fprintf(ap, "%ld & ", Test_getImproveFO(test, i));

    for(int i=1; i<=test->nNeigh; ++i)
        fprintf(ap, "%.4f & ", Test_getTimeNeigh(test, i));

    fprintf(ap, "%.4f & ", Test_getTotalTime(test));

    fprintf(ap, "%d & ", first);
    fprintf(ap, "%d & ", type);
    fprintf(ap, "%d & ", nMoves);
    fprintf(ap, "%d & ", nSizeSamplingShake);
    fprintf(ap, "%d & ", itRNA);
    fprintf(ap, "%d & ", itLAHC);
    fprintf(ap, "%d & ", lfa);
    fprintf(ap, "%f & ", perc);
    fprintf(ap, "%f & ", percRS);
    fprintf(ap, "%d & ", tm);
    fprintf(ap, "%d & ", tj);
    fprintf(ap, "%d & ", rm);
    fprintf(ap, "%d & ", rj);
    fprintf(ap, "%d ", nSol);

    fprintf(ap, "\n");

    fclose(ap);

}

void Test_writeResultLAHC_thread(Test *test, char **argv, Cost initFO, Neighborhood *neigh, Solution* sol, int lfa, int nSol, int seed, int thread)
{

    int i=0;

    char name[20] = {"resultsLAHC_tx.txt"};

    if(thread == 0)
        name[13] = '0';
    else if(thread == 1)
        name[13] = '1';
    else if(thread == 2)
        name[13] = '2';
    else if(thread == 3)
        name[13] = '3';

    //FILE* ap = fopen("resultsLAHC_t.txt", "a");
    FILE* ap = fopen(name, "a");

    fprintf(ap, "%s & %ld & %ld & %ld & %ld & ", argv[2], initFO, Sol_getTPD(sol), Sol_getTMS(sol), Sol_getCost(sol));

    for(i=0; i<test->nNeigh; ++i)
        fprintf(ap, "%.3Lf & ", Neighbor_getIdIntensity(neigh, i));

    for(i=1; i<=test->nNeigh; ++i)
        fprintf(ap, "%d & ", Test_getImproveNeigh(test, i));

    for(i=1; i<=test->nNeigh; ++i)
        fprintf(ap, "%d & ", Test_getVisitNeigh(test, i));

    for(i=1; i<=test->nNeigh; ++i)
        fprintf(ap, "%ld & ", Test_getImproveFO(test, i));

    for(i=1; i<=test->nNeigh; ++i)
        fprintf(ap, "%.4f & ", Test_getTimeNeigh(test, i));


    fprintf(ap, "%.4f & %d & ", Test_getTotalTime(test),lfa);


    for(i=0; i<4; ++i)
        fprintf(ap, "%d & ", Neighbor_getMinK(neigh, i));

    for(i=0; i<4; ++i)
        fprintf(ap, "%d & ", Neighbor_getMaxK(neigh, i));


    fprintf(ap, "%d & %d & %d \n", nSol, seed, thread);



    fclose(ap);


}

void Test_writeResultLAHC(Test *test, char **argv, Cost initFO, Neighborhood *neigh, Solution* sol, int nSol, int lfa, int seed, int nCostList, int nDiversification, int nStayDiversification, int nWOImprove, double perc, int sw, double lr, int it, float psw, int tm, int tj, int rm, int rj )
{
    int i=0;

    char fil[256] = "resultsLAHC_";
    char sideway[256] ="";
    char lalr[256] ="";
    char itUp[256] ="";
    char pfsw[256] ="";
    char slfa[256] ="";
    char div[256]  ="";
    char stay[256] ="";
    char perd[256] ="";
    char wimp[256] ="";
    char ctm[256] ="";
    char ctj[256] ="";
    char crm[256] ="";
    char crj[256] ="";

    sprintf(sideway, "%d", sw);
    strcat (fil,sideway);

    strcat (fil,"_");
    sprintf(perd, "%.2f", perc);
    strcat (fil,perd);
    strcat (fil,"_");
    sprintf(lalr, "%.8f", lr);
    strcat (fil,lalr);
    strcat (fil,"_");
    sprintf(pfsw, "%.8f", psw);
    strcat (fil,pfsw);
    strcat (fil,"_");
    sprintf(itUp, "%d", it);
    strcat (fil,itUp);

    strcat (fil,".txt");

    FILE* ap = fopen(fil, "a");

    fprintf(ap, "%s & %ld & %ld & %ld & %ld & ", argv[2], initFO, Sol_getTPD(sol), Sol_getTMS(sol), Sol_getCost(sol));

    for(i=0; i<test->nNeigh; ++i)
        fprintf(ap, "%.3Lf & ", Neighbor_getIdIntensity(neigh, i));

    for(i=1; i<=test->nNeigh; ++i)
        fprintf(ap, "%d & ", Test_getImproveNeigh(test, i));

    for(i=1; i<=test->nNeigh; ++i)
        fprintf(ap, "%d & ", Test_getVisitNeigh(test, i));

    for(i=1; i<=test->nNeigh; ++i)
        fprintf(ap, "%ld & ", Test_getImproveFO(test, i));

    for(i=1; i<=test->nNeigh; ++i)
        fprintf(ap, "%.4f & ", Test_getTimeNeigh(test, i));


    fprintf(ap, "%.4f & %d & ", Test_getTotalTime(test),lfa);


    for(i=0; i<4; ++i)
        fprintf(ap, "%d & ", Neighbor_getMinK(neigh, i));

    for(i=0; i<4; ++i)
        fprintf(ap, "%d & ", Neighbor_getMaxK(neigh, i));

    fprintf(ap, "%d & ", nDiversification);
    fprintf(ap, "%d & ", nStayDiversification);
    fprintf(ap, "%d & ", nWOImprove);
    fprintf(ap, "%d & ", nCostList);
    fprintf(ap, "%d & ", lfa);
    fprintf(ap, "%f & ", perc);
    //    fprintf(ap, "%f & ", percRS);
    fprintf(ap, "%d & ", tm);
    fprintf(ap, "%d & ", tj);
    fprintf(ap, "%d & ", rm);
    fprintf(ap, "%d & ", rj);
    fprintf(ap, "%d \n", nSol);

    fclose(ap);

}

void Test_writeResultSA(Test *test, char **argv, Cost initFO, int *assortment, Solution* sol)
{

    FILE* ap = fopen("resultsSA.txt", "a");

    fprintf(ap, "%s & %ld & %ld & %ld & %ld & ", argv[2], initFO, Sol_getTPD(sol), Sol_getTMS(sol), Sol_getCost(sol));

    for(int i=0; i<test->nNeigh; ++i)
        fprintf(ap, "%d & ", assortment[i]);

    for(int i=1; i<=test->nNeigh; ++i)
        fprintf(ap, "%d & ", Test_getImproveNeigh(test, i));

    for(int i=1; i<=test->nNeigh; ++i)
        fprintf(ap, "%d & ", Test_getVisitNeigh(test, i));

    for(int i=1; i<=test->nNeigh; ++i)
        fprintf(ap, "%ld & ", Test_getImproveFO(test, i));

    for(int i=1; i<=test->nNeigh; ++i)
        fprintf(ap, "%.4f & ", Test_getTimeNeigh(test, i));

    fprintf(ap, "%d & ", Test_getT(test));

    fprintf(ap, "%d & ", Test_getSAmax(test));

    fprintf(ap, "%.2f & ", Test_getAlpha(test));

    fprintf(ap, "%.4f ", Test_getTotalTime(test));

    //if( lfa>0 )
    //  fprintf(ap, "& %d ", lfa);

    fprintf(ap, "\n");

    fclose(ap);

}

void Test_setImproveNeigh(Test *test, int idxVet, int value)
{
    assert( test != NULL );
    assert( idxVet > 0 );
    assert( value >= 0 );

    test->improveNeigh[idxVet-1] = value;

}

void Test_incrementImproveNeigh(Test *test, int idxVet, int value)
{
    assert( test != NULL );
    assert( idxVet > 0 );
    assert( value >= 0 );

    test->improveNeigh[idxVet-1] = test->improveNeigh[idxVet-1] + value;

}

void Test_setVisitNeigh(Test *test, int idxVet, int value)
{
    assert( test != NULL );
    assert( idxVet > 0 );
    assert( value >= 0 );

    test->visitNeigh[idxVet-1] = value;

}

void Test_incrementVisitNeigh(Test *test, int idxVet, int value)
{
    assert( test != NULL );
    assert( idxVet > 0 );
    assert( value >= 0 );

    test->visitNeigh[idxVet-1] = test->visitNeigh[idxVet-1] + value;

}

void Test_setTimeNeigh(Test *test, int idxVet, double value)
{
    assert( test != NULL );
    assert( idxVet > 0 );
    assert( value >= 0 );

    test->timeNeigh[idxVet-1] = value;

}

void Test_incrementTimeNeigh(Test *test, int idxVet, double value)
{
    assert( test != NULL );
    assert( idxVet > 0 );
    assert( value >= 0 );

    test->timeNeigh[idxVet-1] = test->timeNeigh[idxVet-1] + value;
}

void Test_setTotalTime(Test *test, double value)
{
    assert( test != NULL );
    assert( value >= 0 );

    test->totalTime = value;
}

void Test_setCurrentFO(Test *test, Cost value)
{
    assert( test != NULL );
    assert( value >= 0 );

    test->currentFO = value;
}

void Test_setImproveFO(Test *test, int idxVet, Cost value)
{
    assert( test != NULL );
    assert( idxVet > 0 );
    assert( value >= 0 );

    test->improveFO[idxVet-1] = value;
}

void Test_setCurrentTime(Test *test, double _time)
{
    assert( test != NULL );
    assert( _time >= 0 );

    test->currentTime = _time;
}

void Test_setCurrentNeigh(Test *test, int nNeigh)
{
    assert( test != NULL );
    assert( nNeigh > 0 );

    test->currentNeigh = nNeigh;
}

void Test_setT(Test *test, int t)
{
    assert( test != NULL );
    assert( t > 0 );

    test->T = t;
}

void Test_setSAmax(Test *test, int samax)
{
    assert( test != NULL );
    assert( samax > 0 );

    test->SAmax = samax;
}

void Test_setAlpha(Test *test, double alpha)
{
    assert( test != NULL );
    assert( alpha > 0 );

    test->alpha = alpha;
}

void Test_incrementImproveFO(Test *test, int idxVet, Cost value)
{
    assert( test != NULL );
    assert( idxVet > 0 );
    assert( value >= 0 );

    test->improveFO[idxVet-1] = test->improveFO[idxVet-1] + value;
}

double Test_getAlpha(Test *test)
{
    assert( test != NULL );

    return test->alpha;
}

int Test_getT(Test *test)
{
    assert( test != NULL );

    return test->T;
}

int Test_getSAmax(Test *test)
{
    assert( test != NULL );

    return test->SAmax;
}

Cost Test_getImproveFO(Test *test, int idxVet)
{
    assert( test != NULL );
    assert( idxVet > 0 );

    return test->improveFO[idxVet-1];
}

int Test_getImproveNeigh(Test *test, int idxVet)
{
    assert( test != NULL );
    assert( idxVet > 0 );

    return test->improveNeigh[idxVet-1];

}

int Test_getVisitNeigh(Test *test, int idxVet)
{
    assert( test != NULL );
    assert( idxVet > 0 );

    return test->visitNeigh[idxVet-1];
}

double Test_getTimeNeigh(Test *test, int idxVet)
{
    assert( test != NULL );
    assert( idxVet > 0 );

    return test->timeNeigh[idxVet-1];
}

double Test_getTotalTime(Test *test)
{
    assert( test != NULL );

    return test->totalTime;
}

Cost Test_getCurrentFO(Test * test)
{
    assert( test != NULL );

    return test->currentFO;
}

double Test_getCurrentTime(Test *test)
{
    assert( test != NULL );

    return test->currentTime;
}

int Test_getCurrentNeigh(Test *test)
{
    assert( test != NULL );

    return test->currentNeigh;
}

void Test_callTest(Test *test, int idNeigh, double _timeNeigh, Cost bestFO)
{

    if((Test_getCurrentFO(test)-bestFO) > 0)
        Test_incrementImproveNeigh(test, idNeigh, 1);

    Test_incrementVisitNeigh(test, idNeigh, 1);
    Test_incrementTimeNeigh(test, idNeigh, _timeNeigh);

    if((Test_getCurrentFO(test)-bestFO)>=0)
        Test_incrementImproveFO(test, idNeigh, (Test_getCurrentFO(test)-bestFO));

    //printf("\nTEST: improveFO[%d]: %ld | improveNeigh[]: %d | visitNeigh[]: %d | timeNeigh[]: %.6f |", idNeigh, test->improveFO[idNeigh-1], test->improveNeigh[idNeigh-1], test->visitNeigh[idNeigh-1], test->timeNeigh[idNeigh-1]);
    //getchar();

}

void Test_rebuild(Instance *inst, Neighborhood *neigh, Solution *sol1, Solution *sol2,  Solution *bestSol1,  Solution *bestSol2)
{

    Solution *solRebuildNor = Sol_create(inst);
    Solution *solRebuildOpt = Sol_create(inst);
    Sol_cpy(solRebuildNor, sol1);
    Sol_cpy(solRebuildOpt, sol2);

    const ModeSet *ms1 = Sol_getModeSet(solRebuildNor);
    const ModeSet *ms2 = Sol_getModeSet(solRebuildOpt);

    //   int  improve = 0;
    //   improve = Neighbor_search_changeTwoMode( neigh, sol, bestSol, solRebuildNor, 8, 0, 10);
    //    improve = Neighbor_search_changeTwoMode( neigh, sol, bestSol, solRebuildOpt, 8, 0, 10);

    int nJobs = Inst_nJobs( Sol_inst(sol1) );
    int teste =0;

    int valid1 = 0, valid2 = 0;
    for(int j1 = 0; j1 < nJobs; j1++) {
        const Job* job1 = Inst_job(inst,j1);
        for(int m1 = 0; m1 < Job_nModes(job1); m1++) {
            const Mode * mode1 = Job_mode(job1,m1);
            if(!Mode_isFeasible(inst,mode1)) continue;
            for(int j2 = j1+1; j2 < nJobs; j2++) {
                const Job* job2 = Inst_job(inst,j2);
                for(int m2 = 0; m2 < Job_nModes(job2); m2++) {
                    const Mode * mode2 = Job_mode(job2,m2);
                    if(!Mode_isFeasible(inst,mode2)) continue;


                    valid1 =  Neighbor_changeTwoMode(neigh, solRebuildNor, j1, j2, m1, m2);
                    valid2 =  Neighbor_changeTwoMode(neigh, solRebuildOpt, j1, j2, m1, m2);
                    if(!valid1) {
                        Sol_cpy(solRebuildNor,sol1);
                        valid1 = 0;
                        teste = 1;
                        printf("invalid 1");
                        //continue;
                    }

                    if(!valid2) {
                        Sol_cpy(solRebuildOpt,sol2);
                        valid2 = 0;
                        teste = 1;
                        printf("invalid 2");
                        //continue;
                    }

                    if(teste) {
                        //getchar();
                        teste = 0;
                        continue;
                    }

                    Sol_rebuild(solRebuildNor);
                    Sol_rebuild_opt(solRebuildOpt, sol2);
                    if(Sol_getCost(solRebuildNor) != Sol_getCost(solRebuildOpt) ) {
                        printf("\n-NOR: %ld OPT: %ld\n ", Sol_getCost(solRebuildNor), Sol_getCost(solRebuildOpt));
                        printf("\nNOR Modes\n ");
                        for(int n = 0; n < nJobs; n++)
                            printf("%d ", Modes_job(ms1,n) );
                        printf("\nOPT Modes\n ");

                        for(int n = 0; n < nJobs; n++)
                            printf("%d ", Modes_job(ms2,n)  );

                        printf("\nNOR seq \n ");
                        for(int n = 0; n < nJobs; n++)
                            printf("%d ", Sol_getSequence(solRebuildNor,n) );

                        printf("\nOPT seq\n ");
                        for(int n = 0; n < nJobs; n++)
                            printf("%d ", Sol_getSequence(solRebuildOpt,n) );
                        getchar();
                    } else
                        printf("\nNOR: %ld OPT: %ld\n ", Sol_getCost(solRebuildNor), Sol_getCost(solRebuildOpt));


                    if(Sol_getCost(solRebuildNor) < Sol_getCost(sol1) ) {
                        Sol_cpy(sol1,solRebuildNor);
                        if(Sol_getCost(solRebuildNor) < Sol_getCost(bestSol1)) {
                            //  printf("\nImprovement changeTwoMode %ld to %ld \n ", Sol_getCost(bestSol), Sol_getCost(sol) );
                            printf("\nimp nor %ld\n ", Sol_getCost(solRebuildNor));
                            Sol_cpy( bestSol1, sol1 );
                        }
                        printf("\n melhorou atual NOR de %ld para %ld \n ", Sol_getCost(sol1), Sol_getCost(solRebuildNor));
                        //  improve = 1;
                        // if(firstImp && kN) return improve;
                    } else
                        Sol_cpy(solRebuildNor,sol1);


                    if(Sol_getCost(solRebuildOpt) < Sol_getCost(sol2) ) {
                        Sol_cpy(sol2,solRebuildOpt);
                        if(Sol_getCost(solRebuildOpt) < Sol_getCost(bestSol2)) {
                            //  printf("\nImprovement changeTwoMode %ld to %ld \n ", Sol_getCost(bestSol), Sol_getCost(sol) );
                            printf("\n imp opt %ld\n ", Sol_getCost(solRebuildOpt));
                            Sol_cpy( bestSol2, sol2 );
                        }
                        printf("\n melhorou atual OPT de %ld para %ld \n ", Sol_getCost(sol2), Sol_getCost(solRebuildOpt));
                        //improve = 1;
                        // if(firstImp && kN) return improve;
                    } else
                        Sol_cpy(solRebuildOpt,sol2);

                    // _time = ( (double)( clock()-tStart )/CLOCKS_PER_SEC );
                    // if( _time > timeRem)  return improve;

                    //                  Sol_cpy(solRebuildNor,sol1);
                    //                Sol_cpy(solRebuildOpt,sol2);

                }
            }
        }
    }

}


void Test_calcIntensityNeigh(Test * test, Neighborhood *neigh, int s,  char ** argv)
{

    long double _time = 0;
    long double _currentTime = 0;
    clock_t tStart = clock();

    int valid =0 ;

    int n = 0, cont = 0;
    int pos;
    while(n < test->nNeigh) {
        while(cont < test->it ) {

            _time =  (double)( clock()-tStart ); ///CLOCKS_PER_SEC );
            Solution *current = Sol_create(test->inst);
            if(s == 0)
                Sol_read(current,  argv[6] );
            if(s == 1 )
                Sol_read(current,  argv[7] );

            Solution *sol = Sol_create(test->inst);
            Sol_cpy(sol,current);

            Solution *bestSol = Sol_create(test->inst);
            Sol_cpy(bestSol,current);


            Cost oldFO = Sol_getCost(current);

            unsigned int oldSolutionHash = Sol_getSolutionHash(current);

            valid = Neighbor_callStocChosen(neigh,sol,bestSol,current,n+1, test);
            if(valid)
                Sol_rebuild_opt(current, sol);

            pos = Neighbor_getIdxAssortment( neigh, n)-1;

            _currentTime = (long double)( clock()-tStart );///CLOCKS_PER_SEC );
            if(Sol_getCost(current) < oldFO ) {
                test->I[pos][s] += 1.0;//(double) MAX(0,(double) oldFO-(double) newFo);
                test->TI[pos][s] += (long double) (_currentTime - _time);
            }
            if(Sol_getCost(current) == oldFO ) {
                if(Sol_getSolutionHash(current) != oldSolutionHash) {
                    test->EQ[pos][s] += (long double) 0.1;
                    test->TE[pos][s] += (long double) (_currentTime - _time);
                } else
                    test->TIV[pos][s] += (long double) (_currentTime - _time);
            }
            if(Sol_getCost(current) > oldFO )
                test->TIV[pos][s] += (long double) (_currentTime - _time);
            // } else {
            //     _currentTime = (long double)( clock()-tStart );///CLOCKS_PER_SEC );
            //     test->TIV[pos][s] += (long double) (_currentTime - _time);
            // }

            Sol_free(&current);
            Sol_free(&sol);
            Sol_free(&bestSol);
            Neighbor_setNullLastJobModify(neigh);
            cont++;
        }

        pos = Neighbor_getIdxAssortment( neigh, n)-1;
        double value = (test->TI[pos][s]+ test->TE[pos][s] + test->TIV[pos][s]);
        double valueTime =0;
        double valueI = 0;
        double valueE = 0;

        int log = atoi(argv[9]);

        if(value != 0) {
            if(log)
                valueTime = log2(value);
            else
                valueTime = value;

            valueI = (long double) test->I[pos][s]/(long double) (valueTime);
            valueE = (long double) ( test->I[pos][s]+ test->EQ[pos][s]) / (long double)(valueTime);
        }

        test->FI[pos][s] = valueI;
        //   printf("FI[%d][%d] => %Lf\n",pos,s, test->FI[pos][s]);
        test->FE[pos][s] = valueE;
        //  printf("FE[%d][%d] => %Lf\n",pos,s, test->FE[pos][s]);

        cont = 0;
        n++;
    }


    //Túlio normalizando por cima
    for(int i=0; i<test->nNeigh; ++i) {
        test->maxFI[s] = MAX(test->maxFI[s], test->FI[i][s]);
        test->minFI[s] = MIN(test->minFI[s], test->FI[i][s]);
        test->maxFE[s] = MAX(test->maxFE[s], test->FE[i][s]);
        test->minFE[s] = MIN(test->minFE[s], test->FE[i][s]);

    }
    test->intervalI[s] = test->maxFI[s] - test->minFI[s];
    test->intervalE[s] = test->maxFE[s] - test->minFE[s];

    for(int i=0; i<test->nNeigh; ++i) {
        if(test->intervalI[s] == 0) test->normFI[i][s] = 1;
        else test->normFI[i][s] = (long double)(test->FI[i][s])/(long double)test->maxFI[s];

        if(test->intervalE[s] == 0) test->normFE[i][s] = 1;
        else test->normFE[i][s] = (long double)(test->FE[i][s])/(long double)test->maxFE[s];
        //   printf("normFI[%d][%d] => %Lf normFE[%d][%d] => %Lf \n",i,s, test->normFI[i][s], i,s, test->normFE[i][s]);
    }


    /*Janniele normalizando por baixo
        for(int i=0; i<test->nNeigh; ++i) {
            if(test->intervalI[s] == 0) test->normFI[i][s] = 0;
            else test->normFI[i][s] = (double)(test->FI[i][s]- test->minFI[s])/(double)test->intervalI[s];

            if(test->intervalE[s] == 0) test->normFE[i][s] = 0;
            else test->normFE[i][s] = (double)(test->FE[i][s]- test->minFE[s])/(double)test->intervalE[s];

            //printf("normFI[%d][%d] => %f normFE[%d][%d] => %f \n",i,s, test->normFI[i][s], i,s, test->normFE[i][s]);
        }
    */

}


void Test_neigh(Instance *inst, Neighborhood *neigh, char ** argv, int argc, Test *test)
{

    test->it = atoi(argv[8]);

    int s = 0;

    while(s < test->nStages) {

        Test_calcIntensityNeigh(test, neigh, s, argv);

        s++;

    }


    //Test_readAnalisysNeigh_offline(test, argv);


}

void Test_writeAnalisysNeigh(Test *test,  char **argv)
{
    int i=0;

    char fi[12] = "intensities_";

    strcat (fi, argv[4]);
    strcat (fi,".txt ");

    FILE* ap = fopen(fi, "a");

    for(int s = 0 ; s < test->nStages ; s++) {
        for(i=0; i<test->nNeigh; ++i)
            fprintf(ap, " %s %d %d %d %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf\n",  argv[2], test->it, s, i, test->I[i][s], test->TI[i][s], test->EQ[i][s], test->TE[i][s], test->TIV[i][s], test->FI[i][s], test->normFI[i][s], test->FE[i][s], test->normFE[i][s]);
    }

    fclose(ap);

}

/*
void Test_readAnalisysNeigh_offline(Test *test,  char **argv)
{

    char inst[20], virgula, ecomercial;
    int it=0, s, n;
    float I, TI, EQ, TE, TIV, FI, NFI, FE, NFE;
    //double I, TI, EQ, TE, TIV, FI, NFI, FE, NFE;
    double sumFI[2], sumFE[2];
    sumFI[0] = 0;
    sumFI[1] = 0;
    sumFE[0] = 0;
    sumFE[1] = 0;


    for(int i = 1 ; i <=10 ; i++) {

        char fi[256] = "intensities_";

        strcat (fi, argv[4]);
        strcat (fi,".txt ");

        FILE* ap = fopen(fi, "r");
        if(ap == NULL)
            printf("Erro, nao foi possivel abrir o arquivo\n");
        else {
            while( (fscanf(ap,"%s %d %d %d %f %f %f %f %f %f %f %f %f\n", inst, &it, &s, &n, &I, &TI, &EQ, &TE, &TIV, &FI, &NFI, &FE, &NFE))!=EOF ) {
                //printf("%s %d %d %d %f %f %f %f %f %f %f %f %f \n", inst, it, s, n,  I,  TI, EQ, TE, TIV, FI, NFI, FE, NFE);
                test->I[n][s] += (double) I;
                test->EQ[n][s] += (double) EQ;
                test->TI[n][s] += (double) TI;
                test->TE[n][s] += (double) TE;
                test->TIV[n][s] += (double) TIV;
                test->FI[n][s] += (double) FI;
                test->FE[n][s] += (double) FE;
                test->normFI[n][s] += NFI;
                test->normFE[n][s] += NFE;
            }

        }
        fclose(ap);
    }





    for(int s=0; s <2 ; s++) {

        for(int i=0; i<test->nNeigh; ++i) {
            test->maxFI[s] = MAX(test->maxFI[s], test->normFI[i][s]);
            test->minFI[s] = MIN(test->minFI[s], test->normFI[i][s]);
            test->maxFE[s] = MAX(test->maxFE[s], test->normFE[i][s]);
            test->minFE[s] = MIN(test->minFE[s], test->normFE[i][s]);

        }
        test->intervalI[s] = test->maxFI[s] - test->minFI[s];
        test->intervalE[s] = test->maxFE[s] - test->minFE[s];

        for(int i=0; i<test->nNeigh; ++i) {
            if(test->intervalI[s] == 0) test->normFI[i][s] = 1;
            else test->normFI[i][s] = (double)(test->normFI[i][s])/(double)test->maxFI[s];

            if(test->intervalE[s] == 0) test->normFE[i][s] = 1;
            else test->normFE[i][s] = (double)(test->normFE[i][s])/(double)test->maxFE[s];


            sumFI[s] += test->normFI[n][s];
            sumFE[s] += test->normFE[n][s];
        }

    }


    FILE* app = fopen("ResultIntensities.txt", "a");
    char* sigla;
    char* cor;



    for(int ii=0; ii <2 ; ii++) {
        for(int i=0; i<test->nNeigh; i++) {
            switch(i) {
                case seqInvert-1:
                    cor="lightgray";
                    sigla = "ISJ";
                    break;
                case seqShiftJob-1:
                    cor="teal";
                    sigla = "OJ";
                    break;
                case seqSwapJob-1:
                    sigla = "STJ";
                    cor="magenta";
                    break;
                case seqShiftProj-1:
                    sigla = "OP";
                    cor="purple";
                    break;
                case seqSwapProj-1:
                    sigla = "SCTP";
                    cor= "cyan";
                    break;
                case seqCompactProj-1:
                    sigla = "CPP";
                    cor="orange";
                    break;
                case changeOneMode-1:
                    sigla = "COM";
                    cor="green";
                    break;
                case changeTwoMode-1:
                    sigla = "CTM";
                    cor="brown";
                    break;
                case changeThreeMode-1:
                    sigla = "CThM";
                    cor= "black";
                    break;
                case changeFourMode-1:
                    sigla = "CFM";
                    cor="pink";
                    break;
                case seqSwapJobFILS-1:
                    sigla = "SSJW";
                    cor="gray";
                    break;
                case seqInsertJobFILS-1:
                    sigla = "SSIW";
                    cor = "yellow";
                    break;
                case seqCompOnExtrem-1:
                    sigla = "SPE";
                    cor= "red";
                    break;
                case seqMoveProj-1:
                    cor= "blue";
                    sigla = "CSP";
                    break;
            }

            fprintf(app, "\n%d %s %Lf %Lf", ii,  sigla, test->normFI[i][ii],  test->normFE[i][ii] );
        }
        fprintf(app, "\n");
    }


    fclose(app);


    for(int ii=0; ii <2 ; ii++) {

        char fil[12] = "pieChart_I_";

        char bufferII[3];
        sprintf(bufferII, "%d", ii);
        //char bufferI[3];
        //sprintf(bufferI, "%d", i);

        strcat (fil,bufferII);
        //strcat (fil,"_");
        //strcat (fil, bufferI);
        strcat (fil,".tex ");

        FILE* appp = fopen(fil, "a");


        fprintf(appp,"\n\n \\def\\angle\{0\} \n \\def\\radius\{3\} \n \\def\\cyclelist\{\{\"lightgray\",\"teal\",\"magenta\",\"purple\",\"cyan\",\"orange\",\"green\",\"brown\",\"black\",\"pink\",\"gray\",\"yellow\",\"red\",\"blue\"\}\} \n \\newcount\\cyclecount \\cyclecount=-1 \n \\newcount\\ind \\ind=-1 \n \\begin\{tikzpicture\}[nodes = \{font=\\tiny\}] \n  \\foreach \\percent\/\\name in \{", cor);


        char filE[12] = "pieChart_E_";

        char bufferIIE[3];
        sprintf(bufferIIE, "%d", ii);
        //char bufferI[3];
        //sprintf(bufferI, "%d", i);

        strcat (filE,bufferIIE);
        //strcat (fil,"_");
        //strcat (fil, bufferI);
        strcat (filE,".tex ");

        FILE* apppE = fopen(filE, "a");


        fprintf(apppE,"\n\n \\def\\angle\{0\} \n \\def\\radius\{3\} \n \\def\\cyclelist\{\{\"lightgray\",\"teal\",\"magenta\",\"purple\",\"cyan\",\"orange\",\"green\",\"brown\",\"black\",\"pink\",\"gray\",\"yellow\",\"red\",\"blue\"\}\} \n \\newcount\\cyclecount \\cyclecount=-1 \n \\newcount\\ind \\ind=-1 \n \\begin\{tikzpicture\}[nodes = \{font=\\tiny\}] \n  \\foreach \\percent\/\\name in \{", cor);

        for(int i=0; i<test->nNeigh; i++) {
            switch(i) {
                case seqInvert-1:
                    cor="lightgray";
                    sigla = "ISJ";
                    break;
                case seqShiftJob-1:
                    cor="teal";
                    sigla = "OJ";
                    break;
                case seqSwapJob-1:
                    sigla = "STJ";
                    cor="magenta";
                    break;
                case seqShiftProj-1:
                    sigla = "OP";
                    cor="purple";
                    break;
                case seqSwapProj-1:
                    sigla = "SCTP";
                    cor= "cyan";
                    break;
                case seqCompactProj-1:
                    sigla = "CPP";
                    cor="orange";
                    break;
                case changeOneMode-1:
                    sigla = "COM";
                    cor="green";
                    break;
                case changeTwoMode-1:
                    sigla = "CTM";
                    cor="brown";
                    break;
                case changeThreeMode-1:
                    sigla = "CThM";
                    cor= "black";
                    break;
                case changeFourMode-1:
                    sigla = "CFM";
                    cor="pink";
                    break;
                case seqSwapJobFILS-1:
                    sigla = "SSJW";
                    cor="gray";
                    break;
                case seqInsertJobFILS-1:
                    sigla = "SSIW";
                    cor = "yellow";
                    break;
                case seqCompOnExtrem-1:
                    sigla = "SPE";
                    cor= "red";
                    break;
                case seqMoveProj-1:
                    cor= "blue";
                    sigla = "CSP";
                    break;
            }


            //fprintf(appp,"\n%.3f/%s,",test->normFI[i][ii]*100,sigla);
            if(i==test->nNeigh)
                fprintf(appp,"\n%.3Lf/%s",(long double) (test->normFI[i][ii]/sumFI[ii])*100, sigla); // *100,sigla);
            else
                fprintf(appp,"\n%.3Lf/%s,",(long double) (test->normFI[i][ii]/sumFI[ii])*100, sigla); // *100,sigla);
            // printf("\nI %.3f/%s,",  test->normFI[i][ii]/ sumFi[ii], sigla);
            if(i==test->nNeigh)
                fprintf(apppE,"\n%.3Lf/%s",(long double) (test->normFE[i][ii]/sumFE[ii])*100, sigla); // *100,sigla);
            else fprintf(apppE,"\n%.3Lf/%s,",(long double) (test->normFE[i][ii]/sumFE[ii])*100, sigla); // *100,sigla);
            //printf("\nE %.3f/%s,",  test->normFE[i][ii]/ sumFE[ii], sigla);


            //            fprintf(appp,"\n%.3f/others",100-(test->normFI[i][ii]*100));

        }
        fprintf(appp, "\n \} \{ \n \\ifx\\percent\\empty\\else              \n \\global\\advance\\cyclecount by 1    \n \\global\\advance\\ind by 1           \n \\ifnum14<\\cyclecount                \n   \\global\\cyclecount=0             \n   \\global\\ind=0                     \n \\fi       \n \\pgfmathparse\{\\cyclelist[\\the\\ind]\}\n \\edef\\color\{\\pgfmathresult\}         \n \\draw[fill=\{\\color!50\},draw=\{\\color\}] (0,0) -- (\\angle:\\radius)\n  arc (\\angle:\\angle+\\percent*3.6:\\radius) -- cycle;\n \\node at (\\angle+0.5*\\percent*3.6:0.7*\\radius) \{\\percent\\,\\%\};\n \\node[pin=\\angle+0.5*\\percent*3.6:\\name]\n  at (\\angle+0.5*\\percent*3.6:\\radius) \{\};\n \\pgfmathparse\{\\angle+\\percent*3.6\} \n  \\xdef\\angle\{\\pgfmathresult\}        \n \\fi \n \}; \n \\end\{tikzpicture\}");
        fprintf(apppE, "\n \} \{ \n \\ifx\\percent\\empty\\else              \n \\global\\advance\\cyclecount by 1    \n \\global\\advance\\ind by 1           \n \\ifnum14<\\cyclecount                \n   \\global\\cyclecount=0             \n   \\global\\ind=0                     \n \\fi       \n \\pgfmathparse\{\\cyclelist[\\the\\ind]\}\n \\edef\\color\{\\pgfmathresult\}         \n \\draw[fill=\{\\color!50\},draw=\{\\color\}] (0,0) -- (\\angle:\\radius)\n  arc (\\angle:\\angle+\\percent*3.6:\\radius) -- cycle;\n \\node at (\\angle+0.5*\\percent*3.6:0.7*\\radius) \{\\percent\\,\\%\};\n \\node[pin=\\angle+0.5*\\percent*3.6:\\name]\n  at (\\angle+0.5*\\percent*3.6:\\radius) \{\};\n \\pgfmathparse\{\\angle+\\percent*3.6\} \n  \\xdef\\angle\{\\pgfmathresult\}        \n \\fi \n \}; \n \\end\{tikzpicture\}");
        printf( "\n\n\n");
        fclose(appp);
        fclose(apppE);
    }



}

*/
int Test_fileExists(char fileName[])
{

    FILE *fp;

    //char mainPath[500] = "";
    //strcat(mainPath, fileName);
    //printf("\nmainPath: %s", mainPath);

    fp=fopen(fileName,"r");

    if(fp) {
        fclose(fp);
        return 1;
    } else
        return 0;

    printf("\nErro in fileExists!");
    return -1;

}
