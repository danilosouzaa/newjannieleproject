#include "cut_gpu.h"

Cut_gpu *AllocationStructCut(int cont, int nConstrains, int nVariables)
{
    size_t size_cut = sizeof(Cut_gpu) +
                      sizeof(TCoefficients)*(cont) +
                      sizeof(TElements)*(cont) +
                      sizeof(TElementsConstraints)*(nConstrains+1) +
                      sizeof(TRightSide)*(nConstrains) +
                      sizeof(TXAsterisc)*(nVariables)+
                      sizeof(TTypeConstraints)*(nConstrains);


    Cut_gpu *cut = (Cut_gpu*)malloc(size_cut);
    assert(cut!=NULL);
    memset(cut,0,size_cut);
    cut->Coefficients = (TCoefficients*)(cut+1);
    cut->Elements = (TElements*)(cut->Coefficients + cont);
    cut->ElementsConstraints = (TElementsConstraints*)(cut->Elements + cont);
    cut->rightSide = (TRightSide*)(cut->ElementsConstraints + (nConstrains+1));
    cut->xAsterisc = (TXAsterisc*)(cut->rightSide + (nConstrains));
    cut->typeConstraints = (TTypeConstraints*)(cut->xAsterisc + (nVariables));
    cut->numberVariables = nVariables;
    cut->numberConstrains = nConstrains;
    cut->cont = cont;
    return cut;
}


Cut_gpu_aux *AllocationStructCutAux(int nConstrains, int nVariables, int nCont)
{
    size_t size_cut_aux = sizeof(Cut_gpu_aux) +
                          sizeof(TInterval)*(nConstrains)+
                          sizeof(TInterval)*(nConstrains)+
                          sizeof(TNames)*(nCont)+
                          sizeof(TNames)*(nConstrains);


    Cut_gpu_aux *cut_aux = (Cut_gpu_aux*)malloc(size_cut_aux);
    assert(cut_aux!=NULL);
    memset(cut_aux,0,size_cut_aux);
    cut_aux->intervalMin = (TInterval*)(cut_aux + 1);
    cut_aux->intervalMax = (TInterval*)(cut_aux->intervalMin + (nConstrains));
    cut_aux->nameElements = (TNames*)(cut_aux->intervalMax + (nConstrains));
    cut_aux->nameConstraints = (TNames*)(cut_aux->nameElements+(nCont));
    cut_aux->numberVariables = nVariables;
    cut_aux->numberConstrains = nConstrains;
    return cut_aux;

}

listNeigh *AllocationListNeigh(int nConstrains, int nList)
{
    size_t size_list = sizeof(listNeigh) +
                       sizeof(TList)*(nList)+
                       sizeof(TPosList)*(nConstrains+1);
    listNeigh *list_t = (listNeigh*)malloc(size_list);
    assert(list_t!=NULL);
    memset(list_t, 0,size_list);
    list_t->list_n = (TList*)(list_t + 1);
    list_t->pos = (TPosList*)(list_t->list_n + nList);
    list_t->nList = nList;
    list_t->nPos = nConstrains + 1;
    return list_t;
}

int returnIndVector(TNames *v,char *nome, int sz)
{
    int i;
    for(i=0; i<sz; i++)
    {
        if(strcmp(v[i].name,nome)==0)
            return i;
    }
    return -1;
}

Cut_gpu* fillStruct(CutCG *ccg_r2,int precision, int numberVariables, int numberConstrains, int cont)
{
   Cut_gpu *ccg =  AllocationStructCut(cont, numberConstrains, numberVariables);
   int i,pos = 0, el, elem;
   ccg->ElementsConstraints[0] = 0;
   for(i = 0; i < VStr_size(ccg_r2->rname); i++)
   {
       ccg->rightSide[i] = VDbl_get(ccg_r2->rowrhs,i);
       ccg->typeConstraints[i] = VInt_get(ccg_r2->rowtype,i);
       if(i<numberConstrains-1)
           ccg->ElementsConstraints[i+1] = ccg->ElementsConstraints[i] + VDbl_size(ccg_r2->rowCoef[i]);
       else
           ccg->ElementsConstraints[i+1] = ccg->cont;
       for( el = 0 ; el < VDbl_size(ccg_r2->rowCoef[i]) ; el++)
       {
           elem = VInt_get(ccg_r2->rowElem[i],el);
           ccg->Coefficients[pos] = VDbl_get(ccg_r2->rowCoef[i],el);
           ccg->Elements[pos] = elem;
           pos++;
       }

   }
   for(i=0; i<ccg->numberVariables; i++)
   {
       ccg->xAsterisc[i] = precision * VDbl_get(ccg_r2->xfElemPP,i);
   }
   return ccg;
}
void setParameters_ccg(parameters_ccg *parCCG, int mode)
{
    parCCG->precision=1000;
    if(mode==0)
    {
        parCCG->numberMaxConst = 8; //m'
        parCCG->nRuns  = 7000; //it
        parCCG->maxDenominator = 100; //max denom
        parCCG->nSizeR1 = 6; //n restricao de cover res
    }
    else
    {
        parCCG->numberMaxConst = 16;
        parCCG->nRuns  = 50000;
        parCCG->maxDenominator = 100;
        parCCG->nSizeR1 = 12;
    }
//    printf("parameters: \n precision: %d, numberMax: %d, nRuns: %d, MaxDen: %d, sizeR1: %d\n",parCCG->precision,parCCG->numberMaxConst,parCCG->nRuns,parCCG->maxDenominator,parCCG->nSizeR1);
//    getchar();

}

