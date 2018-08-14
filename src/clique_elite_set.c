#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include "clique_elite_set.h"
#include "node_heap.h"
#include "memory.h"
#include "macros.h"

#define INFTY (INT_MAX/2)

struct _CliqueEliteSet {
    int n, cap, numCols;
    double simRate;
    const double *x;

    int **clq, *clqSize, *w, *nzs;
};

int clq_get_nz( const CliqueEliteSet *clqES, const int clqSize, const int* clqEl );
void clqES_insert( CliqueEliteSet *clqES, const int clqSize, const int* clqEl, const double weight );
int similarity(const double *x, const int n1, const int* clq1, const int n2, const int* clq2);

CliqueEliteSet *clqES_create( const CliqueSet *_cs, const int numCols, const double x[], const int capacity, const double simRate, int nCols )
{
    int i;
    CliqueEliteSet *ces = xmalloc(sizeof(CliqueEliteSet));

    ces->numCols = numCols;
    ces->x = x;
    ces->cap = capacity;
    ces->simRate = simRate;
    ces->n = 0;

    ces->w = xmalloc(capacity * sizeof(int));
    ces->nzs = xmalloc(capacity * sizeof(int));
    ces->clqSize = xmalloc(capacity * sizeof(int));
    ces->clq = xmalloc(capacity * sizeof(int*));
    for(i = 0; i < capacity; i++)
        ces->clq[i] = xmalloc(numCols * sizeof(int));

    for(i = 0; i < clq_get_number_of_cliques(_cs); i++)
        //clqES_insert(ces, clq_set_clique_size(_cs, i), clq_set_clique_elements(_cs, i), clq_set_weight(_cs, i));
        clqES_insert_inc(ces, clq_set_clique_size(_cs, i), clq_set_clique_elements(_cs, i), clq_set_weight(_cs, i), nCols);

    return ces;
}

void clqES_free( CliqueEliteSet **_clqES )
{
    CliqueEliteSet *clqES = *_clqES;
    int i;

    for(i = 0; i < clqES->cap; i++)
        free( clqES->clq[i] );
    free( clqES->clq );

    free(clqES->clqSize);
    free(clqES->nzs);
    free(clqES->w);

    free( *_clqES );
}

int clq_get_nz( const CliqueEliteSet *clqES, const int clqSize, const int* clqEl )
{
    int nNzClq = 0;
    int j;
    for ( j=0 ; (j<clqSize) ; ++j )
        if (fabs(clqES->x[clqEl[j]])>=0.001)
            nNzClq++;

    return nNzClq;
}

int clqES_size( const CliqueEliteSet *clqES )
{
    return clqES->n;
}

void clqES_insert( CliqueEliteSet *clqES, const int clqSize, const int* clqEl, const double weight )
{
    int nz = clq_get_nz(clqES, clqSize, clqEl);

    if(clqES->n < clqES->cap) {
        memcpy(clqES->clq[clqES->n], clqEl, sizeof(int) * clqSize);
        clqES->clqSize[clqES->n] = clqSize;
        clqES->w[clqES->n] = weight;
        clqES->nzs[clqES->n] = nz;
        clqES->n++;
        return;
    }

    int i;
    double maxSim = -1.0;
    int idxMaxSim = 0, minW = INT_MAX, idxMinW = 0;
    for(i = 0; i < clqES->n; i++) {
        int sim = similarity(clqES->x, clqSize, clqEl, clqES->clqSize[i], clqES->clq[i]);
        double rSim = ((double)sim) / ((double)clqES->nzs[i]);
        if(rSim > maxSim + 1E-6) {
            maxSim = rSim;
            idxMaxSim = i;
        }
        if(clqES->w[i] < minW) {
            minW = clqES->w[i];
            idxMinW = i;
        }
    }

    int idxToReplace = (maxSim + 1E-6 < clqES->simRate) ? idxMinW : idxMaxSim;

    if(weight < clqES->w[idxToReplace])
        return;
    else if(weight == clqES->w[idxToReplace]) {
        if( (nz < clqES->nzs[idxToReplace]) ||
            (nz == clqES->nzs[idxToReplace] && ((clqSize - nz) < (clqES->clqSize[idxToReplace] - clqES->nzs[idxToReplace]))) )
            return;
    }

    memcpy(clqES->clq[idxToReplace], clqEl, sizeof(int) * clqSize);
    clqES->clqSize[idxToReplace] = clqSize;
    clqES->w[idxToReplace] = weight;
    clqES->nzs[idxToReplace] = nz;
}


void clqES_insert_inc( CliqueEliteSet *clqES, const int clqSize, const int* clqEl, const double weight, int nCols )
{
    int nz = clq_get_nz(clqES, clqSize, clqEl);

    if(clqES->n < clqES->cap) {
        memcpy(clqES->clq[clqES->n], clqEl, sizeof(int) * clqSize);
        clqES->clqSize[clqES->n] = clqSize;
        clqES->w[clqES->n] = weight;
        clqES->nzs[clqES->n] = nz;
        clqES->n++;
        return;
    }


    int continci[clqES->n];
    int i;
    double maxSim = -1.0;
    int idxMaxSim = 0, minW = INT_MAX, idxMinW = 0;
    int *incElements;
    ALLOCATE_VECTOR_INI(incElements,int,nCols);
    int *elements;
    ALLOCATE_VECTOR_INI(elements,int,nCols);

    for(i = 0; i < clqES->n; i++) {
        int nElem = 0;
        continci[i] = 0;
        //printf("continci[%d] %d clqES->n %d \n", i, continci[i], clqES->n);fflush(stdout);
        for(int nc =  0 ; nc < nCols ; nc ++) {
            if(nc < clqSize) {
                const int idxA = clqEl[nc];
                // printf(" idxA %d \n", idxA); fflush(stdout);
                if(incElements[idxA]==0) {
                    incElements[idxA] = 1;
                    elements[nElem] = idxA;
                    nElem++;
                } else {
                    //      printf(" inc idxA %d \n", idxA); fflush(stdout);
                    continci[i] +=1;
                }
            }
            if(nc < clqES_get_clique_size(clqES,i) ) {
                int elem = clqES->clq[i][nc];
                //    printf(" idxB %d \n", elem); fflush(stdout);
                if(incElements[elem]==0) {
                    incElements[elem] = 1;
                    elements[nElem] = elem;
                    nElem++;
                } else {
                    //          printf(" inci idxB %d \n", elem);fflush(stdout);
                    continci[i] +=1;
                }
            }

            if((nc>clqSize) && (nc>clqES_get_clique_size(clqES,i))) break;
        }
        // printf("clqSize %d clqES->cap %d, clqES->n %d continici[%d] = %d \n",clqSize, clqES->cap,  clqES->n, i, continci[i]);
        //if(continci[i]>0) getchar();
        for(int ie =0 ; ie <nElem ; ie++) {
            incElements[elements[ie]] = 0;
            elements[ie] = 0;
        }
        int sim = continci[i];
        double rSim = ((double)sim) / ((double)clqES->nzs[i]);
        if(rSim > maxSim + 1E-6) {
            maxSim = rSim;
            idxMaxSim = i;
        }
        if(clqES->w[i] < minW) {
            minW = clqES->w[i];
            idxMinW = i;
        }
    }

    free(elements);
    free(incElements);


    int idxToReplace = (maxSim + 1E-6 < clqES->simRate) ? idxMinW : idxMaxSim;

    if(weight < clqES->w[idxToReplace])
        return;
    else if(weight == clqES->w[idxToReplace]) {
        if( (nz < clqES->nzs[idxToReplace]) ||
            (nz == clqES->nzs[idxToReplace] && ((clqSize - nz) < (clqES->clqSize[idxToReplace] - clqES->nzs[idxToReplace]))) )
            return;
    }

    memcpy(clqES->clq[idxToReplace], clqEl, sizeof(int) * clqSize);
    clqES->clqSize[idxToReplace] = clqSize;
    clqES->w[idxToReplace] = weight;
    clqES->nzs[idxToReplace] = nz;
}


char has_element( const int *v, const int n, const int key )
{
    int l = 0;
    int r = n-1;
    int m;

    while(l <= r) {
        m = (l + r) / 2;
        if(v[m] == key)
            return 1;
        if(key < v[m])
            r = m - 1;
        else
            l = m + 1;
    }

    return 0;
}

/* it returns how many elements clq1 have in common with clq2 */
int similarity(const double *x, const int n1, const int* clq1, const int n2, const int* clq2)
{
    int i, common = 0;
    for(i = 0; i < n1; i++) {
        const int element = clq1[i];
        if(x[element] >= 0.001 && has_element(clq2, n2, element))
            common++;
    }
    return common;
}

const int* clqES_get_clique_elements( const CliqueEliteSet *clqES, const int i )
{
    assert(i >= 0 && i < clqES->n);
    return clqES->clq[i];
}

const int clqES_get_clique_size( const CliqueEliteSet *clqES, const int i )
{
    assert(i >= 0 && i < clqES->n);
    return clqES->clqSize[i];
}

const int clqES_get_weight( const CliqueEliteSet *clqES, const int i )
{
    assert(i >= 0 && i < clqES->n);
    return clqES->w[i];
}
