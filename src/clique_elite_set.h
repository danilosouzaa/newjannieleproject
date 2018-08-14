#ifndef CLIQUE_ELITE_SET_H
#define CLIQUE_ELITE_SET_H

#include "clique.h"

typedef struct _CliqueEliteSet CliqueEliteSet;

CliqueEliteSet *clqES_create( const CliqueSet *_cs, const int numCols, const double x[], const int capacity, const double simRate, int nCols );
void clqES_free( CliqueEliteSet **_clqES );
int clqES_size( const CliqueEliteSet *clqES );
const int* clqES_get_clique_elements( const CliqueEliteSet *clqES, const int i );
const int clqES_get_clique_size( const CliqueEliteSet *clqES, const int i );
const int clqES_get_weight( const CliqueEliteSet *clqES, const int i );
void clqES_insert_inc( CliqueEliteSet *clqES, const int clqSize, const int* clqEl, const double weight, int nCols );
#endif

