/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#ifndef TOKENIZER_H_INCLUDED
#define TOKENIZER_H_INCLUDED

typedef struct _Tokenizer Tokenizer;

Tokenizer *Tok_create();

void Tok_parse( Tokenizer *tok, const char *line, const char delim );

int Tok_nTokens( const Tokenizer *tok );

const char *Tok_token( const Tokenizer *tok, int i );

void Tok_free( Tokenizer **tok );

#endif

