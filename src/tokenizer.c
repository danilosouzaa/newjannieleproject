
/* GOAL PSP Solver II
 * Multi-Mode Resource-Constrained Multi-Project Scheduling Problem (MRCMPSP)
 *
 * Develop as part of the D.Sc. thesis of Soares, Janniele A., with collaboration
 *                                   of Santos, H.G., Toffolo, T. and Baltar, D.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "tokenizer.h"
#include "macros.h"
#include "vec_str.h"
#include "vec_char.h"
#include "vec_int.h"

struct _Tokenizer {
    VecChar *str;
    VecInt *pos;

    VecStr *tokens;
};

Tokenizer *Tok_create()
{
    Tokenizer *tok;
    ALLOCATE( tok, Tokenizer );
    tok->str = VChar_create();
    tok->pos = VInt_create();
    tok->tokens = VStr_create( STR_SIZE );

    return tok;
}

void Tok_parse( Tokenizer *tok, const char* line, const char delim )
{
    {
        char sf[2] = { delim, '\0' };
        if ( !strstr( line, sf ) ) {
            VStr_pushBack( tok->tokens, line );
            return;
        }
    }

    VChar_clear( tok->str );
    VInt_clear( tok->pos );
    VStr_resize( tok->tokens, 0 );

    const int len = strlen( line );
    /* filling string */
    {
        int i;
        for ( i=0; (i<len); ++i )
            VChar_pushBack( tok->str, line[i] );
    }
    VChar_pushBack( tok->str, '\0' );

    /* positions of tokens */
    const char *s = VChar_getPtr( tok->str );
    const char *sp;
    char stDelim[2];
    sprintf( stDelim, "%c", delim );
    int p = 0;
    while ( (sp=strstr(s,stDelim)) ) {
        p += (sp-s);
        VInt_pushBack( tok->pos, p );
        s = sp+1;
        ++p;
    }
    VInt_pushBack( tok->pos, len );

    s = VChar_getPtr( tok->str );

    /* copying contents */
    {
        int p = 0;
        int i;
        for ( i=0; (i<VInt_size( tok->pos )); ++i ) {
            char cell[STR_SIZE];
            int k=0;
            int end = VInt_get( tok->pos, i );
            for (; (p<end); ++p,++k )
                cell[k] = s[p];
            cell[k] = 0;
            VStr_pushBack( tok->tokens, cell );
            p++;
        }
    }
}

int Tok_nTokens(const Tokenizer* tok)
{
    return VStr_size( tok->tokens );
}

const char *Tok_token(const Tokenizer* tok, int i)
{
    return VStr_get( tok->tokens, i );
}

void Tok_free(Tokenizer** tok)
{
    VChar_free( &((*tok)->str) );
    VStr_free( &((*tok)->tokens) );
    VInt_free( &((*tok)->pos) );
    free( *tok );
    *tok = NULL;
}

