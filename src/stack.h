#ifndef STACK_H
#define STACK_H

typedef struct _Stack Stack;

/* returns an initialized object stack */
Stack *Stk_create( int capacity );
int Stk_push(Stack *p, int value);
int Stk_pop (Stack *p);
int Stk_top(Stack *p );
int Stk_isFull(Stack *p);
int Stk_isEmpty(Stack *p);
int Stk_restack(Stack *p, Stack *p_par, Stack *p_impar);
int Stk_print(Stack *p);
void Stk_free( Stack **_p );

#endif // STACK_H

