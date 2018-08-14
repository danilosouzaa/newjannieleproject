#include <stdio.h>
#include "stack.h"
#include "macros.h"

struct _Stack {
    int top;
    int capacity;
    int *item;
};

Stack *Stk_create( int capacity )
{

    Stack* stack;

    ALLOCATE_INI( stack, Stack );

    stack->top = -1;
    stack->capacity = capacity;

    ALLOCATE_VECTOR_INI( stack->item, int, stack->capacity );

    return stack;
}


int Stk_isFull(Stack *p)
{
    if (p->top == (p->capacity-1)) {
        //   printf ("\n\n\t\tA Stack is Full!!!");
        return 1;
    } else
        return 0;
}

int Stk_isEmpty(Stack *p)
{
    if (p->top == -1) {
        //   printf ("\n\n\t\tA Stack is Empty!!!");
        return 1;
    } else {
        //printf ("\n\n\t\tA Stack isn`t Empty!!!");
        return 0;
    }
}

int Stk_push(Stack *p, int value)
{
    p->top++;
    p->item[p->top] = value;
    return 1;
}

int Stk_pop (Stack *p)
{
    int aux;
    aux = p->item[(p->top)--];
    return aux;
}

int Stk_restack(Stack *p, Stack *p_par, Stack *p_impar)
{
    if (Stk_isEmpty(p))
        return 1;
    else {
        int aux;
        if (p->item[p->top]%2) {
            aux = Stk_pop(p);
            Stk_push(p_impar, aux);
        } else {
            aux = Stk_pop(p);
            Stk_push(p_par, aux);
        }
        Stk_restack(p, p_par, p_impar);
        return 0;
    }

}

int Stk_print(Stack *p)
{
    if (Stk_isEmpty(p)) {
        printf("-1,");
        return 1;
    } else {
        int aux;
        aux = Stk_pop(p);
        printf("%d,", aux);
        Stk_print(p);
        return 0;
    }
}

int Stk_top(Stack *p )
{

    return p->top;
}

void Stk_free( Stack **_p )
{

    Stack *p = *_p;

    free( p->item );

    free( p );

    *_p = NULL;
}

