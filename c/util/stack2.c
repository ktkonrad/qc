/*
dynamically sized array implementation of a stack
data of each node is two integers x and y

Kyle Konrad
3/30/2011
*/

#include "stack2.h"
#include <stdlib.h>
#include <stdio.h>

/*
create a new stack
*/
stack *newStack() {
  stack *s = (stack *)malloc(sizeof(stack));
  s->size = SIZE;
  s->x = (int *)malloc(SIZE * sizeof(int));
  s->y = (int *)malloc(SIZE * sizeof(int));
  s->top = -1;
  return s;
}

/*
destroy a stack
*/
void destroyStack(stack *s) {
  free(s->x);
  free(s->y);
  free(s);
}
/*
push a point with x-coordinate x and y-coordinate y onto s
*/
void push(stack *s, int x, int y) {

  // double size of array if it is more full than the loadfactor LF (defined in stack2.h
  if (s->top >= (s->size) * LF_HIGH) {
    s->x = (int *)realloc(s->x, (size_t)(s->size * 2 * sizeof(int)));
    s->y = (int *)realloc(s->y, (size_t)(s->size * 2 * sizeof(int)));
    s->size *= 2;
    if (s->x == NULL || s->y == NULL) {
      fprintf(stderr, "stack2: FATAL ERROR:failed to allocate memory\n");
      exit(1);
    }
  }

  s->x[++s->top] = x;
  s->y[s->top] = y;
}

/*
pop top of stack
put values into x and y

output:
       return value: boolean indicating whether something was popped
*/
int pop(stack *s, int *x, int *y) {
  if (s->top < 0)
    return 0; // empty stack

  *x = s->x[s->top];
  *y = s->y[s->top--];

  if (s->top > SIZE && s->top <= s->size * LF_LOW) {
    s->x = (int *)realloc(s->x, s->size / 2 * sizeof(int));
    s->y = (int *)realloc(s->y, s->size / 2 * sizeof(int));
    s->size /= 2;
  }

  return 1;
}
