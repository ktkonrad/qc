/*
singly-linked list implementation of a stack
data of each node is two integers x and y

Kyle Konrad
3/30/2011
*/

#include "stack.h"
#include <stdlib.h>

/*
create a new stack
*/
stack *newStack() {
  stack *s = malloc(sizeof(stack));
  s->head = NULL;
  return s;
}

/*
destroy a stack
*/
void destroyStack(stack *s) {
  free(s);
}

/*
push a point with x-coordinate x and y-coordinate y onto s
*/
void push(stack *s, int x, int y) {
  node *n = malloc(sizeof(node));
  n->next = s->head;
  n->x = x;
  n->y = y;
  s->head = n;
}

/*
pop top of stack
put values into x and y

output:
       return value: boolean indicating whether something was popped
*/
int pop(stack *s, int *x, int *y) {
  if (s->head == NULL)
    return 0;
  node *n = s->head;
  s->head = n->next;
  *x = n->x;
  *y = n->y;
  free(n);
  return 1;
}
