#ifndef _STACK_H_
#define _STACK_H_

typedef struct node {
  int x;
  int y;
  struct node *next;
} node;

typedef struct {
  node *head;
} stack;

stack *newStack(void);

void destroyStack(stack *s);

void push(stack *s, int x, int y);

int pop(stack *s, int *x, int *y);

#endif
