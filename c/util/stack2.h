#ifndef _STACK_H_
#define _STACK_H_

#define SIZE 100 // initial size of array
#define LF_HIGH 0.75 // upper load factor
#define LF_LOW 0.25 // lower load factor

typedef struct {
  int *x;
  int *y;
  int top; // index of top of stack
  int size; // current size of arrays
} stack;

stack *newStack(void);
void destroyStack(stack *s);
void push(stack *s, int x, int y);
int pop(stack *s, int *x, int *y);

#endif
