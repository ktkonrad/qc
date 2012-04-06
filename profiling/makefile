
CC=gcc
FLAGS=-lm -lgcc -lblas -lgsl

profile: hankel_profiling.c
	$(CC) $(FLAGS) hankel_profiling.c -o profile

hank103.c: hank103.f
	f2c hank103.f