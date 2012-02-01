CC=gcc
CFLAGS = -g -lgsl -lblas -lm

interp: interp_matrix.o
	$(CC) $(CFLAGS) -o interp interp_matrix.o

interp_matrix.o: interp_matrix.c interp_matrix.h
	$(CC) $(CFLAGS) -c interp_matrix.c