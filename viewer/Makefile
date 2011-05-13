# Makefile for VIEWER interactive 1d and 2d multiple array displayer.
#
# Linux Alpha version.
#
# Barnett 00/3/28

CFLAGS =

#GLUTLIBS = -lglut -lMesaGLU -lMesaGL -L/usr/X11R6/lib -lXmu -lXext -lX11 -lX
GLUTLIBS = -lglut -lGLU -lGL -L/usr/X11R6/lib -lXmu -lXext -lX11 -lXi

MATHLIBS = -lm

CC = gcc

OBJS = viewer.o nrutil.o useful.o

default:	viewer

viewer:	$(OBJS)
	$(CC) $(OBJS) -o viewer_linux $(CFLAGS) $(MATHLIBS) $(GLUTLIBS)
	cp -f viewer_linux $(HOME)/bin/

viewer.o:	viewer.c constants.h nrutil.h useful.h
	$(CC) viewer.c -c

nrutil.o:	nrutil.c nrutil.h
	$(CC) nrutil.c -c

useful.o:	useful.c useful.h
	$(CC) useful.c -c

clean:
	rm -f core *.o viewer viewer_linux

