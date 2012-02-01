OS=linux

# Makefile for VERGINI boundary quantum billiards package.

#       Barnett 2/9/98
#	other architectures compilation 	7/26/99
#	Tidied up, failed to do fancy things, but works		99/8/5
# 	Added solve.o 99/8/16
#	Altered so expects useful_*.o and nrutil_*.o
#           to be in $HOME/lib/ 99/10/12
# 	Simplified for non-machine-specific .o files. 99/10/14
#	Changed to C++ on Linux and DEC, 99/10/28
#       Split off the bdry_def deformation program. 00/1/16
#       Reorg of header files & dependencies, 00/4/26.
#       Laptop redhat i386, 3/12/02
#	8/20/03 simplified, for redhat i386.
# 	9/8/03 Cleaned out for simple object-oriented design.
#       1/23/04 Added architecture switches for other processors (again!)
#	5/20/04 added Mac OSX
#	1/10/05 updated for Fedora Core 3 / ATLAS 3.7.8 (arch=LinuxP3)
#	5/30/05 auto linux switch for intel, or opteron w/ GOTO blas.

# architecture-specific switches......................................

# i386/linux with GNU compiler: optim flags from Jon Stott tMCimg.
# regular (redhat-supplied) gnu math libraries in glibc-devel ...
# MATHLIBS_linux = -lm -lblas -llapack -lg2c
#
# GSL + ATLAS libs, requires links by hand (3/5/06) (or lib64 version)
# sudo ln -s libatlas.so.3 /usr/lib/atlas/libatlas.so
# sudo ln -s libcblas.so.3 /usr/lib/atlas/libcblas.so
# sudo ln -s libf77blas.so.3 /usr/lib/atlas/libf77blas.so
# sudo ln -s liblapack.so.3 /usr/lib/atlas/liblapack.so
MATHLIBS_linux := -lgsl -llapack -lblas -lf77blas -latlas -lm
#
# non-GSL working lib choice:
# ATLAS libs (faster dgemm, same or slower dsyev)...
# ATLASlibinc = -I/usr/local/ATLAS/include/LinuxP3 -L/usr/local/ATLAS/lib/LinuxP3
# MATHLIBS_linux := $(ATLASlibinc) -llapack -lcblas -lf77blas -latlas -lg2c -lm
#
# add for debugging / profiling...     -ggdb -pg
# CXXFLAGS_linux = -ggdb -pg -v -Wall -pedantic
CXXFLAGS_linux = -O2 -funroll-loops -malign-double -fstrict-aliasing
CXX_linux = g++

# i386/linux with Intel compiler (on fricka)
ATLASlibinc = -I/usr/local/ATLAS/include/LinuxP3 -L/usr/local/ATLAS/lib/LinuxP3
MATHLIBS_lnx_int = $(ATLASlibinc) -llapack -lcblas -lf77blas -latlas -lg2c -lm
CXXFLAGS_lnx_int = -v -O2 -xK -ipo  -static
CXX_lnx_int = icpc

# Opteron/linux (on oliver)
MATHLIBS_lnx_amd_smp = -lgsl -llapack $(HOME)/lib/libgoto.so -lpthread -lm
CXXFLAGS_lnx_amd_smp = -O2
CXX_lnx_amd_smp = g++

# Macintosh/OSX on G4 (Apple)
MATHLIBS_osx = -lm -framework vecLib
CXXFLAGS_osx = -O2 -faltivec
CXX_osx = g++

# sparc/Solaris (Sun)
MATHLIBS_solaris = -lm -xlic_lib=sunperf
CXXFLAGS_solaris = -O2 -dalign
CXX_solaris = CC

# IRIX/MIPS (SGI)
MATHLIBS_irix = -lm -lcomplib.sgimath
CXXFLAGS_irix = -O2
CXX_irix = f77
# default rule to work under IRIX make: (default rules don't work?)
#CC = $(CXX_irix)
#CFLAGS = $(CXXFLAGS_irix)

# given you're in Linux arch, decide which options to use...
#ifeq ($(OS), linux)
#	ifeq ($(shell uname -i), i386)
#		OS := linux
#	else
#		OS := lnx_opt
#	endif
#endif

# copy flag into generic names...
MATHLIBS = $(MATHLIBS_$(OS)) # use my OS bash env variable as suffix
CXXFLAGS = $(CXXFLAGS_$(OS))
CXX = $(CXX_$(OS))
# end switches .........................................................

default: .showarch verg

all:	test low verg cf2sta bindata binmulti

VERG = verg.o matrix.o billiard.o basis.o colloc.o grid.o nrutil.o useful.o

verg:	cleancore $(VERG)
	$(CXX) -o verg $(VERG) $(CXXFLAGS) $(MATHLIBS)
	cd $(HOME)/bin; rm -f verg; ln -sf $(PWD)/verg verg

TEST = test.o billiard.o basis.o colloc.o grid.o nrutil.o useful.o

test:	cleancore $(TEST)
	$(CXX) -o test $(TEST) $(CXXFLAGS) $(MATHLIBS)

LOW = low.o matrix.o billiard.o basis.o colloc.o grid.o nrutil.o useful.o

low:	cleancore $(LOW)
	$(CXX) -o low $(LOW) $(CXXFLAGS) $(MATHLIBS)
	cd $(HOME)/bin; rm -f low; ln -sf $(PWD)/low low

CF2STA = cf2sta.o basis.o billiard.o grid.o nrutil.o useful.o

cf2sta:	cleancore $(CF2STA)
	$(CXX) -o cf2sta $(CF2STA) $(CXXFLAGS) $(MATHLIBS)
	cd $(HOME)/bin; rm -f cf2sta; ln -sf $(PWD)/cf2sta cf2sta

bindata: bindata.o
	$(CXX) $(CXXFLAGS) -o bindata bindata.o

binm: binmulti.cc
	$(CXX) $(CXXFLAGS) binmulti.cc -o binm


# object file dependency and rule...  should be rid of compilation line though
test.o:	test.cc basis.h billiard.h colloc.h grid.h useful.h nrutil.h
	$(CXX) $(CXXFLAGS) test.cc -c

cf2sta.o: cf2sta.cc basis.h billiard.h colloc.h grid.h
	$(CXX) $(CXXFLAGS) cf2sta.cc -c 

billiard.o: billiard.cc billiard.h
	$(CXX) $(CXXFLAGS) billiard.cc -c

colloc.o: colloc.cc colloc.h basis.h billiard.h
	$(CXX) $(CXXFLAGS) colloc.cc -c

basis.o: basis.cc basis.h billiard.h
	$(CXX) $(CXXFLAGS) basis.cc -c

grid.o: grid.cc grid.h basis.h billiard.h nrutil.h
	$(CXX) $(CXXFLAGS) grid.cc -c

matrix.o: matrix.cc nrutil.h useful.h cxml.h
	$(CXX) $(CXXFLAGS) matrix.cc -c

bindata.o: bindata.cc
	$(CXX) $(CXXFLAGS) bindata.cc -c

nrutil.o: nrutil.cc
	$(CXX) $(CXXFLAGS) nrutil.cc -c

useful.o: useful.cc
	$(CXX) $(CXXFLAGS) useful.cc -c

verg.o: verg.cc basis.h billiard.h colloc.h grid.h nrutil.h useful.h
	$(CXX) $(CXXFLAGS) verg.cc -c $(MATHLIBS)

# pack it up...
tar:
	cd ..; tar zcvf vergini.tar.gz vergini/*.cc vergini/*.h vergini/*.m \
	vergini/Makefile

# tell us the architecture
.showarch:
	echo $(OS) operating system

# clean up all .o and executables...

clean:
	rm -f core *.o test bindata binm verg
cleancore:
	rm -f core
