CC = g++
CFLAGS = -g -lgsl -lgslcblas -lm -funroll-loops -malign-double -fstrict-aliasing -m64

objects = count_nodal_domains.o random_percolation.o stack2.o util.o interp_matrix.o
dependencies = ../vergini/billiard.o ../vergini/useful.o ../vergini/nrutil.o util/util_verg.o


VERG = ../vergini/matrix.o ../vergini/billiard.o ../vergini/basis.o ../vergini/colloc.o ../vergini/grid.o ../vergini/nrutil.o ../vergini/useful.o util/util_verg.o
MATHLIBS = -I/usr/include/atlas -L/usr/lib64/atlas -lgsl -llapack -lblas -lf77blas -latlas -lm

count: nodal_domain_driver.o $(objects) $(dependencies)
	$(CC) $(CFLAGS) -o count nodal_domain_driver.o $(objects) $(dependencies)

perc: percolation_driver.o $(objects)
	$(CC) $(CFLAGS) -o perc percolation_driver.o $(objects)

verg_and_count:	verg_and_count.o nodal_domain_driver_no_main.o ../vergini/verg_no_main.o $(objects) $(dependencies)
	$(CC) $(CFLAGS) -o vc verg_and_count.o nodal_domain_driver_no_main.o ../vergini/verg_no_main.o $(objects) $(VERG) $(MATHLIBS)

nodal_domain_driver.o: nodal_domain_driver.c
	$(CC) $(CFLAGS) -c nodal_domain_driver.c

nodal_domain_driver_no_main.o: nodal_domain_driver_no_main.c
	$(CC) $(CFLAGS) -c nodal_domain_driver_no_main.c

percolation_driver.o: percolation_driver.c
	$(CC) $(CFLAGS) -c percolation_driver.c

verg_and_count.o: verg_and_count.c
	$(CC) $(CFLAGS) -c verg_and_count.c

count_nodal_domains.o: count_nodal_domains.c count_nodal_domains.h
	$(CC) $(CFLAGS) -c count_nodal_domains.c

random_percolation.o: random_percolation.c random_percolation.h
	$(CC) $(CFLAGS) -c random_percolation.c

stack2.o: util/stack2.c util/stack2.h
	$(CC) $(CFLAGS) -c util/stack2.c

util.o: util/util.c util/util.h
	$(CC) $(CFLAGS) -c util/util.c

util_verg.o: util/util_verg.c util/util_verg.h util.o
	$(CC) $(CFLAGS) -c util/util_verg.c

interp_matrix.o: util/interp_matrix.c util/interp_matrix.h util/util.h
	$(CC) $(CFLAGS) -c util/interp_matrix.c

.PHONY : clean
clean:
	-rm count perc vc *~ *.o