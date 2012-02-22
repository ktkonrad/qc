CC = g++
CFLAGS = -g -lgsl -lgslcblas -lm

objects = count_nodal_domains.o random_percolation.o stack2.o util.o interp_matrix.o
dependencies = ../vergini/billiard.o ../vergini/useful.o ../vergini/nrutil.o util_verg.o

count: nodal_domain_driver.o $(objects) $(dependencies)
	$(CC) $(CFLAGS) -o count nodal_domain_driver.o $(objects) $(dependencies)

perc: percolation_driver.o $(objects)
	$(CC) $(CFLAGS) -o perc percolation_driver.o $(objects)

verg_and_count:	verg_and_count.o
	gcc -g -o vc verg_and_count.o

nodal_domain_driver.o: nodal_domain_driver.c
	$(CC) $(CFLAGS) -c nodal_domain_driver.c

percolation_driver.o: percolation_driver.c
	$(CC) $(CFLAGS) -c percolation_driver.c

verg_and_count.o: verg_and_count.c
	gcc -g -c verg_and_count.c

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
	-rm count perc $(objects) interp_matrix_k=*.dat *~